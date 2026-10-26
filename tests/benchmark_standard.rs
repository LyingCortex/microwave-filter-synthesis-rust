//! Conformance of mfs to the cross-language filter benchmark standard.
//!
//! Three jobs, in order of importance:
//!
//! 1. **Anchor verification** — the printed Cameron matrices in
//!    `golden/twoport_anchors.json` are real published data, so feeding an
//!    anchor's *transversal* matrix and its *folded* matrix through mfs's own
//!    response engine and requiring point-wise agreement below the contract's
//!    `response_invariance` tolerance proves that engine agrees with the
//!    standard (and, transitively, with fuqtgui's `response_invariance` gate).
//! 2. **Coverage** — every case mfs's regression tests rely on must exist in the
//!    standard, so `tests/folded_response_invariance.rs` can read its specs
//!    instead of hand-copying them.
//! 3. **Drift guard** — the vendored `bench-data/` copy must equal the
//!    authoritative `filter-bench/data/` when the sibling repo is present, and
//!    otherwise must match its `MANIFEST.sha256`.

#[path = "bench_support/mod.rs"]
mod bench_support;

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use mfs::design::FilterDesign;
use mfs::freq::FrequencyGrid;
use mfs::matrix::CouplingMatrix;
use mfs::response::ResponseSolver;

use bench_support::data;
use bench_support::gates;
use bench_support::sha256;

/// Tolerance of the standard's `response_invariance` gate (1e-10).
const ANCHOR_TOL: f64 = 1e-10;

// ─── 1. anchor verification ──────────────────────────────────────────────────

/// Reads a nested row-major JSON array into a matrix, preserving entry sign.
fn matrix_from_anchor(value: &serde_json::Value) -> Vec<Vec<f64>> {
    value
        .as_array()
        .expect("anchor matrix must be an array of rows")
        .iter()
        .map(|row| {
            row.as_array()
                .expect("anchor matrix row must be an array")
                .iter()
                .map(|v| v.as_f64().expect("anchor matrix entry must be numeric"))
                .collect()
        })
        .collect()
}

/// Builds an mfs `CouplingMatrix` from a `(N+2)×(N+2)` anchor matrix.
///
/// `CouplingMatrix::new(order, data)` takes flattened **row-major** storage with
/// `order = side - 2`; the anchor matrices are stored row-major exactly as
/// printed (row/column order `S, 1..N, L`), so no transposition is involved.
/// Signs are preserved verbatim: the zeros of a coupling matrix carry the
/// coupling signs, and feeding them through `abs()` would be a different filter.
fn coupling_matrix_from_anchor(name: &str, value: &serde_json::Value) -> CouplingMatrix {
    let rows = matrix_from_anchor(value);
    let side = rows.len();
    for row in &rows {
        assert_eq!(row.len(), side, "{name}: anchor matrix must be square");
    }
    let order = side - 2;
    let flat: Vec<f64> = rows.into_iter().flatten().collect();
    CouplingMatrix::new(order, flat)
        .unwrap_or_else(|e| panic!("{name}: cannot build CouplingMatrix (order {order}): {e}"))
}

fn evaluate(matrix: &CouplingMatrix, grid: &FrequencyGrid) -> mfs::response::SParameterResponse {
    ResponseSolver
        .evaluate_normalized(matrix, grid)
        .expect("response evaluation")
}

/// The contract's `response_invariance` metric between two responses.
fn invariance_metric(
    a: &mfs::response::SParameterResponse,
    b: &mfs::response::SParameterResponse,
) -> f64 {
    gates::response_invariance_metric(a, b).expect("equal-length responses")
}

/// Book matrices are transcribed to 4 decimal places and rounded independently
/// from one another, so the *printed* transversal and folded pair cannot agree
/// below ~1e-4 in response terms: 1e-4 of matrix rounding moves the response by
/// ~1.5e-3 (`f(Ω) ~ 1e-4 × 10^Ω/λmin`), which is ~10^7 × the 1e-10 gate.
/// 2e-3 is therefore a *printing-quantum* bound, not the gate tolerance: it still
/// catches a wrong topology, a transposed/reshaped matrix, or dropped signs,
/// each of which perturbs the response by O(1).
const ANCHOR_PRINTED_PAIR_TOL: f64 = 2e-3;

/// Point-wise invariance between an anchor's printed transversal and folded
/// matrices, evaluated by mfs's response engine.
///
/// The gate that matters — `response_invariance` at 1e-10 — is exercised to the
/// full contract tolerance on *book data* by
/// [`anchor_book_matrix_passes_the_response_invariance_gate`], which folds the
/// printed transversal with mfs and compares against the printed folded matrix.
/// This test additionally uses the second printed matrix as an independent probe
/// of the response engine, at the printing-quantum bound above.
fn assert_anchor_invariance(key: &str) {
    let (transversal, folded, order) = anchor_pair(key);
    assert_eq!(transversal.order(), order, "{key}: transversal order");
    assert_eq!(folded.order(), order, "{key}: folded order");
    assert_eq!(folded.side(), order + 2, "{key}: folded side must be N+2");

    let grid = contract_grid();
    let metric = invariance_metric(&evaluate(&transversal, &grid), &evaluate(&folded, &grid));

    assert!(
        metric < ANCHOR_PRINTED_PAIR_TOL,
        "{key}: mfs response differs between the printed transversal and folded matrices by \
         max|Δ| = {metric:e}, above the 4-decimal printing quantum ({ANCHOR_PRINTED_PAIR_TOL:e}); \
         the matrices are not two forms of the same filter"
    );
}

/// Reads an anchor's two printed matrices as mfs matrices.
fn anchor_pair(key: &str) -> (CouplingMatrix, CouplingMatrix, usize) {
    let anchors = data::anchors("twoport").expect("twoport anchors");
    let entry = anchors
        .get(key)
        .unwrap_or_else(|| panic!("anchor {key} missing from golden/twoport_anchors.json"));
    let order = entry["order"].as_u64().expect("anchor order") as usize;
    (
        coupling_matrix_from_anchor(key, &entry["transversal"]),
        coupling_matrix_from_anchor(key, &entry["folded"]),
        order,
    )
}

fn contract_grid() -> FrequencyGrid {
    gates::normalized_grid(
        gates::TWOPORT_RESPONSE_BOUND,
        gates::TWOPORT_RESPONSE_POINTS,
    )
    .expect("contract grid")
}

/// `response_invariance` gate at its contract tolerance (1e-10) on book data:
/// mfs's own folding of the printed transversal matrix must preserve the
/// response exactly, and that response must be the one the printed folded matrix
/// describes (`book_matrix` at the printing quantum).
fn assert_book_matrix_and_invariance_gate(key: &str) {
    let (transversal, printed_folded, _) = anchor_pair(key);
    let grid = contract_grid();

    let design =
        mfs::transform::transform_matrix(&transversal, mfs::transform::TopologyKind::Folded)
            .unwrap_or_else(|e| {
                panic!("{key}: folding the printed transversal matrix failed: {e}")
            });
    let folded = design.matrix;

    // Gate: `response_invariance`, tol 1e-10 (harness.json gate_tols_global).
    let metric = invariance_metric(&evaluate(&transversal, &grid), &evaluate(&folded, &grid));
    assert!(
        metric < ANCHOR_TOL,
        "{key}: mfs folding of the printed transversal matrix changed the response: \
         max|Δ| = {metric:e} >= {ANCHOR_TOL:e} (contract response_invariance gate)"
    );

    // Gate: `book_matrix` (soft, tol 0.02 = coupling_abs): the folded form is a
    // plane-rotation equivalence class, so compare magnitude sets, not entries.
    let deviation = gates::coupling_set_deviation(&folded, &printed_folded)
        .unwrap_or_else(|e| panic!("{key}: coupling-set comparison failed: {e}"));
    assert!(
        deviation < 0.02,
        "{key}: mfs's folded matrix leaves the printed folded matrix's equivalence class: \
         {deviation:e} >= 0.02"
    );
}

#[test]
fn anchor_four_four_book_matrix_passes_the_invariance_gate() {
    assert_book_matrix_and_invariance_gate("FOUR_FOUR");
}

#[test]
fn anchor_six_three_book_matrix_passes_the_invariance_gate() {
    assert_book_matrix_and_invariance_gate("SIX_THREE");
}

#[test]
fn anchor_four_four_transversal_and_folded_agree_in_mfs() {
    assert_anchor_invariance("FOUR_FOUR");
}

#[test]
fn anchor_six_three_transversal_and_folded_agree_in_mfs() {
    assert_anchor_invariance("SIX_THREE");
}

#[test]
fn anchors_expose_the_contract_surface() {
    let anchors = data::anchors("twoport").expect("twoport anchors");

    // The brief's explicit invariant: FOUR_FOUR's folded matrix is 6×6 (N+2 with N=4).
    let folded = &anchors["FOUR_FOUR"]["folded"];
    assert_eq!(
        folded.as_array().expect("folded rows").len(),
        6,
        "FOUR_FOUR folded matrix must have N+2 = 6 rows"
    );

    // 7-1-2 has no transcribed book matrix: the anchor must exist (so the case's
    // anchor_id resolves) and must say so with a null matrix rather than a stub.
    assert!(
        anchors.get("SEVEN_ONE_TWO").is_some(),
        "SEVEN_ONE_TWO anchor must exist for twoport_7_1_2"
    );
    assert!(
        anchors["SEVEN_ONE_TWO"]["folded"].is_null(),
        "SEVEN_ONE_TWO folded matrix is expected to be null (not transcribed)"
    );

    // Every case carrying an anchor_id must resolve to an anchor entry.
    for case in data::cases("twoport").expect("twoport cases") {
        if let Some(anchor_id) = &case.anchor_id {
            assert!(
                anchors.get(anchor_id).is_some(),
                "{}: anchor_id {anchor_id} does not resolve",
                case.qualid
            );
        }
    }
}

#[test]
fn standard_defines_the_twoport_pipeline_and_its_gates() {
    // `stages()` yields the contract's stage set. (Ordering is not fixed by the
    // frozen `Vec<String>` signature — the data is a JSON object — so only the
    // set is asserted here.)
    let stages = data::stages("twoport").expect("twoport stages");
    assert_eq!(
        stages.into_iter().collect::<BTreeSet<_>>(),
        ["synthesis", "folding", "response"]
            .into_iter()
            .map(str::to_string)
            .collect::<BTreeSet<_>>(),
        "twoport stages are part of the contract"
    );

    for name in [
        "response_invariance",
        "folded_mask",
        "book_matrix",
        "return_loss",
    ] {
        let def = data::gate(name).unwrap_or_else(|e| panic!("gate {name}: {e}"));
        assert_eq!(def.family, "twoport", "gate {name} family");
        assert!(
            !def.requires_backend.is_empty(),
            "gate {name} must declare the stages it needs"
        );
    }

    // Tolerances the benches compare against must be the contract's, not local copies.
    let tols = data::gate_tols();
    assert_eq!(tols.get("response_invariance").copied(), Some(1e-10));
    assert_eq!(tols.get("folded_mask").copied(), Some(0.0));
    assert_eq!(tols.get("coupling_abs").copied(), Some(0.02));
    assert_eq!(data::harness("twoport_4_4").timing.warmup, 3);
    assert_eq!(data::harness("twoport_4_4").timing.repeats, 9);
}

// ─── 2. coverage ─────────────────────────────────────────────────────────────

/// Cases the folded-invariance regression test reads from the standard.
///
/// Listing them here (rather than trusting the test file) is what turns
/// "the specs moved into the standard" into something the suite enforces: if a
/// case disappears from `cases/twoport.json`, this fails by name.
const REQUIRED_TWOPORT_CASES: [&str; 6] = [
    "twoport_4_4",
    "twoport_6_3",
    "twoport_4_2_unity",
    "twoport_6_2_culdesac",
    "twoport_8_3",
    "twoport_7_1_2",
];

/// `tests/folded_response_invariance.rs` peers, read from its source so the two
/// cannot drift apart silently.
const FOLDED_INVARIANCE_TEST: &str = "tests/folded_response_invariance.rs";

#[test]
fn standard_twoport_cases_cover_the_folded_invariance_specs() {
    let test_source =
        std::fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join(FOLDED_INVARIANCE_TEST))
            .expect("read tests/folded_response_invariance.rs");

    for qualid in REQUIRED_TWOPORT_CASES {
        let case = data::case(qualid).unwrap_or_else(|e| panic!("standard case {qualid}: {e}"));

        // Well-formed spec, and consistent with the published anchor.
        let order = gates::spec_usize(&case.spec, "order").expect("spec.order");
        assert!(order > 0, "{qualid}: order must be positive");
        let rl = gates::spec_f64(&case.spec, "rl").expect("spec.rl");
        assert!(rl > 0.0, "{qualid}: rl must be positive");
        let zeros = gates::spec_f64_vec(&case.spec, "zeros").expect("spec.zeros");
        assert!(
            zeros.len() <= order,
            "{qualid}: {} transmission zeros exceed order {order}",
            zeros.len()
        );

        if let Some(anchor_id) = &case.anchor_id {
            let anchors = data::anchors("twoport").expect("anchors");
            let anchor = &anchors[anchor_id];
            assert_eq!(
                anchor["order"].as_u64().map(|v| v as usize),
                Some(order),
                "{qualid}: order disagrees with anchor {anchor_id}"
            );
            assert_eq!(
                anchor["rl_db"].as_f64(),
                Some(rl),
                "{qualid}: rl disagrees with anchor {anchor_id}"
            );
        }

        // …and the regression test actually consumes it, instead of a hand copy:
        // it must reference the qualid, and must not carry any of the numeric
        // transmission-zero literals that used to be copied by hand.
        assert!(
            test_source.contains(qualid),
            "{FOLDED_INVARIANCE_TEST} no longer references {qualid}; \
             a spec appears to have been hand-copied back into the test"
        );
    }

    // The values that used to be transcribed into the test file. Any of them
    // reappearing means the hand-copy drift is back.
    const RETIRED_HAND_COPIED_ZEROS: [&str; 13] = [
        "3.7431", "1.8051", "-1.5699", "-6.191", // twoport_4_4
        "1.0900", "1.1691", "1.5057", // twoport_6_3
        "1.3217", "1.8082", // twoport_4_2_unity
        "1.5522", "-1.5522", // twoport_6_2_culdesac
        "-1.3553", "1.1093", // twoport_8_3
    ];
    for literal in RETIRED_HAND_COPIED_ZEROS {
        assert!(
            !test_source.contains(literal),
            "{FOLDED_INVARIANCE_TEST} contains the hand-copied transmission zero {literal}; \
             specs must come from the standard (bench-data/cases/twoport.json)"
        );
    }
    // 7-1-2's single zero ("1.2576") is covered by the same rule; checked with a
    // separate literal because it is short enough to appear inside other numbers.
    assert!(
        !test_source.contains("1.2576"),
        "{FOLDED_INVARIANCE_TEST} contains the hand-copied 7-1-2 transmission zero"
    );
}

#[test]
fn standard_specs_match_the_annotated_anchor_values() {
    // Sized to the contract: 13 twoport cases (Cameron book + the 7-1-2 case
    // mfs contributed + the 1.4.0 CST official cases), each with a resolving
    // anchor_id.
    let cases = data::cases("twoport").expect("twoport cases");
    assert_eq!(
        cases.len(),
        13,
        "twoport case count (contract 1.4.0 adds the CST official cases)"
    );
    for case in &cases {
        assert!(
            case.anchor_id.is_some(),
            "{}: every twoport case must point at a publication anchor",
            case.qualid
        );
        assert!(
            !case.provenance.is_null(),
            "{}: provenance is required",
            case.qualid
        );
    }
}

// ─── 3. drift guard ──────────────────────────────────────────────────────────

/// Walks `root` and returns every file's path relative to it, with `/` separators.
fn relative_files(root: &Path) -> Vec<PathBuf> {
    fn walk(dir: &Path, base: &Path, out: &mut Vec<PathBuf>) {
        let entries =
            std::fs::read_dir(dir).unwrap_or_else(|e| panic!("read_dir {}: {e}", dir.display()));
        for entry in entries {
            let path = entry.expect("dir entry").path();
            if path.is_dir() {
                walk(&path, base, out);
            } else {
                out.push(
                    path.strip_prefix(base)
                        .expect("path under base")
                        .to_path_buf(),
                );
            }
        }
    }
    let mut out = Vec::new();
    walk(root, root, &mut out);
    out.sort();
    out
}

fn authoritative_data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("filter-bench")
        .join("data")
}

/// Tier 1: the sibling standard repository is present — compare byte-for-byte.
#[test]
fn vendored_data_matches_the_authoritative_source() {
    let vendored = data::data_dir();
    let source = authoritative_data_dir();

    if !source.join("harness.json").is_file() {
        // Tier 2 (below) is the guard for solo clones / CI.
        eprintln!(
            "skipping byte-for-byte comparison: {} not present (solo clone / CI)",
            source.display()
        );
        return;
    }

    let mut compared = 0usize;
    for rel in relative_files(&source) {
        let source_file = source.join(&rel);
        let vendored_file = vendored.join(&rel);
        assert!(
            vendored_file.is_file(),
            "vendored copy is missing {} (run: python filter-bench/tools/vendor.py --into mfs/bench-data)",
            rel.display()
        );
        let lhs = std::fs::read(&source_file).expect("read source");
        let rhs = std::fs::read(&vendored_file).expect("read vendored");
        assert!(
            lhs == rhs,
            "vendored {} differs from the authoritative copy; \
             re-run: python filter-bench/tools/vendor.py --into mfs/bench-data",
            rel.display()
        );
        compared += 1;
    }
    assert!(
        compared >= 8,
        "expected the standard's data tree to have at least 8 files, compared {compared}"
    );

    // Nothing may be vendored that the source does not have.
    let source_set: BTreeSet<PathBuf> = relative_files(&source).into_iter().collect();
    for rel in relative_files(&vendored) {
        if rel == Path::new("VENDORED.md") || rel == Path::new("MANIFEST.sha256") {
            continue;
        }
        assert!(
            source_set.contains(&rel),
            "vendored file {rel:?} is not part of the authoritative data tree (stale copy?)"
        );
    }
}

/// Tier 2: standalone verification against `MANIFEST.sha256` (runs everywhere).
#[test]
fn vendored_data_matches_its_manifest() {
    let vendored = data::data_dir();
    let manifest_text = std::fs::read_to_string(vendored.join("MANIFEST.sha256"))
        .expect("bench-data/MANIFEST.sha256");
    let manifest = sha256::parse_manifest(&manifest_text);
    assert!(
        !manifest.is_empty(),
        "MANIFEST.sha256 is empty; the drift guard would be vacuous"
    );

    for (rel, expected) in &manifest {
        let path = vendored.join(rel);
        let actual = sha256::digest_file(&path).unwrap_or_else(|e| panic!("hash {rel}: {e}"));
        assert_eq!(
            &actual, expected,
            "vendored {rel} does not match MANIFEST.sha256 (local edit? re-run vendor.py)"
        );
    }

    for rel in relative_files(&vendored) {
        if rel == Path::new("VENDORED.md") || rel == Path::new("MANIFEST.sha256") {
            continue;
        }
        let key = rel.to_string_lossy().replace('\\', "/");
        assert!(
            manifest.contains_key(&key),
            "vendored file {key} is absent from MANIFEST.sha256 (stale copy?)"
        );
    }
}

/// The loader itself must be reading the standard, not a copy of its numbers.
#[test]
fn loader_reads_the_vendored_standard() {
    assert_eq!(data::check_contract().expect("contract"), "1.4.0");

    let cases = data::cases("multiplexer").expect("multiplexer cases");
    assert_eq!(
        cases.len(),
        10,
        "multiplexer case count (contract 1.4.0 adds the CST duplexers + contiguous diplexer)"
    );
    let ids: BTreeSet<&str> = cases.iter().map(|c| c.id.as_str()).collect();
    assert!(
        ids.contains("diplexer"),
        "synthetic diplexer case must exist"
    );
    assert!(
        ids.contains("diplexer_type2"),
        "type-2 diplexer case must exist"
    );
    assert_eq!(
        data::case("multiplexer_diplexer").expect("diplexer").tags,
        vec!["synthetic", "type1", "n8"]
    );

    // Unknown input must be an error, never a silent default.
    assert!(data::cases("nope").is_err());
    assert!(data::case("twoport_nope").is_err());
}

// ─── 4. support-module self-checks ───────────────────────────────────────────
//
// `tests/bench_support/` is shared by two bench targets and two test targets via
// `#[path]`, so it is not a standalone test target. Its own checks live here,
// where they are enforced by `cargo test --test benchmark_standard`.

#[test]
fn sha256_matches_published_nist_vectors() {
    // FIPS 180-4 examples. If this routine ever degrades into a placeholder, the
    // drift guard below would silently pass — so it is pinned to known digests.
    assert_eq!(
        sha256::digest(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
    assert_eq!(
        sha256::digest(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
    assert_eq!(
        sha256::digest(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
        "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
    );
}

/// A >2-block message hashed in odd-sized chunks must equal the one-shot digest
/// (exercises the streaming buffer path, which the manifest check depends on).
#[test]
fn sha256_streaming_matches_one_shot() {
    let data: Vec<u8> = (0..1000u32).map(|i| (i % 251) as u8).collect();
    let mut hasher = sha256::Sha256::new();
    for chunk in data.chunks(37) {
        hasher.update(chunk);
    }
    assert_eq!(sha256::hex(&hasher.finish()), sha256::digest(&data));
}

/// `folded_mask` gate (hard, tolerance 0.0): every matrix mfs folds must sit
/// entirely inside **one** complete Cameron folding convention - the contract
/// criterion `min(violation(canonical), violation(mirrored))` (contract 1.1.0;
/// `data/gates.json:gates.folded_mask`). Cameron's folded form is not unique
/// (§8.3): mfs's `to_folded` emits the mirrored variant for even orders, the
/// book's printed folded matrices the canonical one; satisfying either complete
/// convention passes.
///
/// Cases whose spec has no integer `order` (the 1.4.0 CST `cst_multiband`, an
/// `order: null` Arbitrary-response design) are L3 spec-only: the contract
/// mandates a capability-SKIP, not a FAIL, so they are skipped here.
#[test]
fn folded_mask_holds_for_every_standard_twoport_case() {
    for case in data::cases("twoport").expect("twoport cases") {
        let order = match gates::spec_usize(&case.spec, "order") {
            Ok(order) => order,
            Err(_) => {
                eprintln!(
                    "folded_mask: {} has no integer order (L3 spec-only) -> capability-SKIP",
                    case.qualid
                );
                continue;
            }
        };
        let rl = gates::spec_f64(&case.spec, "rl").expect("rl");
        let zeros = gates::spec_f64_vec(&case.spec, "zeros").expect("zeros");
        let design = mfs::design::FilterDesign::prototype(order, rl)
            .zeros(zeros)
            .synthesize()
            .unwrap_or_else(|e| panic!("{}: synthesis failed: {e}", case.qualid));
        let folded = design
            .to_folded()
            .unwrap_or_else(|e| panic!("{}: folding failed: {e}", case.qualid));
        assert_eq!(
            gates::folded_mask_violations(&folded),
            0.0,
            "{}: folded matrix has nonzero entries outside the contract's folded mask \
             (canonical ∪ mirrored)",
            case.qualid
        );
    }
}

/// The mask predicate must actually reject something - otherwise the check
/// above would be vacuous. `(2, N+1)` on a 4th-order filter is outside both
/// conventions, while the reduction pivot mfs emits at even order (`(2, 4)`,
/// index sum N+2) is accepted through the mirror.
#[test]
fn folded_mask_gate_rejects_a_far_corner() {
    use mfs::matrix::CouplingMatrixBuilder;

    let dirty = CouplingMatrixBuilder::new(4)
        .expect("order")
        .set(2, 5, 0.25)
        .expect("far corner")
        .allow_asymmetric(true)
        .build()
        .expect("build");
    assert!((gates::folded_mask_violations(&dirty) - 0.25).abs() < 1e-15);

    let pivot = CouplingMatrixBuilder::new(4)
        .expect("order")
        .set(2, 4, 0.0342)
        .expect("pivot")
        .allow_asymmetric(true)
        .build()
        .expect("build");
    assert_eq!(
        gates::folded_mask_violations(&pivot),
        0.0,
        "the even-order pivot mfs emits must be inside the union mask"
    );
}

/// Cross-implementation evidence: mfs's *synthesized and folded* matrix must
/// reproduce the book's **printed** matrix - two independent derivations of one
/// filter, so agreement at the printing quantum pins the chain
/// (spec -> synthesis -> folding -> response engine) against the publication.
///
/// Two conventions are absorbed, both measured rather than assumed:
///
/// * **Lattice / zero-sign**: the standard's `zeros` are signed Ω values; mfs
///   and the book resolve the source-vs-load zero assignment oppositely, which
///   reflects the response in Ω (`S(-Ω)`). The reflected comparison is used, and
///   it is what makes the agreement ~1e-4 instead of ~1e-1 (4-4: 3.6e-4 vs 9.1e-2;
///   6-3 is symmetric, so both readings agree at 9.2e-4).
/// * **Printing**: the book is transcribed to 4 decimals, hence
///   `ANCHOR_PRINTED_PAIR_TOL`.
///
/// Element-wise equality is not available at any tolerance: the folded form is a
/// plane-rotation equivalence class and coupling matrices are defined up to a
/// global sign, so the structural check compares the sorted main-chain and
/// cross-coupling magnitude sets (`coupling_abs` tol 0.02).
#[test]
fn mfs_synthesized_folded_matrix_reproduces_the_printed_book_matrix() {
    for (qualid, anchor_id) in [("twoport_4_4", "FOUR_FOUR"), ("twoport_6_3", "SIX_THREE")] {
        let case = data::case(qualid).expect("case");
        let order = gates::spec_usize(&case.spec, "order").expect("order");
        let (_, printed_folded, anchor_order) = anchor_pair(anchor_id);
        assert_eq!(
            anchor_order, order,
            "{qualid}: case and anchor order must agree"
        );

        let design = mfs::design::FilterDesign::prototype(
            order,
            gates::spec_f64(&case.spec, "rl").expect("rl"),
        )
        .zeros(gates::spec_f64_vec(&case.spec, "zeros").expect("zeros"))
        .synthesize()
        .expect("synthesis");
        let folded = design.to_folded().expect("folding");
        assert_eq!(
            gates::folded_mask_violations(&folded),
            0.0,
            "{qualid}: folded matrix must be inside the contract's folded mask"
        );

        let deviation = gates::coupling_set_deviation(&folded, &printed_folded).expect("sets");
        assert!(
            deviation < 0.02,
            "{qualid}: folded matrix leaves the printed matrix's equivalence class: {deviation:e}"
        );

        let grid = contract_grid();
        let ours = evaluate(&folded, &grid);
        let printed = evaluate(&printed_folded, &grid);
        let metric = invariance_metric(&ours, &printed)
            .min(invariance_metric(&ours, &reflected_omega(&printed)));
        assert!(
            metric < ANCHOR_PRINTED_PAIR_TOL,
            "{qualid}: mfs's folded matrix does not reproduce the printed matrix's response: \
             max|Δ| = {metric:e} >= {ANCHOR_PRINTED_PAIR_TOL:e} (4-decimal printing quantum)"
        );
    }
}

/// N=8 regression entry for the documented "order 2–8 fixture band".
///
/// `twoport_8_3` (Cameron 2018 Sec 9.6 / Fig 9.22a, p.347) is the standard's
/// only order-8 two-port case, and its anchor `EIGHT_THREE_FOLDED` carries
/// the printed **folded** matrix only (no transversal transcription), so it
/// cannot ride on [`anchor_pair`]. This test therefore goes straight from
/// the case spec through synthesis and folding, and pins the result against
/// the printed folded matrix at the printing quantum — the same two
/// conventions (zero-sign reflection, 4-decimal printing) are absorbed as in
/// [`mfs_synthesized_folded_matrix_reproduces_the_printed_book_matrix`].
#[test]
fn mfs_synthesized_folded_order8_reproduces_the_printed_book_matrix() {
    let qualid = "twoport_8_3";
    let case = data::case(qualid).expect("case");
    let order = gates::spec_usize(&case.spec, "order").expect("order");
    assert_eq!(order, 8, "{qualid}: this regression entry pins the N=8 band edge");

    let anchors = data::anchors("twoport").expect("twoport anchors");
    let printed_folded = coupling_matrix_from_anchor(
        "EIGHT_THREE_FOLDED",
        &anchors["EIGHT_THREE_FOLDED"]["folded"],
    );
    assert_eq!(printed_folded.order(), order, "anchor order must match the case");

    let design = mfs::design::FilterDesign::prototype(
        order,
        gates::spec_f64(&case.spec, "rl").expect("rl"),
    )
    .zeros(gates::spec_f64_vec(&case.spec, "zeros").expect("zeros"))
    .synthesize()
    .expect("synthesis");
    let folded = design.to_folded().expect("folding");
    assert_eq!(
        gates::folded_mask_violations(&folded),
        0.0,
        "{qualid}: folded matrix must be inside the contract's folded mask"
    );

    let deviation = gates::coupling_set_deviation(&folded, &printed_folded).expect("sets");
    assert!(
        deviation < 0.02,
        "{qualid}: folded matrix leaves the printed matrix's equivalence class: {deviation:e}"
    );

    let grid = contract_grid();
    let ours = evaluate(&folded, &grid);
    let printed = evaluate(&printed_folded, &grid);
    let metric = invariance_metric(&ours, &printed)
        .min(invariance_metric(&ours, &reflected_omega(&printed)));
    assert!(
        metric < ANCHOR_PRINTED_PAIR_TOL,
        "{qualid}: mfs's folded matrix does not reproduce the printed matrix's response: \
         max|Δ| = {metric:e} >= {ANCHOR_PRINTED_PAIR_TOL:e} (4-decimal printing quantum)"
    );
}

/// The same response read in the opposite Ω direction (`ω -> -ω`), which is how
/// the two sign conventions for the lattice transmission zeros differ.
fn reflected_omega(
    response: &mfs::response::SParameterResponse,
) -> mfs::response::SParameterResponse {
    let mut reflected = response.clone();
    reflected.samples.reverse();
    for sample in &mut reflected.samples {
        sample.normalized_omega = -sample.normalized_omega;
    }
    reflected
}

#[test]
fn response_invariance_metric_is_zero_for_identical_responses() {
    let (transversal, _, _) = anchor_pair("FOUR_FOUR");
    let response = evaluate(&transversal, &contract_grid());
    assert_eq!(
        gates::response_invariance_metric(&response, &response).expect("metric"),
        0.0
    );
}

/// The response-invariance metric must be sensitive to an actual difference —
/// otherwise the anchor checks above would be vacuous. The printed transversal
/// and folded matrices of a book anchor are exactly such a difference: two
/// transcriptions of one filter, differing by the 4-decimal printing quantum.
#[test]
fn response_invariance_metric_detects_a_real_difference() {
    let grid = contract_grid();
    let (transversal, printed_folded, _) = anchor_pair("FOUR_FOUR");
    let metric = gates::response_invariance_metric(
        &evaluate(&transversal, &grid),
        &evaluate(&printed_folded, &grid),
    )
    .expect("metric");
    assert!(
        metric > 1e-6,
        "metric failed to notice the printed transversal/folded difference: {metric:e}"
    );
    assert!(
        metric < ANCHOR_PRINTED_PAIR_TOL,
        "metric is implausibly large for two transcriptions of one filter: {metric:e}"
    );
}

// ─── 5. cross-backend equivalence (rust-native vs pinned python reference) ──

/// The reference's eigenvalues are real to within 1e-9 (folded resonator
/// subblock is real symmetric). The mfs metric is bounded by the
/// cross_backend_equivalence gate tolerance 1e-6, but the *measured* deviation
/// between the two backends is at machine precision (≤ ~3e-12 across all 12
/// synthesizable twoport cases, verified by the rust-native bench). The conformance
/// test pins the gate at its contract tolerance (1e-6) so a regression that
/// breaks the eigenvalue match is caught, while the 1e-9 realness bound pins
/// the reference against being silently replaced by a complex-valued matrix.
///
/// The reference is the python backend's (fuqtgui `domain.*`) folded
/// resonator-subblock eigenvalues, generated 2026-10-07 and pinned in
/// `bench-data/golden/cross_backend_reference.json`. The two backends use
/// *different* synthesis and folding algorithms; agreement at machine precision
/// confirms both produce the same filter.
#[test]
fn cross_backend_equivalence_matches_python_reference() {
    let xref = data::CrossBackendReference::load().expect("load cross-backend reference");
    let cases = data::cases("twoport").expect("twoport cases");
    assert!(
        !xref.cases.is_empty(),
        "cross_backend_reference.json has no cases"
    );
    // Every *synthesizable* twoport case (one with an integer `order`) must
    // have a reference entry. `order: null` cases are L3 spec-only - the
    // reference covers exactly the synthesizable set (contract: the file
    // requires coverage of every synthesizable twoport case, no more).
    for case in &cases {
        // `order: null` (L3 spec-only) cases are not in the reference set.
        if gates::spec_usize(&case.spec, "order").is_err() {
            continue;
        }
        assert!(
            xref.cases.contains_key(case.qualid.as_str()),
            "case {} missing from cross_backend_reference.json",
            case.qualid
        );
    }
    // The contract's gate tolerance (1e-6).
    let tol = data::gate_tols()
        .get("cross_backend_equivalence")
        .copied()
        .unwrap_or(1e-6);
    // Verify each case: synthesize with mfs, fold, compare eigenvalues.
    let mut worst = 0.0_f64;
    let mut checked = 0;
    for case in &cases {
        // `order: null` (L3 spec-only) cases have no reference entry.
        if gates::spec_usize(&case.spec, "order").is_err() {
            eprintln!(
                "cross_backend_equivalence: {} has no integer order (L3 spec-only) -> capability-SKIP",
                case.qualid
            );
            continue;
        }
        let entry = &xref.cases[case.qualid.as_str()];
        let order = gates::spec_usize(&case.spec, "order").expect("order");
        assert_eq!(
            entry.order, order,
            "{}: reference order {} != case order {}",
            case.qualid,
            entry.order,
            order
        );
        // Reference eigenvalues must be real (|im| < 1e-9).
        for e in &entry.eig {
            assert!(
                e.im.abs() < 1e-9,
                "{}: non-real reference eigenvalue im={:.3e}",
                case.qualid,
                e.im
            );
        }
        let design = FilterDesign::prototype(order, gates::spec_f64(&case.spec, "rl").expect("rl"))
            .zeros(gates::spec_f64_vec(&case.spec, "zeros").expect("zeros"))
            .synthesize()
            .expect("synthesis");
        let folded = design.to_folded().expect("folding");
        let dev = gates::cross_backend_equivalence_metric(&folded, entry).expect("eigenvalue metric");
        assert!(
            dev <= tol,
            "{}: mfs eigenvalue deviation {:.3e} > contract tolerance {:.1e}",
            case.qualid,
            dev,
            tol
        );
        worst = worst.max(dev);
        checked += 1;
    }
    eprintln!(
        "cross_backend_equivalence: worst={worst:.3e} over {} twoport cases (tol {:.1e})",
        checked,
        tol
    );
}
