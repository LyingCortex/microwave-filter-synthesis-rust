//! Loader for the vendored benchmark standard (`bench-data/`) plus the shared
//! Report / StageTiming / GateResult types of the cross-language contract.
//!
//! Shapes mirror the frozen `filterbench` Rust crate API. Nothing here is
//! mfs-specific; the mfs-facing helpers live in [`super::gates`].

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Loader error type (crate-mirroring alias).
pub type Error = Box<dyn std::error::Error + Send + Sync>;
/// Loader result type.
pub type Result<T, E = Error> = std::result::Result<T, E>;

/// Directory (relative to `CARGO_MANIFEST_DIR`) holding the vendored standard.
pub const DATA_DIR_REL: &str = "bench-data";

/// Absolute path of the vendored standard data.
pub fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(DATA_DIR_REL)
}

/// Local platform string, e.g. `linux-x86_64`.
pub fn platform() -> String {
    format!("{}-{}", std::env::consts::OS, std::env::consts::ARCH)
}

// ─── Case ────────────────────────────────────────────────────────────────────

/// One benchmark case of a family. Global id is `qualid == "<family>_<id>"`.
#[derive(Debug, Clone)]
pub struct Case {
    /// Family-local id, e.g. `4_4`.
    pub id: String,
    /// Family, e.g. `twoport`.
    pub family: String,
    /// Global id, e.g. `twoport_4_4`.
    pub qualid: String,
    /// Provenance record (source document / note).
    pub provenance: serde_json::Value,
    /// Capability tags, e.g. `type1`, `synthetic`, `contiguous`.
    pub tags: Vec<String>,
    /// Family-specific algorithm input (data, not interpretation).
    pub spec: serde_json::Value,
    /// Physical invariants the result must satisfy.
    pub invariants: serde_json::Value,
    /// Key into `golden/<family>_anchors.json`, when the case has a publication anchor.
    pub anchor_id: Option<String>,
}

impl Case {
    /// True when the case carries `tag`.
    pub fn has_tag(&self, tag: &str) -> bool {
        self.tags.iter().any(|t| t == tag)
    }
}

#[derive(Deserialize)]
struct CasesDoc {
    schema_version: String,
    cases: Vec<CaseRaw>,
}

#[derive(Deserialize)]
struct CaseRaw {
    id: String,
    #[serde(default)]
    provenance: serde_json::Value,
    #[serde(default)]
    tags: Vec<String>,
    spec: serde_json::Value,
    #[serde(default)]
    invariants: serde_json::Value,
    #[serde(default)]
    anchor_id: Option<String>,
}

fn read_json(path: &Path) -> Result<serde_json::Value> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    Ok(serde_json::from_str(&text)
        .map_err(|e| format!("invalid JSON in {}: {e}", path.display()))?)
}

/// All cases of `family`, in file order.
pub fn cases(family: &str) -> Result<Vec<Case>> {
    let path = data_dir().join("cases").join(format!("{family}.json"));
    let doc: CasesDoc = serde_json::from_value(read_json(&path)?)
        .map_err(|e| format!("invalid case document {}: {e}", path.display()))?;
    if doc.schema_version != "1.0" {
        return Err(format!("{path:?}: unexpected schema_version {}", doc.schema_version).into());
    }
    Ok(doc
        .cases
        .into_iter()
        .map(|raw| Case {
            qualid: format!("{family}_{}", raw.id),
            id: raw.id,
            family: family.to_string(),
            provenance: raw.provenance,
            tags: raw.tags,
            spec: raw.spec,
            invariants: raw.invariants,
            anchor_id: raw.anchor_id,
        })
        .collect())
}

/// Single case by global id `<family>_<id>`.
pub fn case(qualid: &str) -> Result<Case> {
    let (family, id) = qualid
        .split_once('_')
        .ok_or_else(|| format!("qualid {qualid:?} must be <family>_<id>"))?;
    cases(family)?
        .into_iter()
        .find(|c| c.id == id)
        .ok_or_else(|| format!("unknown case {qualid:?}").into())
}

/// Golden anchors of `family` (`golden/<family>_anchors.json`), raw JSON.
pub fn anchors(family: &str) -> Result<serde_json::Value> {
    read_json(
        &data_dir()
            .join("golden")
            .join(format!("{family}_anchors.json")),
    )
}

// ─── Harness (timing / iteration policy) ─────────────────────────────────────

/// Timing policy of one case.
#[derive(Debug, Clone, PartialEq)]
pub struct Timing {
    /// Untimed warm-up calls before measuring.
    pub warmup: usize,
    /// Measured calls; the reported statistic is taken over these.
    pub repeats: usize,
    /// Statistic name, e.g. `median`.
    pub stat: String,
}

/// Merged harness policy for one case (`defaults` <- per-case overrides).
#[derive(Debug, Clone, PartialEq)]
pub struct Harness {
    /// Timing policy.
    pub timing: Timing,
    /// Iteration bound, when the case sets one.
    pub max_iter: Option<usize>,
    /// Gate tolerances: global table overlaid with defaults and per-case entries.
    pub gate_tols: HashMap<String, f64>,
    /// Free-form note (e.g. "N=300; use --stages poly").
    pub note: Option<String>,
}

#[derive(Deserialize)]
struct HarnessDoc {
    contract_version: String,
    defaults: HarnessEntryRaw,
    #[serde(default)]
    gate_tols_global: BTreeMap<String, f64>,
    #[serde(default)]
    cases: BTreeMap<String, HarnessEntryRaw>,
}

#[derive(Deserialize, Default, Clone)]
struct HarnessEntryRaw {
    #[serde(default)]
    timing: Option<TimingRaw>,
    #[serde(default)]
    max_iter: Option<usize>,
    #[serde(default)]
    gate_tols: BTreeMap<String, f64>,
    #[serde(default)]
    note: Option<String>,
}

#[derive(Deserialize, Clone)]
struct TimingRaw {
    warmup: usize,
    repeats: usize,
    stat: String,
}

fn harness_doc() -> Result<HarnessDoc> {
    let path = data_dir().join("harness.json");
    serde_json::from_value(read_json(&path)?)
        .map_err(|e| format!("invalid harness document {}: {e}", path.display()).into())
}

impl From<TimingRaw> for Timing {
    fn from(raw: TimingRaw) -> Self {
        Self {
            warmup: raw.warmup,
            repeats: raw.repeats,
            stat: raw.stat,
        }
    }
}

/// Harness policy for `qualid`; unknown cases fall back to `defaults`.
pub fn harness(qualid: &str) -> Harness {
    let doc = harness_doc().unwrap_or(HarnessDoc {
        contract_version: String::new(),
        defaults: HarnessEntryRaw::default(),
        gate_tols_global: BTreeMap::new(),
        cases: BTreeMap::new(),
    });
    let per = doc.cases.get(qualid).cloned().unwrap_or_default();

    let mut gate_tols: HashMap<String, f64> = doc.gate_tols_global.into_iter().collect();
    gate_tols.extend(doc.defaults.gate_tols.clone());
    gate_tols.extend(per.gate_tols.clone());

    Harness {
        timing: per
            .timing
            .or(doc.defaults.timing)
            .map(Timing::from)
            .unwrap_or(Timing {
                warmup: 3,
                repeats: 9,
                stat: "median".to_string(),
            }),
        max_iter: per.max_iter.or(doc.defaults.max_iter),
        gate_tols,
        note: per.note.or(doc.defaults.note),
    }
}

/// Global gate tolerances (`harness.json:gate_tols_global`).
pub fn gate_tols() -> HashMap<String, f64> {
    harness_doc()
        .map(|d| d.gate_tols_global.into_iter().collect())
        .unwrap_or_default()
}

// ─── Gates ───────────────────────────────────────────────────────────────────

/// Canonical definition of one gate (`gates.json:gates`).
#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct GateDef {
    /// Family the gate belongs to (`*` = all families).
    pub family: String,
    /// What the gate verifies, in words.
    pub semantics: String,
    /// `hard` (failure blocks timing) or `soft`.
    pub kind: String,
    /// Metric definition.
    pub metric: String,
    /// Stages the gate needs; a backend lacking them records `SKIP`.
    pub requires_backend: Vec<String>,
    /// Key into the tolerance tables, when the gate has one.
    #[serde(default)]
    pub tol_key: Option<String>,
}

#[derive(Deserialize)]
struct GatesDoc {
    contract_version: String,
    schema_version: String,
    stages: BTreeMap<String, BTreeMap<String, String>>,
    gates: BTreeMap<String, GateDef>,
}

fn gates_doc() -> Result<GatesDoc> {
    let path = data_dir().join("gates.json");
    serde_json::from_value(read_json(&path)?)
        .map_err(|e| format!("invalid gates document {}: {e}", path.display()).into())
}

/// Canonical stage names of `family`, in contract order.
pub fn stages(family: &str) -> Result<Vec<String>> {
    let doc = gates_doc()?;
    doc.stages
        .get(family)
        .map(|s| s.keys().cloned().collect())
        .ok_or_else(|| format!("unknown family {family:?}").into())
}

/// Canonical gate definition by name.
pub fn gate(name: &str) -> Result<GateDef> {
    let doc = gates_doc()?;
    doc.gates
        .get(name)
        .cloned()
        .ok_or_else(|| format!("unknown gate {name:?}").into())
}

/// All canonical gate definitions.
pub fn gates() -> Result<BTreeMap<String, GateDef>> {
    Ok(gates_doc()?.gates)
}

/// Validates cross-file contract/schema versions; returns the contract version.
pub fn check_contract() -> Result<String> {
    let h = harness_doc()?;
    let g = gates_doc()?;
    if h.contract_version != g.contract_version {
        return Err(format!(
            "contract_version drift: harness={} gates={}",
            h.contract_version, g.contract_version
        )
        .into());
    }
    if g.schema_version != "1.0" || h.contract_version.is_empty() {
        return Err(format!(
            "unexpected schema/contract versions: {}/{}",
            g.schema_version, h.contract_version
        )
        .into());
    }
    Ok(g.contract_version)
}

// ─── Timing ──────────────────────────────────────────────────────────────────

/// Result of timing one stage (report-schema shape).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct StageTiming {
    /// Median wall time in ms.
    pub median_ms: f64,
    /// 95th percentile wall time in ms.
    pub p95_ms: Option<f64>,
    /// Fastest sample in ms.
    pub min_ms: Option<f64>,
    /// Number of measured samples.
    pub n: Option<usize>,
}

/// Times `f` after `warmup` untimed calls and `repeats` measured calls.
///
/// Mirrors the standard's reference implementation (`filterbench.time_call`):
/// samples sorted, `p95` at index `round(0.95 * (n - 1))`, median over all samples.
pub fn time_call<F: FnMut()>(mut f: F, warmup: usize, repeats: usize) -> StageTiming {
    for _ in 0..warmup {
        f();
    }
    let mut samples: Vec<f64> = Vec::with_capacity(repeats.max(1));
    for _ in 0..repeats.max(1) {
        let start = Instant::now();
        f();
        samples.push(start.elapsed().as_secs_f64() * 1e3);
    }
    samples.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = samples.len();
    let p95_index = ((0.95 * (n as f64 - 1.0)).round() as usize).min(n - 1);
    StageTiming {
        median_ms: median(&samples),
        p95_ms: Some(samples[p95_index]),
        min_ms: Some(samples[0]),
        n: Some(n),
    }
}

fn median(sorted: &[f64]) -> f64 {
    match sorted.len() {
        0 => f64::NAN,
        n if n % 2 == 1 => sorted[n / 2],
        n => 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]),
    }
}

// ─── Report ──────────────────────────────────────────────────────────────────

/// Outcome of one gate on one case.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GateResult {
    /// `PASS` / `FAIL` / `SKIP`.
    pub status: String,
    /// Observed metric, when applicable.
    pub metric: Option<f64>,
    /// Tolerance the metric was compared against, when applicable.
    pub tol: Option<f64>,
    /// Human-readable explanation (required for `SKIP`).
    pub detail: Option<String>,
}

impl GateResult {
    /// Passing gate with an observed metric.
    pub fn pass(metric: f64, tol: Option<f64>) -> Self {
        Self {
            status: "PASS".to_string(),
            metric: Some(metric),
            tol,
            detail: None,
        }
    }

    /// Failing gate with an observed metric.
    pub fn fail(metric: f64, tol: Option<f64>) -> Self {
        Self {
            status: "FAIL".to_string(),
            metric: Some(metric),
            tol,
            detail: None,
        }
    }

    /// Skipped gate; `detail` must say why the backend cannot evaluate it.
    pub fn skip(detail: impl Into<String>) -> Self {
        Self {
            status: "SKIP".to_string(),
            metric: None,
            tol: None,
            detail: Some(detail.into()),
        }
    }

    /// Fully explicit constructor.
    pub fn new(
        status: impl Into<String>,
        metric: Option<f64>,
        tol: Option<f64>,
        detail: Option<String>,
    ) -> Self {
        Self {
            status: status.into(),
            metric,
            tol,
            detail,
        }
    }

    /// Compares `metric` against `tol` and builds the matching result.
    pub fn from_tol(metric: f64, tol: f64) -> Self {
        if metric <= tol {
            Self::pass(metric, Some(tol))
        } else {
            Self::fail(metric, Some(tol))
        }
    }
}

/// One case's row of the report.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CaseResult {
    /// Global case id.
    pub case: String,
    /// Aggregate status (`PASS` / `FAIL` / `SKIP`), derived from `gates`.
    pub status: String,
    /// Per-gate outcomes, keyed by canonical gate name.
    pub gates: BTreeMap<String, GateResult>,
    /// Per-stage timings, keyed by canonical stage name.
    pub stages: BTreeMap<String, StageTiming>,
    /// Wall-clock seconds spent on the whole case, when tracked.
    pub wall_s: Option<f64>,
}

impl CaseResult {
    /// Creates an empty row for `qualid` (status starts as `SKIP` until gates land).
    pub fn new(qualid: impl Into<String>) -> Self {
        Self {
            case: qualid.into(),
            status: "SKIP".to_string(),
            gates: BTreeMap::new(),
            stages: BTreeMap::new(),
            wall_s: None,
        }
    }

    /// Aggregate status: no gates -> `SKIP`; any `FAIL` -> `FAIL`; all `SKIP` -> `SKIP`; else `PASS`.
    pub fn derive_status(&self) -> &'static str {
        if self.gates.is_empty() {
            "SKIP"
        } else if self.gates.values().any(|g| g.status == "FAIL") {
            "FAIL"
        } else if self.gates.values().all(|g| g.status == "SKIP") {
            "SKIP"
        } else {
            "PASS"
        }
    }

    /// Writes [`Self::derive_status`] into the `status` field.
    pub fn sync_status(&mut self) -> &str {
        self.status = self.derive_status().to_string();
        &self.status
    }

    /// Records one gate outcome.
    pub fn set_gate(&mut self, name: impl Into<String>, gate: GateResult) -> &mut Self {
        self.gates.insert(name.into(), gate);
        self
    }

    /// Records a stage timing.
    pub fn set_stage(&mut self, stage: impl Into<String>, timing: StageTiming) -> &mut Self {
        self.stages.insert(stage.into(), timing);
        self
    }
}

/// Backend identity block of the report.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HarnessInfo {
    /// Backend name, e.g. `mfs`.
    pub backend: String,
    /// Implementation language, e.g. `rust`.
    pub language: String,
    /// Platform string.
    pub platform: String,
}

/// Cross-language benchmark report (matches `data/schema/report.schema.json`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Report {
    /// Report schema version.
    pub schema_version: String,
    /// Standard contract version.
    pub contract_version: String,
    /// Backend identity.
    pub harness: HarnessInfo,
    /// One row per case.
    pub results: Vec<CaseResult>,
}

/// Report schema version mandated by the standard.
pub const REPORT_SCHEMA_VERSION: &str = "1.0";

impl Report {
    /// Creates a report for one backend, failing loudly on contract drift.
    pub fn new(backend: &str, language: &str, platform: &str) -> Self {
        let contract_version =
            check_contract().unwrap_or_else(|e| panic!("benchmark contract: {e}"));
        Self {
            schema_version: REPORT_SCHEMA_VERSION.to_string(),
            contract_version,
            harness: HarnessInfo {
                backend: backend.to_string(),
                language: language.to_string(),
                platform: platform.to_string(),
            },
            results: Vec::new(),
        }
    }

    /// Appends a case row and returns it for filling.
    pub fn add_case(&mut self, qualid: &str) -> &mut CaseResult {
        self.results.push(CaseResult::new(qualid));
        self.results.last_mut().expect("just pushed")
    }

    /// Looks up an already-added case row.
    pub fn result(&self, qualid: &str) -> Option<&CaseResult> {
        self.results.iter().find(|r| r.case == qualid)
    }

    /// Writes the report as pretty JSON, refreshing aggregate statuses first.
    pub fn write(&self, path: &Path) -> Result<()> {
        let mut view = self.clone();
        for result in &mut view.results {
            result.sync_status();
        }
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        std::fs::write(path, format!("{}\n", serde_json::to_string_pretty(&view)?))
            .map_err(|e| format!("cannot write {}: {e}", path.display()).into())
    }

    /// `2` when any case failed (contract: failing cases are reported, not timed).
    pub fn exit_code(&self) -> i32 {
        if self
            .results
            .iter()
            .any(|r| r.status == "FAIL" || r.derive_status() == "FAIL")
        {
            2
        } else {
            0
        }
    }
}
