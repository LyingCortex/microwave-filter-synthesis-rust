# Cameron 2003 Coupling-Matrix Algorithm Notes

## Purpose

This note summarizes the coupling-matrix synthesis and reconfiguration flow
described by Richard J. Cameron in the 2003 paper on advanced coupling matrix
synthesis techniques for microwave filters, and supplements it with the 2011
ZTE Communications overview article "Advanced Synthesis Techniques for
Microwave Filters."

The goal of this document is not to restate every derivation from the paper,
but to capture the algorithmic steps clearly enough that we can:

- review the Rust implementation against the intended procedure
- identify which matrix form is required before each transformation
- design tests that verify topology changes without changing the response

## Scope

This note focuses on the parts of the paper that matter most to this codebase:

- the `N+2` normalized coupling-matrix model
- synthesis of the canonical transversal matrix
- reconfiguration through similarity transforms
- conversion to folded and arrow canonical forms
- generation and shifting of trisections
- merging adjacent trisections into higher-order sections such as quartets

It does not attempt to fully reproduce the polynomial-generation part of the
paper. For the approximation side, Cameron's generalized Chebyshev recurrence
and the earlier 1999 paper remain the primary references.

The ZTE article is still useful here because it restates several points in a
more implementation-oriented way:

- why coupling-matrix methods displaced classical element extraction in
  microwave filter work
- how the canonical transversal matrix is obtained by equating two short-
  circuit admittance descriptions of the same network
- why folded and arrow forms are better treated as realizable descendants of
  the transversal matrix rather than as unrelated matrix templates
- how trisections, quartets, and box sections fit into a single dependency
  chain

## Historical and Engineering Context

The ZTE article adds important background for why this whole workflow exists.

Before the 1970s, most filter synthesis methods were based on extracting
lumped elements or transmission-line sections directly from filtering
polynomials. Cameron's summary explains that this became inadequate once
satellite and terrestrial telecom systems pushed into more crowded spectrum
with tighter requirements on:

- close-to-band rejection
- in-band insertion loss and group delay linearity
- physically realizable microwave cross-couplings

The coupling-matrix method became the preferred tool because it gives a direct
correspondence between matrix entries and physical coupling elements, while
also allowing the matrix to be reconfigured through similarity transforms
without restarting the synthesis from scratch.

## Notation and Conventions

The paper uses the normalized low-pass bandpass-prototype viewpoint:

- filter degree: `N`
- matrix size: `(N + 2) x (N + 2)`
- nodes: `S, 1, 2, ..., N, L`
- source node index: `S`
- load node index: `L`
- diagonal terms: resonator frequency offsets / detunings
- off-diagonal terms: couplings between source, load, and resonators

In the ZTE article, the underlying prototype is described as a generalized
multicoupled bandpass network. In its original form this gives an `N x N`
coupling matrix for the resonator loops alone. Replacing the transformer model
with immittance inverters and adding explicit input/output inverters extends
the model to the familiar `(N + 2) x (N + 2)` microwave form.

The canonical `N+2` matrix contains:

- source couplings in the first row/column
- load couplings in the last row/column
- resonator self-couplings on the diagonal
- optional direct source-load coupling `M_SL` for fully canonical responses

For asymmetric filtering characteristics, the article also notes that each
resonator loop may include a series frequency-invariant reactance (FIR). In
software terms, this is a reminder that diagonal terms are not just nuisance
detunings: they are part of the mechanism that lets the matrix represent
asymmetric responses.

For practical coding work, the most important distinction is the matrix form:

- `transversal`: resonators are independent branches between source and load
- `folded`: nearest-neighbor main line plus a structured set of cross-couplings
- `arrow` or `wheel`: main line plus couplings from each resonator to the load
- `trisection cascade`: an arrow-derived arrangement where each trisection
  realizes one transmission zero

## Overall Algorithm

The paper's workflow can be summarized as:

1. Generate the filtering polynomials.
2. Build a canonical `N+2` transversal coupling matrix from those polynomials.
3. Apply analytically chosen similarity rotations to convert the transversal
   matrix into a realizable canonical form, usually folded or arrow.
4. If an advanced configuration is desired, continue from the correct
   canonical form:
   - from `folded` or `transversal` to `arrow`
   - from `arrow` to cascaded trisections
   - from adjacent trisections to quartets, quintets, or box sections
5. Realize the final sparse topology in the physical microwave structure.

The ZTE article makes the physical interpretation explicit:

1. generate `S21(s)` and `S11(s)` polynomials from the electrical
   specifications
2. synthesize a canonical coupling matrix
3. reconfigure the matrix until the non-zero entries match the couplings that
   the chosen structure can actually realize
4. compute the dimensions of the physical coupling elements from those matrix
   values

One key point from the paper is procedural:

- advanced sections are not generated from an arbitrary matrix
- they are generated from a specific canonical form reached by similarity
  transforms

That precondition matters when reviewing code.

## Step 1: Build the `N+2` Transversal Matrix

The transversal matrix is the easiest canonical form to synthesize directly.
It consists of `N` first-order branches in parallel between source and load.
There are no inter-resonator couplings between different branches.

Structurally, the matrix contains:

- `M_Sk` for source-to-resonator couplings
- `M_kL` for resonator-to-load couplings
- diagonal terms `M_kk`
- optionally `M_SL` for the fully canonical case
- all other entries zero

In Cameron's procedure, the matrix is obtained by writing the network
admittance matrix in two ways:

- from the prescribed scattering-function polynomials
- from the circuit representation of the transversal array

Then the two descriptions are equated, yielding the coupling values of the
transversal matrix.

The ZTE article gives a more concrete statement of this construction:

- build the two-port short-circuit admittance matrix `[Y_N]` from the target
  scattering-function polynomials
- build the same `[Y_N]` from the transversal circuit made of parallel
  first-order sections
- equate the two forms to solve for `M_Sk`, `M_kL`, diagonal terms, and
  optional `M_SL`

For implementation work, the important takeaway is:

- the transversal matrix is a synthesis target, not merely a placeholder shape
- if we synthesize it through residues, the result still needs to obey the
  exact canonical sparsity pattern

The article also re-emphasizes the minimum-path implication of `M_SL`:

- the maximum number of finite-position transmission zeros is
  `nfz_max = N - n_min`
- `n_min` is the number of resonator nodes in the shortest source-to-load path
- in a fully canonical network with direct source-load coupling,
  `n_min = 0`, so up to `N` finite transmission zeros are realizable

## Step 2: Similarity Transformations

Once the transversal matrix is available, Cameron uses similarity transforms
to preserve the network response while changing the topology.

The transformed matrix is:

```text
M' = R * M * R^T
```

where `R` is an orthogonal rotation matrix that differs from identity only in
the two pivot rows/columns involved in the rotation.

For a two-pivot rotation between indices `i` and `j`, the non-trivial entries
of `R` are:

```text
R_ii = cos(theta)
R_jj = cos(theta)
R_ij = -sin(theta)
R_ji =  sin(theta)
```

The angle `theta` is chosen so that a selected entry of `M'` is annihilated.
This is the fundamental mechanism used throughout the paper:

- clear an unwanted coupling
- create a wanted coupling elsewhere
- preserve the electrical response

This preservation property is essential. A correct topology transform should
change sparsity, but not the ideal `S11` and `S21` produced by the matrix.

The ZTE article phrases the same point in eigenstructure terms:

- similarity transforms preserve the eigenvalues and eigenvectors of `M`
- therefore the transformed matrix has the same ideal transfer and reflection
  characteristics as the original one

## Step 3: Transversal to Folded

The folded matrix is one of the standard canonical forms used in microwave
filter realization.

Conceptually, the transform proceeds by a formal sequence of rotations that:

- remove the excess source couplings of the transversal form
- remove the excess load couplings in a mirrored way
- concentrate the remaining non-zero couplings into a folded pattern

The folded form keeps:

- the main-line couplings
- the input and output couplings
- selected cross-couplings that make finite transmission zeros realizable

Implementation notes:

- folded conversion is a sequence problem, not a one-shot formula
- the order of annihilations matters
- sign normalization may be applied afterward so nearest-neighbor couplings
  match the preferred realization convention

## Step 4: Transversal to Arrow

The arrow or wheel matrix is another canonical form described in the paper.
It retains the main-line chain and places the cross-couplings into the final
row and column, visually forming an arrow toward the lower-right corner.

The arrow form is especially important because it is the starting point for
the trisection synthesis procedure.

The paper's logic is:

1. synthesize the canonical transversal matrix
2. transform it into the arrow canonical matrix with a formal rotation series
3. generate one trisection per transmission zero inside that arrow matrix
4. shift each trisection leftward to form a cascade

That means:

- trisection generation is not a direct transform from an arbitrary matrix
- the arrow matrix is the required intermediate form

The ZTE article also gives a useful geometric description:

- the main-line couplings form the rim of the "wheel"
- source/load cross-couplings form the spokes
- in matrix form, the cross-couplings gather in the last row and column,
  producing the visual "arrow" shape

## Step 5: Trisection Synthesis in the Arrow Matrix

A trisection is a three-node section able to realize one transmission zero.
It may be:

- internal
- input-coupled
- output-coupled
- or, in the degree-1 canonical case, equivalent to direct source-load
  coupling

The paper states that the basis of trisection synthesis is a zero-determinant
condition evaluated at the target transmission-zero frequency `omega_0`.

In practical terms:

- choose one prescribed transmission zero
- condition the arrow matrix by a rotation at pivot `[N-1, N]`
- choose the angle so the trisection condition is satisfied for that zero
- this creates the first trisection at the tail end of the arrow matrix
- shift that trisection leftward by additional rotations

Then repeat for the next transmission zero until the full cascade is formed.

Algorithmically, the pattern is:

1. start from arrow matrix `M^(0)`
2. condition the tail section for the first zero
3. shift the created trisection left to its desired location
4. return to the tail region for the next zero
5. repeat until every finite transmission zero is assigned to one trisection

The paper's high-level rule is clear even without reproducing every page-level
equation:

- first create the trisection in the arrow tail
- then move it through the network by rotations

The ZTE article makes two additional points worth keeping in mind:

- a trisection can be internal, input-coupled, output-coupled, or degenerate
  to the direct `M_SL` case in a degree-1 canonical network
- each trisection realizes exactly one transmission zero by the minimum-path
  rule, which is why trisection cascades are especially useful for asymmetric
  characteristics

## Step 6: Merge Adjacent Trisections into Quartets or Higher Sections

Once a cascade of trisections exists, adjacent trisections can be merged into
larger sections.

The paper explicitly discusses:

- quartets formed by merging two trisections
- quintets formed by merging three trisections

The reason to do this is practical:

- some technologies realize quartets or higher sections more naturally than
  isolated trisections

The algorithmic idea is:

- build the trisections first
- then use further rotations to eliminate one internal coupling and create the
  desired merged section

So the dependency chain is:

```text
transversal -> arrow -> trisection cascade -> quartet / quintet / box variants
```

That ordering should be preserved in implementation.

## Box and Extended-Box Sections

The paper also describes how a synthesized trisection can be transformed into
box or extended-box sections.

The relevant pattern is:

- start from an already synthesized trisection
- apply a cross-pivot rotation to annihilate one of the trisection's main-line
  couplings
- this creates a different sparse topology that may be easier to realize in
  some technologies, including dual-mode structures

The key point is the same as above:

- the box section is derived from a trisection already embedded in the matrix
- it is not synthesized as a standalone primitive from scratch

The ZTE article adds two implementation-relevant details:

- the basic box section is obtained by a cross-pivot rotation that annihilates
  the second main-line coupling of the trisection
- in the resulting box section, one coupling is always negative, regardless of
  the sign of the original trisection cross-coupling

It also notes why box sections are attractive in practice:

- the basic box realizes one transmission zero without a diagonal cross-
  coupling
- this can be convenient for dual-mode technologies
- complementary responses can be obtained by conjugating the self-couplings,
  which in hardware corresponds to retuning resonators rather than rebuilding
  couplings

For extended box sections, the article gives a direct minimum-path rule:

- 4th-degree extended box: up to 1 finite transmission zero
- 6th-degree extended box: up to 2
- 8th-degree extended box: up to 3
- in general: up to `(N - 2) / 2`

## Minimum-Path Rule

The paper repeatedly uses the minimum-path rule to reason about how many
finite-position transmission zeros a topology can realize.

In implementation terms, this rule is useful as a sanity check:

- a topology with too short a source-to-load path can realize more zeros
- a topology with too long a minimum path cannot realize the requested number
  of zeros without extra couplings

Examples:

- a fully canonical network with direct source-load coupling can realize up to
  `N` finite transmission zeros
- each trisection realizes one transmission zero
- a basic box section realizes one transmission zero
- an extended box realizes up to `(N - 2) / 2` transmission zeros
- larger sections inherit their capability from the number of embedded
  trisections and the resulting minimum path

## What Must Be Preserved During Reconfiguration

Every similarity-transform-based reconfiguration should preserve:

- the ideal transfer and reflection characteristics
- the eigenstructure of the matrix
- reciprocity and symmetry of the real-valued coupling matrix

For software verification, this suggests three classes of tests:

1. structural tests
   Check that the intended couplings are annihilated or created.
2. invariance tests
   Compare ideal responses before and after the rotation sequence.
3. precondition tests
   Check that algorithms expecting `arrow` or `trisection` input reject or
   document other starting forms.

## Implementation Checklist

When comparing code to the Cameron procedure, the most important questions are:

- Is the transversal matrix genuinely canonical, or only approximately shaped
  like one?
- Does each topology transform start from the canonical form assumed by the
  paper?
- Are trisections created from an arrow matrix first and only then shifted?
- Are quartets formed by merging existing trisections rather than by skipping
  the trisection stage?
- Do rotation sequences preserve the response numerically?
- Are sign conventions and indexing conventions documented explicitly?

## Suggested Mapping to This Repository

For future review work in this repo, a practical mental map is:

- approximation code:
  builds the filtering polynomials
- matrix synthesis code:
  should build the canonical transversal matrix
- topology transforms:
  should convert transversal into folded or arrow by similarity rotations
- advanced section extraction:
  should start from the required canonical precursor, especially arrow for
  trisections
- response solver:
  should be used to verify response invariance after each transform sequence

## References

1. Richard J. Cameron, "Advanced coupling matrix synthesis techniques for
   microwave filters," IEEE Transactions on Microwave Theory and Techniques,
   vol. 51, no. 1, pp. 1-10, Jan. 2003. DOI:
   `10.1109/TMTT.2002.806937`
2. Richard J. Cameron, "General coupling matrix synthesis methods for
   Chebyshev filtering functions," IEEE Transactions on Microwave Theory and
   Techniques, vol. 47, no. 4, pp. 433-442, Apr. 1999. DOI:
   `10.1109/22.754877`
3. Richard J. Cameron, "Advanced Synthesis Techniques for Microwave Filters,"
   ZTE Communications, 2011. This is a useful secondary summary of the
   historical motivation, `N+2` matrix interpretation, short-circuit
   admittance construction of the transversal matrix, and the folded, arrow,
   trisection, quartet, and box-section workflow.
