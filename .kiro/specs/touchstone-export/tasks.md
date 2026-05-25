# Implementation Tasks

## Task 1: Create touchstone module with types
- [ ] Create `src/touchstone/mod.rs` with public exports
- [ ] Create `src/touchstone/types.rs` with `FreqUnit`, `DataFormat`, `TouchstoneVersion`, `TouchstoneConfig`, `TouchstoneFile`
- [ ] Add `Default` impl for `TouchstoneConfig` (GHz, RI, 50Ω, V1)
- [ ] Register `pub mod touchstone` in `lib.rs`
- [ ] Add `FreqUnit` scaling factor method: `fn scale(&self) -> f64`

## Task 2: Implement Touchstone v1.0 writer
- [ ] Create `src/touchstone/writer.rs`
- [ ] Implement option line formatting: `# {unit} S {format} R {impedance}`
- [ ] Implement RI data row formatting (8 values per line for 2-port)
- [ ] Implement MA data row formatting (magnitude + angle in degrees)
- [ ] Implement DB data row formatting (dB + angle in degrees)
- [ ] Implement comment line formatting (prefixed with `!`)
- [ ] Implement `to_touchstone_string(response, config) -> Result<String>`
- [ ] Add auto-generated header comments (order, RL, center, BW, MFS version)
- [ ] Test: round-trip RI format preserves values within 1e-10

## Task 3: Implement Touchstone v2.0 writer
- [ ] Add `[Version] 2.0` header
- [ ] Add `[Number of Ports] 2` declaration
- [ ] Add `[Number of Frequencies] N` declaration
- [ ] Wrap data in `[Network Data]` / `[End]` markers
- [ ] Test: v2.0 output contains all required sections

## Task 4: Implement file write support
- [ ] Implement `write_touchstone(response, config, path) -> Result<()>`
- [ ] Validate file extension is `.s2p`
- [ ] Map IO errors to `MfsError::PreconditionViolation`
- [ ] Test: write and read back produces identical content

## Task 5: Implement Touchstone parser
- [ ] Create `src/touchstone/parser.rs`
- [ ] Parse option line (`#` line) extracting unit, format, impedance
- [ ] Parse v1.0 data rows (handle 2-port 8-column format)
- [ ] Parse v2.0 headers and `[Network Data]` section
- [ ] Skip comment lines (`!` prefix)
- [ ] Convert parsed data back to internal Hz + (re, im) representation
- [ ] Return descriptive errors with line numbers for malformed input
- [ ] Implement `parse_touchstone(content) -> Result<TouchstoneFile>`
- [ ] Test: parse(export(response)) ≈ response (round-trip property)

## Task 6: Integrate with FilterDesign API
- [ ] Add `TouchstoneBuilder` struct with fluent configuration
- [ ] Add `FilterDesign::to_touchstone() -> TouchstoneBuilder`
- [ ] Add `FilterDesign::save_touchstone(path) -> Result<()>` (default config)
- [ ] `TouchstoneBuilder::build() -> Result<String>`
- [ ] `TouchstoneBuilder::save(path) -> Result<()>`
- [ ] Auto-populate comments from design metadata (order, RL, center, BW)
- [ ] If bandpass design: auto-compute response for export
- [ ] If prototype design: require explicit frequency range or error
- [ ] Test: `design.save_touchstone("test.s2p")` produces valid file

## Task 7: Integrate with Python bindings
- [ ] Add `to_touchstone` method to `PyFilterDesign`
- [ ] Add `save_touchstone` method to `PyFilterDesign`
- [ ] Accept keyword args: `freq_unit`, `format`, `impedance`, `version`
- [ ] Map string args to enums: "GHz" → `FreqUnit::GHz`, "DB" → `DataFormat::DB`
- [ ] Raise `ValueError` for invalid arguments
- [ ] Test: Python round-trip produces valid output (manual verification)

## Task 8: Add to prelude and update docs
- [ ] Export `FreqUnit`, `DataFormat` from prelude (if needed by users)
- [ ] Update `docs/api-reference.md` with Touchstone section
- [ ] Update `README.md` with Touchstone example
- [ ] Update `CHANGELOG.md` for v0.2.0
- [ ] Add example: `examples/touchstone_export.rs`

## Task 9: Final validation
- [ ] Run all existing tests (ensure no regressions)
- [ ] Verify exported .s2p file imports correctly in a reference tool
- [ ] Verify round-trip: export → parse → compare
- [ ] Verify Python bindings compile and work
- [ ] Tag v0.2.0
