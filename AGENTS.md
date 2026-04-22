# Repository Guidelines

## Project Structure & Module Organization

This repository is a Rust 2024 library crate named `mfs`. Core code lives in `src/`, with public exports coordinated from `src/lib.rs` and common user-facing imports in `src/prelude.rs`. Major modules are organized by synthesis stage: `spec/` for filter specifications, `approx/` for polynomial and generalized Chebyshev approximation, `synthesis/` for orchestration and matrix recovery, `matrix/` for coupling-matrix artifacts, `transform/` for topology transforms, `response/` and `verify/` for evaluation checks, and `output/` for reports. Integration tests and JSON fixtures live in `tests/`. Runnable examples live in `examples/`. Design notes and roadmaps live in `docs/`.

## Build, Test, and Development Commands

- `cargo check --all-targets`: type-check the library, tests, and examples.
- `cargo test --all-targets`: run all unit, integration, and example-target tests.
- `cargo fmt --all -- --check`: verify Rust formatting matches `rustfmt`.
- `cargo clippy --all-targets --all-features -- -D warnings`: run lint checks with warnings treated as errors.
- `cargo doc --no-deps`: build crate documentation without dependency docs.
- `cargo run --example quickstart_report`: run a representative example locally.

## Coding Style & Naming Conventions

Use idiomatic Rust formatted by `rustfmt`; do not hand-align code. Module and file names use `snake_case`, public types use `UpperCamelCase`, and functions use `snake_case`. Keep domain boundaries explicit: prefer validated constructors and typed artifacts over loosely structured data. Comments should explain mathematical intent or non-obvious synthesis steps, not restate the code.

## Testing Guidelines

Add or update tests whenever behavior, numerical output, topology shape, or public API contracts change. Place cross-module behavior tests in `tests/*.rs`; keep fixture data in `tests/*.json` or reusable fixture helpers under `src/fixtures/`. Follow existing descriptive test filenames such as `transform_invariance.rs` and `section_synthesis_facade.rs`. For focused work, run the relevant test file first, then `cargo test --all-targets` before submitting.

## Commit & Pull Request Guidelines

Recent commits use short imperative summaries, for example `Refactor generalized Chebyshev synthesis flow` or `Expand synthesis pipeline and literature fixture coverage`. Keep commits focused and avoid mixing unrelated cleanup with feature work. Pull requests should include a short summary, motivation, testing notes, limitations or follow-up work, and any public API or matrix artifact implications. Use the checklist in `.github/PULL_REQUEST_TEMPLATE.md`.

## Security & Configuration Tips

Do not commit generated build output or local experiment artifacts. Keep `Cargo.lock` updated for reproducible CI. Report security issues using `SECURITY.md` rather than public issue threads.
