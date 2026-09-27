# Contributing

Thanks for contributing to `mfs`.

## Before You Start

- Check [CHANGELOG.md](CHANGELOG.md) and [docs/architecture.md](docs/architecture.md)
  for the current module layout and known limitations.
- Prefer small, focused changes over mixed refactors

## Toolchain & Gates

- Rust **≥ 1.86** (see `rust-version` in `Cargo.toml`; CI tests stable + 1.86).
- CI gates that must stay green:
  ```
  cargo clippy --all-targets -- -D warnings
  RUSTDOCFLAGS="-D warnings" cargo doc --no-deps
  cargo test --all-targets
  ```
- Formatting: the codebase is **intentionally not `cargo fmt`-gated** and is not
  rustfmt-clean. Do NOT run a repo-wide `cargo fmt` — keep your diff's local
  style consistent with the surrounding code instead. (A fmt baseline + gate is
  a candidate follow-up; it would be a dedicated reformat PR, not drive-by.)

## Development Workflow

1. Create a branch for your change.
2. Make the smallest coherent change that solves one problem.
3. Add or update tests.
4. Run:

```powershell
cargo check --all-targets
cargo test
```

5. Update docs when behavior, architecture, or public API changes.

## Coding Expectations

- Keep domain semantics explicit.
- Prefer validated constructors over loose structs.
- Keep solver and backend details behind internal boundaries.
- Avoid leaking temporary implementation choices into the public API.
- Add comments only when mathematical intent is not obvious from the code.

## Pull Requests

A good pull request usually includes:

- a short problem statement
- the chosen approach
- tests or fixtures covering the change
- notes on limitations or follow-up work

## Scope Guidance

Good fits for contributions:

- bug fixes
- test improvements
- benchmark fixtures
- frequency mapping improvements
- generalized Chebyshev helper integration
- response-solver verification
- documentation improvements

Please avoid mixing unrelated cleanup into the same change set.

## Questions

If a change will significantly affect public API shape, polynomial
representation, or matrix artifact semantics, document the tradeoff clearly in
the pull request before expanding scope.
