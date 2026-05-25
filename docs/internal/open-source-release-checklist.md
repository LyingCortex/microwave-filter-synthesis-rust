# Open Source Release Checklist

This checklist tracks the remaining work needed before or immediately after
publishing `mfs` on GitHub as a public open-source library.

## Repository Setup

- Create the public GitHub repository
- Add the real repository URL to `Cargo.toml`
- Add the real repository URL to `README.md`
- Enable GitHub Actions
- Configure branch protection for the default branch if desired

## Licensing And Governance

- Confirm the dual-license choice is acceptable
- Add maintainers or organization ownership on GitHub
- Decide how security disclosures should be reported

## Release Quality

- Keep `cargo test` green
- Keep `cargo check --all-targets` green
- Keep docs build green in CI
- Keep packaging dry-run green in CI
- Freeze a minimum benchmark/reference set for public validation
- Decide whether the current `0.1.0` surface should be published to `crates.io`

## Metadata

- Add `repository` in `Cargo.toml`
- Optionally add `homepage` and `documentation`
- Add badges to `README.md` once URLs are known

Current status:

- `cargo package --allow-dirty` already succeeds
- CI now covers lint, tests, docs, packaging, and a scheduled dependency audit
- repository and homepage URLs are configured
- docs.rs URL is configured for crate documentation once published

## Public Communication

- Write the first GitHub release notes
- Summarize current capabilities and limitations clearly
- Mention that generalized Chebyshev integration is still in progress

## Nice-To-Have

- Add issue labels and project boards
- Add benchmark workflow or nightly validation
- Add example result snapshots or reference fixtures
