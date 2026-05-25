use mfs::prelude::*;

#[test]
fn prelude_supports_new_filter_design_api() -> Result<()> {
    let design = FilterDesign::prototype(4, 20.0)
        .zeros([-1.5])
        .synthesize()?;

    assert_eq!(design.order(), 4);
    assert_eq!(design.matrix().order(), 4);
    Ok(())
}

#[test]
fn prelude_supports_bandpass_design() -> Result<()> {
    let design = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6)
        .zeros_hz([6.5e9, 7.0e9])
        .synthesize()?;

    assert_eq!(design.order(), 4);
    let response = design.response(6.5e9, 7.0e9, 11)?;
    assert_eq!(response.samples.len(), 11);
    Ok(())
}

#[test]
fn prelude_supports_topology_transforms() -> Result<()> {
    let design = FilterDesign::prototype(4, 20.0)
        .zeros([-1.5, 1.5])
        .synthesize()?;

    let folded = design.to_folded()?;
    assert_eq!(folded.order(), 4);
    Ok(())
}
