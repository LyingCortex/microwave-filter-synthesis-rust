use mfs::prelude::*;
use mfs::touchstone::{self, FreqUnit, DataFormat, TouchstoneConfig};

#[test]
fn export_bandpass_to_touchstone_string() -> Result<()> {
    let design = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6)
        .zeros_hz([6.5e9, 7.0e9])
        .synthesize()?;

    let response = design.response(6.5e9, 7.0e9, 11)?;
    let config = TouchstoneConfig::default();
    let content = touchstone::to_touchstone_string(&response, &config)?;

    // Verify format
    assert!(content.contains("# GHZ S RI R 50"));
    assert!(content.lines().count() > 11); // header + 11 data lines
    println!("{content}");
    Ok(())
}

#[test]
fn export_db_format() -> Result<()> {
    let design = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
    let response = design.response(6.5e9, 7.0e9, 5)?;

    let config = TouchstoneConfig {
        freq_unit: FreqUnit::MHz,
        format: DataFormat::DB,
        impedance: 50.0,
        version: mfs::touchstone::TouchstoneVersion::V1,
        comments: vec!["Test filter".to_string()],
    };
    let content = touchstone::to_touchstone_string(&response, &config)?;

    assert!(content.contains("# MHZ S DB R 50"));
    assert!(content.contains("! Test filter"));
    println!("{content}");
    Ok(())
}

#[test]
fn save_and_read_roundtrip() -> Result<()> {
    let design = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
    let response = design.response(6.5e9, 7.0e9, 21)?;

    let tmp = std::env::temp_dir().join("mfs_test_roundtrip.s2p");
    let config = TouchstoneConfig::default();
    touchstone::write_touchstone(&response, &config, &tmp)?;

    // Read back
    let loaded = touchstone::read_touchstone(&tmp)?;
    assert_eq!(loaded.samples.len(), 21);

    // Compare S21 magnitudes
    for (orig, read) in response.samples.iter().zip(loaded.samples.iter()) {
        let diff = (orig.s21_mag() - read.s21_mag()).abs();
        assert!(diff < 1e-6, "S21 mismatch: orig={}, read={}", orig.s21_mag(), read.s21_mag());
    }

    std::fs::remove_file(tmp).ok();
    Ok(())
}

#[test]
fn filter_design_save_touchstone() -> Result<()> {
    let design = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6)
        .zeros_hz([6.5e9, 7.0e9])
        .synthesize()?;

    let tmp = std::env::temp_dir().join("mfs_test_design.s2p");
    design.save_touchstone(&tmp)?;

    assert!(tmp.exists());
    let content = std::fs::read_to_string(&tmp).unwrap();
    assert!(content.contains("# GHZ S RI R 50"));
    assert!(content.contains("MFS v"));

    std::fs::remove_file(tmp).ok();
    Ok(())
}

#[test]
fn filter_design_to_touchstone_builder() -> Result<()> {
    let design = FilterDesign::bandpass(6, 23.0, 6.75e9, 300e6)
        .zeros_hz([6.4e9, 7.0e9])
        .synthesize()?;

    let content = design.to_touchstone()
        .freq_unit(FreqUnit::GHz)
        .format(DataFormat::DB)
        .comment("Custom comment")
        .build()?;

    assert!(content.contains("# GHZ S DB R 50"));
    assert!(content.contains("! Custom comment"));
    assert!(content.contains("! MFS v"));
    Ok(())
}
