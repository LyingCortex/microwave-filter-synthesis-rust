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
fn build_network_roundtrip_via_file() -> Result<()> {
    // Synthesize -> build_network (NetworkBuilder) -> save .s2p ->
    // read_touchstone (Network::new) -> compare. This round-trip exercises
    // both the builder's serialization and the 0.14.5 parser.
    let design = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6)
        .zeros_hz([6.5e9, 7.0e9])
        .synthesize()?;
    let response = design.response(6.5e9, 7.0e9, 21)?;
    let net = mfs::touchstone::build_network(&response, 50.0)?;

    // Live network sanity (NetworkBuilder built something usable, not text).
    assert_eq!(net.rank, 2);
    assert_eq!(net.f.len(), 21);
    let center = net.try_s_ri_at(10, 2, 1).expect("network point 10 S21");
    let expected = response.samples[10];
    assert!((center.re - expected.s21_re).abs() < 1e-9);
    assert!((center.im - expected.s21_im).abs() < 1e-9);

    let tmp = std::env::temp_dir().join("mfs_test_build_network.s2p");
    net.save(tmp.to_str().unwrap()).unwrap();

    let loaded = touchstone::read_touchstone(&tmp)?;
    assert_eq!(loaded.samples.len(), 21);
    for (orig, read) in response.samples.iter().zip(loaded.samples.iter()) {
        assert!((orig.frequency_hz - read.frequency_hz).abs() < 1.0); // GHz scaling round-trip
        assert!((orig.s11_re - read.s11_re).abs() < 1e-9);
        assert!((orig.s11_im - read.s11_im).abs() < 1e-9);
        assert!((orig.s21_re - read.s21_re).abs() < 1e-9);
        assert!((orig.s21_im - read.s21_im).abs() < 1e-9);
    }

    std::fs::remove_file(tmp).ok();
    Ok(())
}

#[test]
fn build_network_roundtrip_in_memory() -> Result<()> {
    // build_network -> to_touchstone_string -> Network::from_str -> compare:
    // a serialization round-trip with no filesystem involved.
    let design = FilterDesign::bandpass(6, 23.0, 6.75e9, 300e6).synthesize()?;
    let net = design.build_network()?;
    assert_eq!(net.rank, 2);
    assert!(net.f.len() > 100); // default_response is a 201-point sweep

    let text = net.to_touchstone_string().expect("serialize network to text");
    let reparsed = ::touchstone::Network::from_str("roundtrip.s2p", &text)
        .expect("parse network text back");
    assert_eq!(reparsed.rank, 2);
    assert_eq!(reparsed.f.len(), net.f.len());
    for i in [0, reparsed.f.len() / 2, reparsed.f.len() - 1] {
        let a = net.try_s_ri_at(i, 1, 1).expect("net S11");
        let b = reparsed.try_s_ri_at(i, 1, 1).expect("reparsed S11");
        assert!((a.re - b.re).abs() < 1e-9, "S11 re at {i}: {a:?} vs {b:?}");
        assert!((a.im - b.im).abs() < 1e-9, "S11 im at {i}");
        let c = net.try_s_ri_at(i, 2, 1).expect("net S21");
        let d = reparsed.try_s_ri_at(i, 2, 1).expect("reparsed S21");
        assert!((c.re - d.re).abs() < 1e-9, "S21 re at {i}: {c:?} vs {d:?}");
        assert!((c.im - d.im).abs() < 1e-9, "S21 im at {i}");
    }
    Ok(())
}

#[test]
fn build_network_rejects_empty_response() {
    use mfs::response::SParameterResponse;
    let empty = SParameterResponse { samples: Vec::new() };
    let err = mfs::touchstone::build_network(&empty, 50.0).unwrap_err();
    assert!(matches!(err, mfs::MfsError::PreconditionViolation(_)));
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
