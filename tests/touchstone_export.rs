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
fn cascade_with_identity_is_identity() -> Result<()> {
    let design = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
    let net = design.build_network()?;

    // Identity 2-port (S = [[0,1],[1,0]]) on the same frequency grid:
    // cascading it onto `net` must leave `net` unchanged (ABCD = I).
    let mut id_builder = ::touchstone::NetworkBuilder::new("identity.s2p", 2).z0(50.0);
    for &freq in net.f.iter() {
        let zero = ::touchstone::Complex { re: 0.0, im: 0.0 };
        let one = ::touchstone::Complex { re: 1.0, im: 0.0 };
        let matrix = ::touchstone::SMatrix {
            rank: 2,
            data: vec![vec![zero, one], vec![one, zero]],
        };
        id_builder = id_builder.point(freq, matrix);
    }
    let identity = id_builder.build().expect("identity network builds");

    let cascaded = mfs::touchstone::cascade(&net, &identity)?;
    assert_eq!(cascaded.rank, 2);
    assert_eq!(cascaded.f.len(), net.f.len());
    for i in [0, net.f.len() / 2, net.f.len() - 1] {
        for (port, label) in [(1, "S11"), (2, "S21")] {
            let orig = net.try_s_ri_at(i, port, 1).unwrap();
            let casc = cascaded.try_s_ri_at(i, port, 1).unwrap();
            assert!(
                (orig.re - casc.re).abs() < 1e-6 && (orig.im - casc.im).abs() < 1e-6,
                "{label} at point {i}: {orig:?} vs {casc:?}"
            );
        }
    }
    Ok(())
}

#[test]
fn cascade_rejects_mismatched_grids() -> Result<()> {
    // Different band centers -> different default response grids. The wrapper
    // must return a typed error, not silently pair index-wise (upstream
    // Network::cascade does the latter).
    let a = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
    let b = FilterDesign::bandpass(4, 20.0, 7.0e9, 300e6).synthesize()?;
    let net_a = a.build_network()?;
    let net_b = b.build_network()?;
    let err = mfs::touchstone::cascade(&net_a, &net_b).unwrap_err();
    assert!(
        matches!(err, mfs::MfsError::PreconditionViolation(_)),
        "mismatched grids must be rejected, got {err:?}"
    );
    Ok(())
}

#[test]
fn cascade_rejects_1_port_network() {
    // A 1-port network must be rejected before it can hit upstream's panic.
    let one_port = ::touchstone::NetworkBuilder::new("one_port.s1p", 1)
        .point(
            1.0e9,
            ::touchstone::SMatrix {
                rank: 1,
                data: vec![vec![::touchstone::Complex { re: 0.5, im: -0.1 }]],
            },
        )
        .build()
        .expect("1-port network builds");
    let two_port = ::touchstone::NetworkBuilder::new("two_port.s2p", 2)
        .point(
            1.0e9,
            ::touchstone::SMatrix {
                rank: 2,
                data: vec![
                    vec![
                        ::touchstone::Complex { re: 0.1, im: 0.0 },
                        ::touchstone::Complex { re: 0.9, im: 0.0 },
                    ],
                    vec![
                        ::touchstone::Complex { re: 0.9, im: 0.0 },
                        ::touchstone::Complex { re: 0.1, im: 0.0 },
                    ],
                ],
            },
        )
        .build()
        .expect("2-port network builds");
    let err = mfs::touchstone::cascade(&one_port, &two_port).unwrap_err();
    assert!(matches!(err, mfs::MfsError::PreconditionViolation(_)));
}

#[test]
fn cascade_rejects_per_port_reference_impedance() {
    // A 2-port v2 file with a per-port `[Reference]` line parses into
    // ReferenceImpedance::PerPort while the scalar `z0` field holds the
    // first port's 50 Ω — so two such networks pass the `a.z0 != b.z0`
    // guard and would panic inside upstream `Network::cascade`. The wrapper
    // must reject them with a typed error instead.
    let text = "[Version] 2.1\n\
                # GHZ S RI R 50\n\
                [Number of Ports] 2\n\
                [Reference] 50 75\n\
                [Network Data]\n\
                1.0 0.1 0.0 0.9 0.0 0.9 0.0 0.1 0.0\n\
                2.0 0.2 0.0 0.8 0.0 0.8 0.0 0.2 0.0\n\
                [End]\n";
    let a = ::touchstone::Network::from_str("a.s2p", text).expect("parses per-port R network");
    assert!(
        matches!(a.reference_impedance(), ::touchstone::ReferenceImpedance::PerPort(_)),
        "test setup: expected PerPort impedance"
    );
    let err = mfs::touchstone::cascade(&a, &a).unwrap_err();
    assert!(
        matches!(err, mfs::MfsError::PreconditionViolation(_)),
        "per-port impedance must be rejected, got {err:?}"
    );
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
