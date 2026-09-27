//! 回归测试: 用**真实合成滤波器**验证折叠保持响应。
//!
//! 背景: `tests/transform_invariance.rs` 里的折叠测试只用了一张全正、无自耦合、
//! 玩具 transversal 矩阵, 它掩盖了一个真实 bug——`CouplingMatrix::flip_sign`
//! 曾把对角线元素第三次变号, 导致偶数阶非对称滤波器折叠后自耦合符号错误、
//! 传输零点丢失。此处用 Cameron 书里的黄金算例(真实合成滤波器)覆盖该缺陷,
//! 防止回归。
//!
//! 算例规格 **不再手抄**: 全部来自跨语言基准标准的 vendored 副本
//! (`bench-data/cases/twoport.json`, 由 `tests/benchmark_standard.rs` 做漂移校验)。
//! 每个测试直接引用标准的全局 id (`twoport_*`), 所以标准增删算例时这里不会漂移。

#[path = "bench_support/mod.rs"]
mod bench_support;

use bench_support::data;
use mfs::error::Result;
use mfs::freq::FrequencyGrid;
use mfs::prelude::FilterDesign;
use mfs::response::ResponseSolver;
use mfs::transform::{transform_matrix, TopologyKind};
use mfs::verify::{compare_responses, ResponseTolerance};

/// 按标准的全局 id 取出算例规格 (order / rl / normalized TZ)。
fn standard_case(qualid: &str) -> (usize, f64, Vec<f64>) {
    let case = data::case(qualid).unwrap_or_else(|e| panic!("standard case {qualid}: {e}"));
    let order = bench_support::gates::spec_usize(&case.spec, "order").expect("spec.order");
    let rl = bench_support::gates::spec_f64(&case.spec, "rl").expect("spec.rl");
    let zeros = bench_support::gates::spec_f64_vec(&case.spec, "zeros").expect("spec.zeros");
    (order, rl, zeros)
}

/// 合成一个 generalized-Chebyshev 滤波器, 折叠, 并断言响应保持不变。
fn assert_folded_preserves_response(
    name: &str,
    order: usize,
    rl_db: f64,
    zeros: &[f64],
) -> Result<()> {
    let design = FilterDesign::prototype(order, rl_db)
        .zeros(zeros.to_vec())
        .synthesize()?;
    let transversal = design.matrix().clone();
    let folded = transform_matrix(&transversal, TopologyKind::Folded)?.matrix;

    // 用更细的网格覆盖到阻带深处的传输零点(那里对符号错误最敏感)。
    let grid = FrequencyGrid::linspace(-4.0, 4.0, 801)?;
    let solver = ResponseSolver {};
    let baseline = solver.evaluate_normalized(&transversal, &grid)?;
    let transformed = solver.evaluate_normalized(&folded, &grid)?;
    let comparison = compare_responses(&baseline, &transformed)?;

    assert!(
        comparison.passes(ResponseTolerance::strict()),
        "{name}: folded response changed from transversal: {comparison:?}"
    );
    Ok(())
}

/// 等同 `assert_folded_preserves_response`, 但规格取自标准。
fn assert_standard_case(qualid: &str) -> Result<()> {
    let (order, rl, zeros) = standard_case(qualid);
    assert_folded_preserves_response(qualid, order, rl, &zeros)
}

/// 偶数阶、非对称 TZ 分布——正是曾丢失传输零点的场景。
#[test]
fn folded_preserves_response_order4_asymmetric_2tz() -> Result<()> {
    // 4-2 (Cameron Fig 8.23b) — spec from the standard, not hand-copied.
    assert_standard_case("twoport_4_2_unity")
}

#[test]
fn folded_preserves_response_order6_asymmetric_3tz() -> Result<()> {
    assert_standard_case("twoport_6_3")
}

#[test]
fn folded_preserves_response_order8_asymmetric_3tz() -> Result<()> {
    assert_standard_case("twoport_8_3")
}

/// 偶数阶、对称 TZ 分布。
#[test]
fn folded_preserves_response_order6_symmetric_2tz() -> Result<()> {
    assert_standard_case("twoport_6_2_culdesac")
}

/// 奇数阶、非对称 TZ 分布。
#[test]
fn folded_preserves_response_order7_asymmetric_1tz() -> Result<()> {
    // 7-1-2 原本只存在于本文件的手抄常量; 现由标准 (twoport_7_1_2) 承载。
    assert_standard_case("twoport_7_1_2")
}

/// 全规范 (nfz == order): 4 阶 4 个 TZ。此前 synthesis 直接报维度错误,
/// 现已可合成并正确折叠。
#[test]
fn folded_preserves_response_order4_fully_canonical() -> Result<()> {
    assert_standard_case("twoport_4_4")
}

/// 全规范偶数阶 6 阶 6 个 TZ —— 标准没有对应的印刷算例, 保留为程序生成的补充。
#[test]
fn folded_preserves_response_order6_fully_canonical() -> Result<()> {
    assert_folded_preserves_response(
        "6-6 fully canonical",
        6,
        22.0,
        [2.1, 1.6, 1.2, -1.1, -1.5, -1.9].as_slice(),
    )
}

/// 高阶、偶数阶、非对称(压力测试)—— 标准无对应算例, 保留为补充。
#[test]
fn folded_preserves_response_order10_asymmetric_4tz() -> Result<()> {
    assert_folded_preserves_response(
        "10-4 asymmetric",
        10,
        23.0,
        [1.2, 1.5, -1.8, -2.0].as_slice(),
    )
}

/// 全极(无 TZ)——纯主链, 最基础的折叠。
#[test]
fn folded_preserves_response_order4_allpole() -> Result<()> {
    assert_folded_preserves_response("4 all-pole", 4, 20.0, [].as_slice())
}

/// 奇数阶全极。
#[test]
fn folded_preserves_response_order5_allpole() -> Result<()> {
    assert_folded_preserves_response("5 all-pole", 5, 20.0, [].as_slice())
}
