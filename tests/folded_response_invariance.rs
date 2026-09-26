//! 回归测试: 用**真实合成滤波器**验证折叠保持响应。
//!
//! 背景: `tests/transform_invariance.rs` 里的折叠测试只用了一张全正、无自耦合、
//! 玩具 transversal 矩阵, 它掩盖了一个真实 bug——`CouplingMatrix::flip_sign`
//! 曾把对角线元素第三次变号, 导致偶数阶非对称滤波器折叠后自耦合符号错误、
//! 传输零点丢失。此处用 Cameron 书里的黄金算例(真实合成滤波器)覆盖该缺陷,
//! 防止回归。

use mfs::error::Result;
use mfs::freq::FrequencyGrid;
use mfs::prelude::FilterDesign;
use mfs::response::ResponseSolver;
use mfs::transform::{transform_matrix, TopologyKind};
use mfs::verify::{compare_responses, ResponseTolerance};

/// 合成一个 generalized-Chebyshev 滤波器, 折叠, 并断言响应保持不变。
fn assert_folded_preserves_response(
    name: &str,
    order: usize,
    rl_db: f64,
    zeros: &[f64],
) -> Result<()> {
    let design = FilterDesign::prototype(order, rl_db).zeros(zeros.to_vec()).synthesize()?;
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

/// 偶数阶、非对称 TZ 分布——正是曾丢失传输零点的场景。
#[test]
fn folded_preserves_response_order4_asymmetric_2tz() -> Result<()> {
    assert_folded_preserves_response("4-2 (Cameron Fig 8.23b)", 4, 22.0, &[1.3217, 1.8082])
}

#[test]
fn folded_preserves_response_order6_asymmetric_3tz() -> Result<()> {
    assert_folded_preserves_response("6-3 (Cameron Fig 8.32)", 6, 22.0, &[1.0900, 1.1691, 1.5057])
}

#[test]
fn folded_preserves_response_order8_asymmetric_3tz() -> Result<()> {
    assert_folded_preserves_response("8-3 (Cameron Fig 9.22a)", 8, 23.0, &[-1.3553, 1.1093, 1.2180])
}

/// 偶数阶、对称 TZ 分布。
#[test]
fn folded_preserves_response_order6_symmetric_2tz() -> Result<()> {
    assert_folded_preserves_response("6-2 (Cameron Fig 9.21)", 6, 23.0, &[1.5522, -1.5522])
}

/// 奇数阶、非对称 TZ 分布。
#[test]
fn folded_preserves_response_order7_asymmetric_1tz() -> Result<()> {
    assert_folded_preserves_response("7-1-2 (Cameron Fig 8.13)", 7, 23.0, &[1.2576])
}

/// 全规范 (nfz == order): 4 阶 4 个 TZ。此前 synthesis 直接报维度错误,
/// 现已可合成并正确折叠。
#[test]
fn folded_preserves_response_order4_fully_canonical() -> Result<()> {
    assert_folded_preserves_response("4-4 fully canonical (Cameron Fig 8.21)", 4, 22.0,
        &[3.7431, 1.8051, -1.5699, -6.1910])
}

/// 全规范偶数阶 6 阶 6 个 TZ。
#[test]
fn folded_preserves_response_order6_fully_canonical() -> Result<()> {
    assert_folded_preserves_response("6-6 fully canonical", 6, 22.0,
        &[2.1, 1.6, 1.2, -1.1, -1.5, -1.9])
}

/// 高阶、偶数阶、非对称(压力测试)。
#[test]
fn folded_preserves_response_order10_asymmetric_4tz() -> Result<()> {
    assert_folded_preserves_response("10-4 asymmetric", 10, 23.0, &[1.2, 1.5, -1.8, -2.0])
}

/// 全极(无 TZ)——纯主链, 最基础的折叠。
#[test]
fn folded_preserves_response_order4_allpole() -> Result<()> {
    assert_folded_preserves_response("4 all-pole", 4, 20.0, &[])
}

/// 奇数阶全极。
#[test]
fn folded_preserves_response_order5_allpole() -> Result<()> {
    assert_folded_preserves_response("5 all-pole", 5, 20.0, &[])
}