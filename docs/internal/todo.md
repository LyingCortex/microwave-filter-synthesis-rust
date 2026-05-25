# MFS 滤波器综合库 — TODO List

## 高优先级（核心功能完善）

- [ ] **实现真正的 Wheel 拓扑变换** — 当前只是 Arrow 的别名，应实现环形耦合结构的专用旋转序列
- [ ] **改进 Placeholder 回退** — 要么移除它并让失败显式报错，要么实现一个基于 Cauchy 方法的备选综合路径
- [ ] **支持通带内传输零点**（`|zero| < 1`）— 放宽 `safe_sqrt_term` 的限制，处理复数传输零点
- [ ] **支持重极点留数展开** — 实现 Laurent 展开或多重极点分解

## 中优先级（工程质量）

- [ ] **增加高阶滤波器测试**（N=8, 10, 12）— 验证 Durand-Kerner 求根器在高阶时的数值稳定性
- [ ] **添加 HighPass / BandStop 频率映射** — 扩展 `freq` 模块支持更多滤波器类型
- [ ] **实现 Touchstone (.s2p) 导出** — 让响应数据可以被外部 EDA 工具读取
- [ ] **添加收敛失败的结构化错误恢复** — 当根求解器不收敛时，提供诊断信息而非泛化错误
- [ ] **`CouplingMatrix` 对称性强制** — 当前 `set_entry` 不自动保持对称，容易引入 bug

## 低优先级（生态扩展）

- [ ] **Python 绑定 (pyo3)** — 架构文档中 Phase 4 的目标
- [ ] **CLI 工具** — 命令行驱动的综合流程
- [ ] **Butterworth / Elliptic 近似族** — 扩展 `approx` 模块
- [ ] **多频带映射** — 支持双通带等复杂规格
- [ ] **性能基准测试** — 建立 criterion benchmarks，特别是矩阵旋转和响应求解热路径
- [ ] **文档完善** — 为公共 API 添加更多 doc-examples，准备 docs.rs 发布

## 代码质量改进

- [ ] **`coupling_matrix.rs` 拆分** — 该文件超过 500 行，建议将截面提取（triplet/quadruplet/trisection）移到独立子模块
- [ ] **统一错误语义** — `MfsError::Unsupported` 被用于太多不同场景（数值失败、未实现功能、无效输入），建议细分
- [ ] **减少 `unwrap_or_default()` 使用** — 矩阵访问中大量使用，可能掩盖越界 bug，考虑在 debug 模式下 panic
