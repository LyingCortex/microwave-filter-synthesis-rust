# CT/CQ 拓扑提取算法设计 (mfs Rust实现)

## 概述

在mfs中实现经典Givens旋转法的CT（Cascaded Trisection）和CQ（Cascaded Quadruplet）拓扑变换。

**输入**：Transversal或Arrow耦合矩阵 + 传输零点列表 + 目标位置  
**输出**：CT或CQ形式的耦合矩阵（特定稀疏模式）

> ✅ **状态 (2026-10-06, 实测已修复 + order-8 两 CQ 平铺已验证)**：CT/CQ 均改走**边缘-capable** 的已验证原语——`to_cascaded_trisections` 用 `extract_triplet`（`validate_triplet_center` 允许 center=1，产生 source 端连接的 S-1-2 trisection），`to_cascaded_quadruplets` 用 `extract_cascade_quadruplets`（与 Python `cq` 绑定同一条路径）。**CQ 平铺规则（关键）**：quadruplet 必须**互不重叠、自 source 侧起**平铺——第 *i* 个零点对占用谐振器 `4i-3..=4i`，因此需要 `order >= 4 × num_pairs`，否则返回 `PreconditionViolation`（密集情形请用 `to_cascaded_trisections`，其 3-谐振器 trisection 更省）。实测：order-4 单对 → 干净 quadruplet；**order-6 单对（2 零点）→ 干净单 quadruplet**（`M[1,4]` 交叉耦合存在、主链 1..=6 连通、source→res1/load→res6、无其他 2-skip、max|ΔS21| 2e-14）；**order-8 两对（4 零点）→ 复现已验证的 `CascadeQuartets` 布局** `[[1,3],[1,4],[5,7],[5,8]]`（source 侧非重叠平铺，max|ΔS21| 1.4e-13，与 CT 结构可区分 dM>0.05）。注意：**order-6 配两对（4 零点）在拓扑上不可能干净平铺**（两个 4-谐振器 quadruplet 无法互不重叠地放进 6 个谐振器），旧版尾部 peel 循环留下的 ~1e-2 skip 耦合不是 bug 而是这种越约束情形的固有残差。CT 按 Fig. 21 (order 4, 22 dB, TZ 1.8082/1.3217) 实测复现。回归门：`ct_synthesis_reproduces_fig21_cameron2011`、`cq_synthesis_single_quadruplet_structure`、`cq_synthesis_order6_single_quadruplet`、`cq_synthesis_order8_two_quadruplets`、`cq_synthesis_rejects_too_small_order`、`ct_synthesis_produces_valid_matrix`。

---

## 算法详述

### 1. Trisection提取 (单个零点)

**输入**: Arrow矩阵M (N+2×N+2), 传输零点Ωz, 目标位置pos  
**输出**: 含一个trisection交叉耦合的矩阵

**步骤**:

```
Step 1: 计算旋转角度θ
  - 末端旋转消除Arrow最后一列
  - θ = atan(M[N, N-1] / (Ωz + M[N, N]))
  - 旋转pivot: (N-1, N)，即最后两个谐振器

Step 2: 应用Givens相似变换
  - R[N-1,N-1] = cos(θ),  R[N-1,N] = -sin(θ)
  - R[N,N-1]   = sin(θ),  R[N,N]   = cos(θ)  
  - M' = R × M × R^T

Step 3: Pull操作（移动trisection到目标位置）
  - 需要 (N-1-pos) 次旋转
  - 每次消零: 将交叉耦合从位置k移到k-1
  - 旋转消除: M'[k, k-2] 元素，pivot=(k-2, k-1)
  - 方向: "row" 旋转，θ = atan(M[k-2, col] / M[k-1, col])
```

**结果矩阵模式** (N=6, 零点在pos=2):
```
     S   1   2   3   4   L
S [  0  Ms1  0   0   0   0  ]
1 [ Ms1  0  M12  0   0   0  ]
2 [  0  M12  0  M23 M24  0  ]  ← M24是交叉耦合(跳过节点3)
3 [  0   0  M23  0  M34  0  ]
4 [  0   0  M24 M34  0  MnL ]
L [  0   0   0   0  MnL  0  ]
```

### 2. Cascade Trisection (多个零点)

**输入**: Arrow矩阵M, 零点列表[Ω1, Ω2, ...], 位置列表[(p1_from, p1_to), ...]  
**输出**: 含多个trisection的矩阵

**步骤**:
```
for each (Ωz, pos) in zip(zeros, positions):
    M = extract_trisection(M, Ωz, pos)
```

每次提取后Arrow结构缩短一列（已提取的trisection不再参与后续旋转）。

### 3. Cascade Quadruplet (零点成对, source 侧非重叠平铺)

**输入**: Arrow矩阵M, 零点对[(Ω1,Ω2), (Ω3,Ω4), ...], 位置列表  
**输出**: 含quadruplet块的矩阵

**平铺规则**（实测已验证）：每个 quadruplet 占 4 个谐振器，**互不重叠、自 source 侧起**平铺——第 *i* 个零点对（1-indexed）占用谐振器 `4i-3..=4i`，其交叉耦合为 `M[4i-3, 4i]`。因此必须满足 `order >= 4 × num_pairs`；否则（如 order-6 配两对）两个 4-谐振器块无法互不重叠地放入链中，旧版尾部 peel 会留下 ~1e-2 skip 耦合（拓扑不可能，非 bug）。

**步骤**:
```
# 前置校验：order >= 4 * num_pairs，否则拒绝（密集情形改用 CT）
for i in 1..=num_pairs:
    s = 4*(i-1) + 1              # 第 i 个 quadruplet 的 source 侧起点
    # trisection partner: span 2, [s, s+2]
    # quadruplet 对角: span 3, [s, s+3]
    z_locs.push([s, s+2])
    z_locs.push([s, s+3])
M = extract_cascade_quadruplets(M, zeros, z_locs)
# 内部按位置从尾向 source 逐个 extract_quadruplet 合并（common_res=4）
```

**Quadruplet 结果模式** (4 个谐振器块, 交叉耦合 `M14` 存在):
```
[  0   M12   0   M14 ]
[ M12   *   M23   *  ]    * = 谐振器对角(detuning, 一般非零)
[  0   M23   0   M34 ]
[ M14   *   M34   0  ]
```
`M14` 是 quadruplet 的交叉耦合，实现一对传输零点。注意：order-4 单对情形谐振器对角被 peel 清零，但 order-8 两对情形对角保留 detuning（如 M[2,2]=-0.021），这是已验证参考布局的固有特征，不作为"必须清零"的门槛。

---

## Rust实现架构

### 文件位置
`src/transform/sections.rs` (已存在，追加新函数)

### 公共API

```rust
/// 从Arrow矩阵提取单个trisection (CT)
pub fn extract_trisection_at(
    matrix: &CouplingMatrix,
    transmission_zero: f64,
    target_position: usize,  // 1-indexed resonator position
) -> Result<CouplingMatrix>;

/// 从Arrow矩阵提取多个trisection (完整CT)
pub fn extract_cascade_trisections(
    matrix: &CouplingMatrix,
    zeros: &[f64],
    positions: &[[usize; 2]],  // [from_node, to_node] per zero
) -> Result<CouplingMatrix>;

/// 从Arrow矩阵提取cascade quadruplets (CQ)
pub fn extract_cascade_quadruplets(
    matrix: &CouplingMatrix,
    zeros: &[f64],
    positions: &[[usize; 2]],
) -> Result<CouplingMatrix>;
```

### 内部辅助函数

```rust
/// 末端旋转锁定零点
fn trisection_end_rotation(matrix: &mut CouplingMatrix, tz: f64);

/// Pull操作：将交叉耦合从末端移到目标位置
fn pull_trisection(matrix: &mut CouplingMatrix, from: usize, to: usize);

/// 合并两个相邻trisection为quadruplet
fn merge_to_quadruplet(matrix: &mut CouplingMatrix, position: usize);
```

### Python绑定扩展

```rust
// FilterDesign新增方法
impl PyFilterDesign {
    fn ct(&self, positions: Vec<[usize; 2]>) -> PyResult<PyObject>;
    fn cq(&self, positions: Vec<[usize; 2]>) -> PyResult<PyObject>;
}
```

---

## 验证计划

1. 4阶2零点CT vs Python `TrisectionExtractList` → folded CM一致
2. 6阶3零点CT vs Python → S参数一致
3. 8阶4零点CQ vs Python `CascadeQuartets` → folded CM一致
4. 对称零点验证矩阵对称性保持
5. S参数验证：CT/CQ矩阵产生的响应与transversal矩阵相同

---

## 依赖
- `apply_givens_similarity_inplace` (已有)
- `to_arrow()` (已有)
- `safe_angle()` (已有)
- 无新外部依赖

## 预计代码量
~200行Rust

---

## 第二阶段：Householder变换框架 (后续实现)

### 核心优势
- 数值稳定性远高于Givens旋转序列
- 最多N步完成任意拓扑（vs Givens需要O(N²)次旋转）
- 统一框架覆盖CT/CQ/BOX/EB/任意自定义拓扑

### 数学基础

```
H = I - 2·u·uᵀ/(uᵀ·u)

矩阵重构: M_N = (∏ Hᵢ)ᵀ · M₀ · (∏ Hᵢ)
```

u的格式: `uᵀ = [0, ũᵀ, 0]` → 保持源/负载不变

### 算法A：解析重构（Folded/Arrow）

```
输入: 横向矩阵 M₀
for i = 1 to N-1:
    x = 当前待消列向量
    z = 目标单位向量方向
    u = x - |x|·z
    H_i = I - 2·u·uᵀ/(uᵀ·u)
    M = Hᵢᵀ · M · Hᵢ
输出: 标准拓扑矩阵
```

### 算法B：指定传输零点 (CT/CQ直接分配)

```
Step 1: h₁ = (Mc + ωz·I)⁻¹·a₁ / ‖(Mc + ωz·I)⁻¹·a₁‖
  - Mc: 核心矩阵区 (N×N内部)
  - a₁: 源端耦合向量
  - ωz: 目标传输零点

Step 2: u = e₁ - h₁

Step 3: M₁ = H₁·M₀·H₁ → 在指定位置建立线性相关性

后续HT步骤直接形成CT/CQ块。
```

### 算法C：优化重构 (EB等复杂拓扑)

```
目标函数: L(U) = ½·‖T ∘ M_N(U)‖²_F
  - T: 掩码矩阵(应为零的位置=1)
  - ∘: Hadamard积

梯度: ∂L/∂uᵢ = 8·[(uᵢᵀRᵢuᵢ)/(uᵢᵀuᵢ)²·uᵢ - Rᵢuᵢ/(uᵢᵀuᵢ)]
  - Rᵢ: 中间矩阵(块结构加速计算)

优化器: L-BFGS 或 Adam, 收敛条件 L < 1e-8
初始化: 拉丁超立方采样(LHS)覆盖高维空间
```

### Rust实现规划

```rust
// src/transform/householder.rs (新文件)

/// 解析HT：Transversal → Folded/Arrow (N-1步)
pub fn householder_to_folded(matrix: &CouplingMatrix) -> Result<CouplingMatrix>;
pub fn householder_to_arrow(matrix: &CouplingMatrix) -> Result<CouplingMatrix>;

/// 零点指定HT：直接分配TZ到CT/CQ位置
pub fn householder_assign_zero(
    matrix: &CouplingMatrix,
    zero: f64,
    target_position: usize,
) -> Result<CouplingMatrix>;

/// 优化HT：任意拓扑mask → 梯度优化
pub fn householder_optimize_topology(
    matrix: &CouplingMatrix,
    topology_mask: &[Vec<bool>],  // true=应为零
    zeros: &[f64],                // 需要保持的TZ
    max_iter: usize,
) -> Result<CouplingMatrix>;

/// WY加速：累积多次HT为单次矩阵运算
fn wy_accumulate(householder_vectors: &[Vec<f64>]) -> (DMatrix<f64>, DMatrix<f64>);
```

### 依赖
- `nalgebra` 或 `faer`: 矩阵求逆(算法B的(Mc+ωzI)⁻¹)
- 优化器: 自研L-BFGS或复用`optimize.rs`

### 预计代码量
~400行Rust

### 优先级
在Givens法CT/CQ验证通过后实现。对于标准CT/CQ Givens法已足够，
HT主要用于：
1. Extended Box等无解析路径的拓扑
2. 需要更好数值稳定性的高阶(N>12)场景
3. 多解搜索（配合LHS初始化）
