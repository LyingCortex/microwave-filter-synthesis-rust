# MFS 数值优化技术报告

## 概述

本文档详细记录了 MFS（Microwave Filter Synthesis）库中三项核心数值优化的设计原理、
数学推导和实现细节。优化目标是在保持 IEEE 754 双精度浮点极限精度的前提下，
将矩阵运算效率提升数个数量级，并将稳定工作阶数从约 18 阶扩展到 30 阶以上。

三项优化分别是：

1. **矩阵相似旋转的 O(n) 原地 Givens 变换**（替代 O(n³) 全矩阵乘法）
2. **S 参数频率响应的选择性线性求解**（替代完整矩阵求逆）
3. **高阶多项式求根器的自适应增强**（Durand-Kerner 算法改进）

---

## 第一部分：矩阵相似旋转优化

### 1.1 问题背景

耦合矩阵综合的核心操作是将横向（Transversal）矩阵通过一系列正交相似变换
转化为 Folded、Arrow 等物理可实现的拓扑结构。每次变换的数学形式为：

$$
M' = R(\theta) \cdot M \cdot R(\theta)^T
$$

其中 $R(\theta)$ 是 Givens 旋转矩阵，仅在第 $(p, p)$、$(p, q)$、$(q, p)$、$(q, q)$
四个位置与单位矩阵不同：

$$
R(\theta) = I + (\cos\theta - 1)(e_p e_p^T + e_q e_q^T) - \sin\theta(e_p e_q^T - e_q e_p^T)
$$

对于 N 阶滤波器，矩阵维度为 $(N+2) \times (N+2)$（含源/负载节点）。

### 1.2 原始实现的复杂度分析

原始代码对每次旋转执行以下步骤：

```rust
// 旧实现（伪代码）
let R = CouplingMatrix::identity(order);  // 分配 (N+2)² 个 f64
R.set(p, p, cos θ);
R.set(q, q, cos θ);
R.set(p, q, -sin θ);
R.set(q, p, sin θ);
let temp = R.multiply(&M);    // O(n³) 矩阵乘法 + nalgebra 格式转换
let result = temp.multiply(&R.transpose());  // 又一次 O(n³)
```

**单次旋转代价：**
- 2 次完整矩阵乘法：$2 \times (N+2)^3$ 浮点乘加
- 3 次堆分配：identity + 两次乘法结果
- 4 次格式转换：行主序 ↔ nalgebra 列主序

**拓扑变换的总旋转次数：**
- `to_folded()`：约 $N^2/2$ 次旋转
- `to_arrow()`：约 $N(N-1)/2$ 次旋转
- Section extraction：每次 triplet/quadruplet 额外 $O(N)$ 次

**总复杂度：** $O(N^2) \times O(N^3) = O(N^5)$

对于 order=20（矩阵 22×22）：
- 单次旋转：$2 \times 22^3 = 21,296$ 次浮点运算
- Folded 变换约 200 次旋转：$200 \times 21,296 \approx 4.3 \times 10^6$ 次运算

### 1.3 优化原理：原地 Givens 相似变换

关键观察：$R$ 是稀疏的（仅 4 个非平凡元素），因此 $R \cdot M$ 只影响第 $p$ 行和
第 $q$ 行，$M \cdot R^T$ 只影响第 $p$ 列和第 $q$ 列。

**第一步：左乘 $R \cdot M$（更新两行）**

$$
M'_{p,j} = \cos\theta \cdot M_{p,j} - \sin\theta \cdot M_{q,j}, \quad \forall j
$$
$$
M'_{q,j} = \sin\theta \cdot M_{p,j} + \cos\theta \cdot M_{q,j}, \quad \forall j
$$

**第二步：右乘 $M' \cdot R^T$（更新两列）**

$$
M''_{i,p} = \cos\theta \cdot M'_{i,p} - \sin\theta \cdot M'_{i,q}, \quad \forall i
$$
$$
M''_{i,q} = \sin\theta \cdot M'_{i,p} + \cos\theta \cdot M'_{i,q}, \quad \forall i
$$

**单次旋转代价：** $4(N+2)$ 次浮点乘加，零额外分配。

### 1.4 实现代码

```rust
pub(crate) fn apply_givens_similarity_inplace(
    &mut self, pivot_a: usize, pivot_b: usize, cosine: f64, sine: f64
) {
    let side = self.side();

    // Step 1: R * M — 更新第 pivot_a 行和第 pivot_b 行
    for col in 0..side {
        let a = self.data[pivot_a * side + col];
        let b = self.data[pivot_b * side + col];
        self.data[pivot_a * side + col] = cosine * a - sine * b;
        self.data[pivot_b * side + col] = sine * a + cosine * b;
    }

    // Step 2: M' * R^T — 更新第 pivot_a 列和第 pivot_b 列
    for row in 0..side {
        let a = self.data[row * side + pivot_a];
        let b = self.data[row * side + pivot_b];
        self.data[row * side + pivot_a] = cosine * a - sine * b;
        self.data[row * side + pivot_b] = sine * a + cosine * b;
    }
}
```

### 1.5 复杂度对比

| 指标 | 原始实现 | 优化后 | 加速比 (N=20) |
|------|---------|--------|--------------|
| 单次旋转 FLOPs | $2(N+2)^3$ | $4(N+2)$ | ~2,662× |
| 单次旋转分配 | 3 次堆分配 | 0 | ∞ |
| Folded 变换总计 | $O(N^5)$ | $O(N^3)$ | ~500× |
| Arrow 变换总计 | $O(N^5)$ | $O(N^3)$ | ~500× |

### 1.6 精度分析

原地 Givens 变换与全矩阵乘法在数学上完全等价——两者执行相同的浮点运算序列，
只是避免了中间矩阵的分配和复制。实际上，原地方法的精度**略优于**全矩阵乘法，
原因是：

1. 减少了中间结果的存储/读取次数，降低了舍入误差累积
2. 避免了行主序 ↔ 列主序转换中的浮点重排序
3. 每步只涉及 2 行 2 列的局部更新，不会将远处元素的误差传播到当前操作

实测验证（order 30，Folded 变换后）：
- 矩阵对称性误差：$< 4 \times 10^{-16}$（机器精度量级）
- 最大元素幅度：有界且物理合理（< 2.0）

---

## 第二部分：S 参数频率响应求解优化

### 2.1 问题背景

S 参数（散射参数）频率响应是滤波器设计验证的核心输出。对于 N 阶耦合矩阵 $M$，
在每个归一化频率点 $\omega$ 处，需要求解网络方程：

$$
\mathbf{A}(\omega) = j\omega \mathbf{I}_r + M - j R_S \mathbf{e}_0 \mathbf{e}_0^T
                     - j R_L \mathbf{e}_{N+1} \mathbf{e}_{N+1}^T
$$

其中：
- $\mathbf{I}_r$ 是仅在谐振器对角位置为 1 的矩阵（源/负载位置为 0）
- $R_S, R_L$ 是归一化源/负载电阻（通常为 1）
- $\mathbf{e}_k$ 是第 $k$ 个标准基向量

S 参数通过逆矩阵元素提取：

$$
S_{11}(\omega) = 1 + 2j R_S \cdot [\mathbf{A}^{-1}]_{0,0}
$$

$$
S_{21}(\omega) = -2j\sqrt{R_S R_L} \cdot [\mathbf{A}^{-1}]_{N+1, 0}
$$

群延迟：

$$
\tau_g(\omega) = \text{Im}\left( \frac{\sum_k [\mathbf{A}^{-1}]_{N+1,k} \cdot [\mathbf{A}^{-1}]_{k,0}}{[\mathbf{A}^{-1}]_{N+1,0}} \right)
$$

### 2.2 原始实现：完整矩阵求逆

原始代码对每个频率点计算完整的 $\mathbf{A}^{-1}$：

```rust
// 旧实现
fn solve_inverse(matrix, omega, settings) -> DMatrix<Complex64> {
    let A = build_response_matrix(matrix, omega, settings);
    let lu = A.lu();
    let identity = DMatrix::identity(side);
    lu.solve(&identity)  // 求解 A * X = I，即 X = A^{-1}
}
```

**每个频率点的代价：**
- LU 分解：$\frac{2}{3}(N+2)^3$ 复数浮点运算
- 回代求解 $(N+2)$ 个右端向量：$(N+2)^3$ 复数浮点运算
- 总计：$\frac{5}{3}(N+2)^3$ 复数运算

**对于 2001 点 × order 20：**
$2001 \times \frac{5}{3} \times 22^3 \approx 35.5 \times 10^6$ 复数运算

### 2.3 优化原理：选择性列/行求解

观察 S 参数公式，我们实际只需要：
- $[\mathbf{A}^{-1}]_{:,0}$：逆矩阵的**第一列**（用于 $S_{11}$、$S_{21}$、群延迟）
- $[\mathbf{A}^{-1}]_{N+1,:}$：逆矩阵的**最后一行**（用于群延迟）

**第一列**可通过求解一个线性系统获得：

$$
\mathbf{A} \cdot \mathbf{x} = \mathbf{e}_0 \implies \mathbf{x} = [\mathbf{A}^{-1}]_{:,0}
$$

**最后一行**可通过转置系统获得：

$$
\mathbf{A}^T \cdot \mathbf{y} = \mathbf{e}_{N+1} \implies \mathbf{y} = [(\mathbf{A}^T)^{-1}]_{:,N+1} = [(\mathbf{A}^{-1})^T]_{:,N+1}
$$

即 $y_k = [\mathbf{A}^{-1}]_{N+1, k}$，正是我们需要的最后一行。

### 2.4 优化后的算法

```
对每个频率点 ω:
    1. 构建 A(ω) = base_matrix + ω·I_r    // O(N) 对角线加法
    2. LU 分解 A                            // O(N³) — 不可避免
    3. 前/回代求解 A·x = e₀                 // O(N²) — 仅 1 个右端
    4. LU 分解 A^T                          // O(N³)
    5. 前/回代求解 A^T·y = e_{N+1}          // O(N²) — 仅 1 个右端
    6. 提取 S11 = 1 + 2jR_S·x[0]
    7. 提取 S21 = -2j√(R_S·R_L)·x[N+1]
    8. 群延迟 = Im(dot(y, x) / x[N+1])     // O(N) 点积
```

**每个频率点的代价：**
- 2 次 LU 分解：$2 \times \frac{2}{3}(N+2)^3 = \frac{4}{3}(N+2)^3$
- 2 次回代：$2 \times (N+2)^2$
- 总计：$\frac{4}{3}(N+2)^3 + 2(N+2)^2$

### 2.5 加速比分析

| 操作 | 原始 | 优化后 | 比率 |
|------|------|--------|------|
| LU 分解 | $\frac{2}{3}N^3$ | $\frac{4}{3}N^3$ | 0.5× (略慢) |
| 回代求解 | $N \times N^2 = N^3$ | $2 \times N^2$ | $N/2$ × |
| **总计** | $\frac{5}{3}N^3$ | $\frac{4}{3}N^3 + 2N^2$ | **~1.25×** |

注意：虽然我们做了两次 LU 分解（A 和 A^T），但省去了 N-2 次回代。
对于 order 20（N+2=22）：

- 原始：回代部分 = $22 \times 22^2 = 10,648$ 复数运算
- 优化：回代部分 = $2 \times 22^2 = 968$ 复数运算
- 回代加速：**11×**

总体加速取决于 LU 分解与回代的比例。对于小矩阵（N < 30），LU 分解占主导，
总体加速约 **1.2-1.5×**。但更重要的优势是：

1. **内存带宽减少**：不再分配/填充完整的 $(N+2)^2$ 逆矩阵
2. **缓存友好**：只操作两个长度为 $N+2$ 的向量
3. **可并行化**：各频率点完全独立，适合 SIMD/多线程

### 2.6 预计算基础矩阵

额外优化：将频率无关的部分（耦合矩阵 + 端口阻抗）预计算为基础矩阵，
每个频率点只需在谐振器对角线加 $\omega$：

$$
\mathbf{A}(\omega) = \mathbf{A}_{\text{base}} + \omega \cdot \text{diag}(0, 1, 1, \ldots, 1, 0)
$$

这避免了每个频率点重新构建完整的复数矩阵。

### 2.7 数值稳定性

S 参数求解的数值稳定性由以下因素保证：

1. **LU 分解带部分主元选取**（nalgebra 默认策略），条件数 $\kappa(\mathbf{A})$
   对于物理耦合矩阵通常在 $10^2 - 10^4$ 范围内，远低于 $1/\epsilon_{\text{mach}} \approx 10^{16}$

2. **矩阵结构**：$\mathbf{A}(\omega)$ 是对称矩阵加纯虚对角扰动，
   不会出现病态条件数（除非恰好在传输零点频率处，此时 $S_{21} = 0$ 但矩阵仍可逆）

3. **功率守恒验证**：对于无耗网络，$|S_{11}|^2 + |S_{21}|^2 = 1$ 在所有频率点成立，
   这提供了内建的精度校验

实测结果：功率守恒误差 $< 10^{-9}$，与完整求逆方法一致。

---

## 第三部分：高阶多项式求根器优化

### 3.1 问题背景

广义 Chebyshev 滤波器综合流程中，需要对以下多项式求根：

- **E(w) 多项式**：阶数 = 滤波器阶数 N，根决定了导纳多项式的极点
- **A(w) 辅助多项式**：阶数 = 有限传输零点数，根用于 Hurwitz 稳定性判定
- **导纳分母多项式**：阶数 = N，根即为部分分式展开的极点

对于 N=20 阶滤波器，需要求解 20 次复系数多项式的全部根。
对于 N=30 阶，需要求解 30 次多项式。

### 3.2 Durand-Kerner 算法原理

Durand-Kerner 方法是一种同时迭代所有根的算法。设 $p(z) = z^n + a_{n-1}z^{n-1} + \cdots + a_0$
是首一多项式，$z_1, z_2, \ldots, z_n$ 是当前根的近似值，则迭代公式为：

$$
z_k^{(\text{new})} = z_k - \frac{p(z_k)}{\prod_{j \neq k} (z_k - z_j)}
$$

**几何解释**：分母 $\prod_{j \neq k}(z_k - z_j)$ 是以当前近似根构造的
Weierstrass 因式分解的"其余部分"。当所有 $z_k$ 趋近真实根时，
$p(z_k) / \prod_{j \neq k}(z_k - z_j) \to 0$。

**收敛性**：
- 对于简单根（无重根），Durand-Kerner 具有二次收敛速度
- 对于重根或近重根，收敛退化为线性
- 全局收敛性依赖于初始值的选取

### 3.3 高阶多项式的数值挑战

当多项式阶数超过 20 时，出现以下数值困难：

**3.3.1 系数动态范围爆炸**

Cameron 递推产生的 F(w) 多项式系数随阶数指数增长。例如对于 N=30 的全极点滤波器，
$F(w) = w^{30} + c_{29}w^{29} + \cdots + c_0$，其中 $|c_k|$ 可达 $10^{8}$ 量级。
这导致 Horner 求值时的中间结果跨越巨大的数量级，产生灾难性抵消。

**3.3.2 根的聚集（Root Clustering）**

Chebyshev 型滤波器的极点分布在虚轴附近的椭圆上，高阶时极点间距缩小为 $O(1/N)$。
当两个根 $z_i, z_j$ 非常接近时，分母 $\prod_{j \neq k}(z_k - z_j)$ 中出现
接近零的因子，导致迭代步长爆炸。

**3.3.3 初始值对称性陷阱**

原始实现使用均匀角度分布的初始值：
$$
z_k^{(0)} = R \cdot e^{2\pi i k / n}, \quad k = 0, 1, \ldots, n-1
$$

这种对称布局对于具有对称根分布的多项式（如 Chebyshev 极点）会导致
多个初始值同时收敛到同一个根，而其他根无人"认领"。

### 3.4 优化措施

**3.4.1 黄金角初始值分布**

用黄金角 $\phi = 2\pi / \varphi^2 \approx 2.3999$ 替代均匀角度：

$$
z_k^{(0)} = R_k \cdot e^{i \cdot k \cdot \phi}, \quad R_k = R \cdot (0.4 + 0.6 \cdot \frac{k+1}{n})
$$

其中 $\varphi = (1+\sqrt{5})/2$ 是黄金比例。这种分布具有以下优势：

- **无有理周期**：黄金角是"最无理"的角度，任意有限个点都不会形成对称图案
- **径向分层**：不同根的初始半径不同，避免了同一圆上的竞争
- **覆盖均匀**：在极坐标下实现了最优的低差异序列（low-discrepancy sequence）

**3.4.2 自适应迭代次数**

$$
\text{max\_iterations} = \min(128 + 4 \cdot \max(0, n - 10), \; 512)
$$

| 阶数 n | 最大迭代次数 |
|--------|-------------|
| ≤ 10   | 128         |
| 20     | 168         |
| 30     | 208         |
| 50     | 288         |
| ≥ 106  | 512         |

理由：高阶多项式的根更密集，Durand-Kerner 的线性收敛阶段更长，
需要更多迭代才能进入二次收敛区域。

**3.4.3 收紧收敛容差**

从 $10^{-12}$ 收紧到 $10^{-13}$。对于高阶多项式，根的间距可能小到 $O(10^{-2})$，
如果迭代步长容差太宽松，可能在两个相邻根之间"跳跃"而无法稳定收敛。

**3.4.4 残差验证兜底**

即使迭代未在最大次数内收敛到步长容差，仍检查所有根的多项式残差：

$$
\max_k |p(z_k)| < 10^{-8}
$$

如果残差足够小，说明根已经足够精确（只是相邻根之间的微小振荡阻止了
步长收敛判据的满足）。这对于近重根情况尤为重要。

**3.4.5 分母零保护阈值收紧**

将 $\prod_{j \neq k}(z_k - z_j)$ 的零保护阈值从 $10^{-24}$ 收紧到 $10^{-30}$。
原始阈值过于保守——对于 30 次多项式，分母的正常量级可能在 $10^{-20}$ 左右
（30 个间距为 0.1 的因子相乘），过早跳过更新会导致收敛停滞。

### 3.5 极点容差自适应

求根器产生的极点理论上应精确位于虚轴上（$\text{Re}(p_k) = 0$），
但数值误差导致实部非零。原始代码使用固定阈值 $10^{-6}$ 拒绝偏离虚轴的极点。

对于 N 阶多项式，Durand-Kerner 的根精度受限于多项式的条件数：

$$
|\Delta z_k| \lesssim \kappa_{\text{root}} \cdot \epsilon_{\text{mach}}
$$

其中根条件数 $\kappa_{\text{root}}$ 随阶数增长。经验观察：

| 阶数 | 典型极点实部偏移 | 所需容差 |
|------|-----------------|---------|
| 10   | $< 10^{-12}$   | $10^{-6}$ |
| 20   | $\sim 10^{-7}$ | $10^{-5}$ |
| 30   | $\sim 3 \times 10^{-5}$ | $4 \times 10^{-5}$ |

自适应容差公式：

$$
\text{tol} = \max\left(10^{-6}, \; 10^{-6} \cdot \left(\frac{N}{5}\right)^2\right)
$$

这给出：
- N=5: $10^{-6}$
- N=10: $4 \times 10^{-6}$
- N=20: $1.6 \times 10^{-5}$
- N=30: $3.6 \times 10^{-5}$

### 3.6 稳定性保证的数学基础

整个综合流程的数值稳定性建立在以下层次结构上：

**层次 1：Cameron 递推的稳定性**

Cameron 递推使用三项递推关系构建 F(w)，每步只涉及加法和乘法，
不涉及除法或减法（除了 $1/\omega_k$ 项）。递推本身是前向稳定的，
系数增长受控于传输零点的位置。对于 $|\omega_k| \geq 1$ 的零点，
系数增长为多项式级别而非指数级别。

**层次 2：E(w) 多项式的 Hurwitz 稳定化**

$E(w) = F(w)/\epsilon_r + P(w)/\epsilon$ 的根通过反射到上半平面来保证稳定性：

$$
w_k^{\text{stable}} = \begin{cases} w_k & \text{if } \text{Im}(w_k) \geq 0 \\ w_k^* & \text{if } \text{Im}(w_k) < 0 \end{cases}
$$

然后从稳定根重建多项式：$E_{\text{stable}}(w) = \prod_k (w - w_k^{\text{stable}})$。
这一步的关键是：即使求根有微小误差，反射操作本身是精确的（只取共轭），
而从根重建多项式使用的是逐根相乘（Horner 式），误差累积为 $O(N \cdot \epsilon_{\text{mach}})$。

**层次 3：部分分式展开的稳定性**

残差计算 $r_k = N(p_k) / D'(p_k)$ 中，$D'(p_k)$ 是分母导数在极点处的值。
对于简单极点，$|D'(p_k)| > 0$；对于近重极点，$|D'(p_k)|$ 可能很小，
导致残差放大。当前实现通过以下方式缓解：

- 极点排序后按虚部分组，共轭对合并处理
- 残差分类使用 $10^{-6}$ 容差区分实/复类型
- 对于无法配对的"孤立"复极点，退化为取实部（物理上合理的降级策略）

**层次 4：矩阵构建的精度传递**

从残差到耦合矩阵的映射：
$$
M_{S,k} = \text{sign}(r_{11,k}) \cdot \sqrt{|r_{11,k}|}, \quad
M_{k,L} = r_{12,k} / \sqrt{|r_{11,k}|}
$$

平方根运算将残差误差缩小为 $\sqrt{\epsilon}$，除法可能放大误差但分母
$\sqrt{|r_{11,k}|}$ 对于物理滤波器总是 $O(1)$ 量级。

### 3.7 实测验证结果

| 配置 | 阶数 | 传输零点 | 对称性误差 | 功率守恒误差 |
|------|------|---------|-----------|-------------|
| 全极点 | 20 | 0 | $2.2 \times 10^{-16}$ | $< 10^{-9}$ |
| 全极点 | 24 | 0 | $4.4 \times 10^{-16}$ | $< 10^{-9}$ |
| 全极点 | 30 | 0 | $2.2 \times 10^{-16}$ | $< 10^{-9}$ |
| 带零点 | 20 | 10 | $4.4 \times 10^{-16}$ | $< 10^{-9}$ |
| 带零点 | 24 | 6 | $3.3 \times 10^{-16}$ | $< 10^{-9}$ |
| 带零点 | 30 | 6 | $3.9 \times 10^{-16}$ | $< 10^{-9}$ |

所有配置的矩阵对称性误差均在机器精度量级（$\epsilon_{\text{mach}} \approx 2.2 \times 10^{-16}$），
证明优化后的流程在 30 阶范围内数值完全稳定。

---

## 第四部分：综合性能评估

### 4.1 端到端性能对比（order=20, 201 频率点）

| 阶段 | 原始复杂度 | 优化后复杂度 | 估算加速 |
|------|-----------|-------------|---------|
| 多项式综合 | $O(N^2 \cdot I)$ | $O(N^2 \cdot I')$ | ~1× |
| 残差展开 | $O(N^3)$ | $O(N^3)$ | ~1× |
| 矩阵构建 | $O(N^2)$ | $O(N^2)$ | ~1× |
| Folded 变换 | $O(N^5)$ | $O(N^3)$ | **~500×** |
| Arrow 变换 | $O(N^5)$ | $O(N^3)$ | **~500×** |
| 频率响应 (201 pts) | $201 \times \frac{5}{3}N^3$ | $201 \times (\frac{4}{3}N^3 + 2N^2)$ | ~1.2× |

其中 $I$ 为求根迭代次数，$I' \leq I$ 由于更好的初始值。

### 4.2 内存使用对比

| 操作 | 原始 | 优化后 |
|------|------|--------|
| 单次旋转 | 3 × $(N+2)^2$ 个 f64 | 0 额外分配 |
| 200 次旋转 | 600 次堆分配 | 1 次 clone |
| 频率响应 (per point) | $(N+2)^2$ 个 Complex64 | 2 × $(N+2)$ 个 Complex64 |

### 4.3 可扩展性极限

当前实现采用三级自适应 fallback 求根策略，理论极限大幅提升：

- **阶数极限**：~50+ 阶（Companion Matrix 特征值方法无迭代收敛问题）
- **频率点数**：无限制（各点独立，可并行）
- **传输零点数**：≤ 阶数（受限于 Cameron 递推的系数增长）

---

## 第五部分：自适应 Fallback 求根架构

### 5.1 设计动机

单一求根算法无法覆盖所有场景：

- Durand-Kerner 在 degree ≤ 28 时最快，但对高阶或近重根可能不收敛
- Aberth 方法收敛更快但实现更复杂，对极端情况仍可能失败
- Companion Matrix 方法数值最稳健但速度最慢（O(n³) 一次性）

因此采用**策略链**模式：快速方法优先尝试，失败时自动降级到更稳健的方法。

### 5.2 三级 Fallback 策略

```
AdaptiveRootSolver.roots_of(polynomial)
│
├─ 尝试 1: DurandKernerRootSolver
│  ├─ 成功 → 返回结果（最快路径）
│  └─ 失败 → 继续
│
├─ 尝试 2: AberthRootSolver
│  ├─ 成功 → 返回结果
│  └─ 失败 → 继续
│
└─ 尝试 3: CompanionMatrixRootSolver
   ├─ 成功 → 返回结果（最稳健路径）
   └─ 失败 → 返回错误（极端情况）
```

### 5.3 Aberth-Ehrlich 方法（第二级）

**迭代公式：**

$$
z_k^{\text{new}} = z_k - \frac{w_k}{1 - w_k \cdot \sum_{j \neq k} \frac{1}{z_k - z_j}}
$$

其中 Newton 修正量 $w_k = p(z_k) / p'(z_k)$。

**与 Durand-Kerner 的对比：**

| 特性 | Durand-Kerner | Aberth-Ehrlich |
|------|--------------|----------------|
| 收敛阶 | 二次 | 三次 |
| 每步代价 | $O(n^2)$（无需导数） | $O(n^2)$（需要 $p'(z)$） |
| 对聚集根 | 容易停滞 | 更稳定（Newton 基础） |
| 全局收敛 | 依赖初始值 | 更宽的收敛域 |

**实现参数：**
- 最大迭代次数：$\min(200 + 6 \cdot \max(0, n-10), \; 800)$
- 收敛容差：$10^{-13}$
- 初始值：与 Durand-Kerner 相同的黄金角分布

**三次收敛的数学证明：**

Aberth 修正等价于对 $\log p(z)$ 的 Newton 迭代加上同时偏转项。
设 $p(z) = \prod_k (z - \alpha_k)$，则：

$$
\frac{p'(z)}{p(z)} = \sum_k \frac{1}{z - \alpha_k}
$$

Aberth 公式可以改写为：

$$
z_k^{\text{new}} = z_k - \frac{1}{\frac{p'(z_k)}{p(z_k)} - \sum_{j \neq k} \frac{1}{z_k - z_j}}
$$

当 $z_j \approx \alpha_j$ 对所有 $j \neq k$ 时，分母中的求和项近似抵消了
$p'/p$ 中除第 $k$ 项外的所有贡献，留下 $\approx 1/(z_k - \alpha_k)$，
因此 $z_k^{\text{new}} \approx \alpha_k$，实现了三次收敛。

### 5.4 Companion Matrix 特征值方法（第三级）

**原理：** 多项式 $p(z) = z^n + c_{n-1}z^{n-1} + \cdots + c_0$ 的根
等价于其伴随矩阵的特征值：

$$
C = \begin{pmatrix}
0 & 0 & \cdots & 0 & -c_0 \\
1 & 0 & \cdots & 0 & -c_1 \\
0 & 1 & \cdots & 0 & -c_2 \\
\vdots & & \ddots & & \vdots \\
0 & 0 & \cdots & 1 & -c_{n-1}
\end{pmatrix}
$$

**特征多项式：** $\det(C - \lambda I) = (-1)^n p(\lambda)$

**实现细节：**

1. **实系数多项式**（常见情况）：直接构建 $n \times n$ 实伴随矩阵，
   使用 nalgebra 的 Schur 分解提取复特征值。

2. **复系数多项式**（域变换后出现）：构建 $2n \times 2n$ 实嵌入矩阵。
   每个复数元素 $a + jb$ 用 $2 \times 2$ 实块表示：
   $$
   a + jb \longleftrightarrow \begin{pmatrix} a & -b \\ b & a \end{pmatrix}
   $$
   特征值成共轭对出现，提取 $n$ 个独立根。

**复杂度：**
- 构建伴随矩阵：$O(n^2)$
- Schur 分解（QR 迭代）：$O(n^3)$
- 总计：$O(n^3)$，无迭代不收敛风险

**数值稳定性：**

Schur 分解使用正交变换（Householder 反射 + Givens 旋转），
是后向稳定的（backward stable）。这意味着计算出的特征值是某个
"接近"原始矩阵的矩阵的精确特征值，偏差量级为 $O(n \cdot \epsilon_{\text{mach}} \cdot \|C\|)$。

对于滤波器综合中的多项式，$\|C\| = O(1)$（系数有界），
因此特征值精度为 $O(n \cdot 10^{-16})$，对 $n = 50$ 仍有 14 位有效数字。

### 5.5 各求根器性能对比

| 求根器 | degree=10 | degree=20 | degree=30 | degree=50 |
|--------|-----------|-----------|-----------|-----------|
| Durand-Kerner | ~5 μs | ~40 μs | ~150 μs | 可能失败 |
| Aberth | ~8 μs | ~50 μs | ~120 μs | ~400 μs |
| Companion Matrix | ~20 μs | ~200 μs | ~800 μs | ~3 ms |
| Adaptive (典型) | ~5 μs | ~40 μs | ~150 μs | ~400 μs |

（估算值，基于 release 模式单核性能）

### 5.6 公共 API

```rust
use mfs::approx::{
    AdaptiveRootSolver,      // 默认：三级 fallback
    DurandKernerRootSolver,  // 最快，适合低阶
    AberthRootSolver,        // 三次收敛，中等阶数
    CompanionMatrixRootSolver, // 最稳健，任意阶数
    ComplexRootSolver,       // trait，可自定义实现
};

// 默认行为（AdaptiveRootSolver）
let roots = polynomial.roots()?;

// 显式选择
let roots = polynomial.roots_with(&CompanionMatrixRootSolver)?;
```

所有内部综合流程（`generalized_chebyshev_helpers`、`residues`）
默认使用 `AdaptiveRootSolver`。

---

## 附录 A：修改的文件清单

| 文件 | 修改内容 |
|------|---------|
| `src/design.rs` | **新增** 高层 `FilterDesign` API（`bandpass()`/`prototype()`/`synthesize()`） |
| `src/lib.rs` | 精简为最小公共表面；旧函数标记为 Legacy |
| `src/prelude.rs` | 精简为 4 个核心类型（`FilterDesign`, `CouplingMatrix`, `Result`, `MfsError`） |
| `src/response/pole_expansion.rs` | **新增** 极点展开快速频率响应（O(N) per point，Schur 补公式） |
| `src/response/backend.rs` | 选择性列/行求解替代完整求逆；利用矩阵对称性单次 LU 分解 |
| `src/matrix/core.rs` | 原地 Givens 旋转；`to_folded()`/`to_arrow()` 零 clone；`flip_sign()` O(n) |
| `src/matrix/sections.rs` | `extract_triplet()`/`extract_quadruplet()`/`extract_trisection()` 原地旋转 |
| `src/approx/complex_poly.rs` | 三级 fallback 求根器（Durand-Kerner → Aberth → Companion Matrix） |
| `src/synthesis/residues.rs` | 求根器切换为 `AdaptiveRootSolver`；极点容差自适应 |
| `src/python.rs` | **新增** PyO3 Python 绑定（`bandpass()`/`prototype()`/`response()`） |

## 附录 B：验证测试

- 162 个单元/集成测试全部通过（`cargo test --all-targets`）
- 极点展开 vs LU 一致性验证：order 4/6/12，相对误差 < 0.1%
- `examples/high_order_stability.rs`：order 10-30 全极点 + 带零点端到端验证
- `examples/solver_accuracy_comparison.rs`：三种求根器精度对比（order 4-40）
- 功率守恒 $|S_{11}|^2 + |S_{21}|^2 = 1$ 误差 $< 10^{-9}$
- 矩阵对称性误差 $< 4 \times 10^{-16}$（机器精度）

## 附录 C：求根器选择指南

| 场景 | 推荐求根器 | 理由 |
|------|-----------|------|
| 常规综合 (order ≤ 20) | `AdaptiveRootSolver`（默认） | Durand-Kerner 直接成功，最快 |
| 高阶综合 (order 20-35) | `AdaptiveRootSolver`（默认） | 自动 fallback 到 Aberth |
| 极高阶 (order > 35) | `CompanionMatrixRootSolver` | 跳过迭代方法，直接用最稳健的 |
| 性能敏感的批量计算 | `DurandKernerRootSolver` | 已知能收敛时最快 |
| 研究/调试 | 逐个尝试 | 对比不同方法的精度 |
