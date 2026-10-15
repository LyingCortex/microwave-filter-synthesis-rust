# MFS 算法参考 (Implementation-Grounded)

本文档描述 `mfs` **代码实际实现的数学**,而非论文里的抽象描述。每一步都给出
公式、变量含义、关键的数值约定/守卫,以及 `file:line` 锚点,便于:

- 把实现对照论文/Cameron 书籍审查;
- 判断某个变换前需要哪种矩阵形式;
- 设计能验证"拓扑改变但响应不变"的测试。

> 与本文档的关系:
> - `docs/cameron-coupling-matrix-algorithms.md` —— Cameron 2003/ZTE 的**论文背景**与流程叙述;
> - `docs/numerical-optimization-report.md` —— 性能优化(O(N) Givens、选择性求解、求根器)的**设计与实测**;
> - 本文档 —— **当前代码到底算了什么**,以源码为准。三者冲突时以本文档 + 源码为准。

---

## 0. 总览:一条流水线

```
FilterSpec / (bands, RL)
   │  approx: 广义 Chebyshev 近似
   ▼
PolynomialSet { E(s), F(s), P(s), eps, eps_r, GeneralizedChebyshevData }
   │  synthesis: 导纳多项式 → 残差展开 → 横向矩阵
   ▼
CouplingMatrix  (Transversal, (N+2)×(N+2), 行主序实值)
   │  transform: 相似变换 (保持响应)
   ▼
Folded / Arrow / CT / CQ
   │  response: S 参数评估 (LU 或极点展开)
   ▼
SParameterResponse { s11, s21, group_delay } @ FrequencyGrid
```

高级功能(复用上面这些原语):
- **双工器 / 多工器**:`synthesis/machiarella.rs` 的 Macchiarella–Tamiazzo 全谱迭代(公共带通域 + 谱分解);
- **多频带**:`approx/multi_band.rs` 的 route-α 双通带等纹波近似;
- **优化 / 提取**:`optimize.rs` 的 Levenberg–Marquardt;
- **谱分解**:`spectral.rs`(cepstrum + 平衡伴随矩阵)。

---

## 1. 符号与全局约定

- 滤波器阶数 `N`(resonator 数);矩阵维度 `side = N+2`(含源节点 row 0、负载节点 row N+1)。
- 归一化低通原型:通带 `|ω| ≤ 1`,归一化频率轴为实数 `ω`(代码里 `normalized_omega`)。
- **s 域约定**:公开多项式用 `s = j·w`(j = √−1, w 为文献里的实辅助变量)。
  - `w → s`:根 `s_root = j·w_root`,即 `(-im, re)`(`generalized_ops.rs:76-78`)。
  - `s → w` 系数:第 k 个系数旋转 `j^k`(`generalized_ops.rs:45-57`)。
- **极点域约定(review P2-12, `response/pole_expansion.rs:81-84`)**:响应求解里存的是
  **实 ω 域极点 `−B_k`**(= 横向矩阵对角 `M[k,k]` 的相反数),`evaluate_at` 用实 `ω` 组
  `denom = ω − pole`。横向矩阵对角存 `−pole.im`(s 域极点 `p = jω₀` ⇒ 对角 `B_k = −ω₀`,
  `transversal.rs:95`),两边一致。
- 矩阵存储:行主序 `data[row*side + col]`。
- **相似度变换**(所有拓扑变换的核心机制):`M' = R·M·Rᵀ`,`R` 正交,只改两个 pivot 行/列。
  它保持特征值 ⇒ 保持理想 S11/S21 ⇒ "改拓扑不改响应"。
- **O(N) 原地 Givens**(`core.rs:328-358`,详见 optimization-report §1):旋转只更新
  2 行 + 2 列,单次 `4·side` 次乘加,零分配。

### 数值阈值(全部出处,便于统一调参)

| 阈值 | 值 | 用途 | 出处 |
|---|---|---|---|
| faer 伴随矩阵后向残差门限 | `1e-3` | 主求根器有效性的硬门限,超了 fall-through | `complex_poly.rs`(FaerRootSolver) |
| 重根聚类(单链) | `1e-5` | 把分裂开的重合根聚回,恢复重数 | `expansion.rs:221` |
| Laurent 主部前导零裁剪 | `max(1e-10·‖·‖, 1e-14)` | 处理分子分母部分相消 | `expansion.rs:368-377` |
| 残差实/复判定 | `1e-6` | classify 的 tolerance | `transversal.rs:46` |
| 极点"实轴噪声"门限 | `1e-10` | `|Im(p)|≤1e-10` 的极点一律当 Real | `expansion.rs:98` |
| 共轭极点对匹配 | `1e-8` | 匹配 `Im(p_a)+Im(p_b)≈0` 且 `Re` 相等 | `expansion.rs:139` |
| 极点实部自适应容差 | `max(1e-6, 1e-6·(N/5)²)` | 高阶 DK 精度退化 O(N²ε) | `transversal.rs:53` |
| 结构零(折叠/箭头) | `1e-7` | to_folded 判定"该不该旋转" | `core.rs:386,417` |
| 结构零(模式校验) | `1e-6` | transform facade 的 pattern 容差 | `verify/mod.rs:51` |
| 小值清零 | `1e-10` | `clean_small_values()` | `core.rs:279-284` |
| 谱分解 para-product 校验 | 相对误差 `< 0.01`(8 点) | 验证 `|D(jω)|² = R(jω)` | `spectral.rs:367-385` |
| 非物理 p0² 拒绝 | `< −1e-12` | machiarella 的 p0² 解 | `machiarella.rs:154` |
| p0² 下溢地板 | `1e-15` | sqrt 良定义 | `machiarella.rs:144` |

## 2. Stage 1 — 广义 Chebyshev 近似 (`approx/`)

入口 `generalized_chebyshev_polynomials(spec)` (`generalized_chebyshev.rs:9`)。
**入口守卫**:先 `validated_transmission_zeros`,再逐个拒绝 `|z| < 1`(通带内)和
`|z| = 1`(带边)的零点 —— typed `InvalidTransmissionZero`。原因:特征函数
`K(s) = F(s)/P(s)` 必须在虚轴上取实值才能满足 **Feldtkeller 恒等式**
`|E(jω)|² = |F(jω)|² + |P(jω)|²/eps²`;`|z|<1` 使 Cameron 递推的
`√(1−1/z²)` 变虚 ⇒ `F(s)` 非 Hermitian ⇒ 恒等式 O(1) 破坏;`|z|=1` 使其为 0
⇒ `F(j)=0` ⇒ eps 无定义(`generalized_chebyshev.rs:16-24`)。

内部 `synthesize_generalized_chebyshev_data` (`generalized_chebyshev_helpers.rs:467`) 依次:

### 2.1 零点补齐
`pad_transmission_zeros`:有限零点 + 用 `∞` 补齐到 N 个(`generalized_chebyshev_helpers.rs:104-134`)。
`reciprocal_or_zero(∞)=0`,`safe_sqrt_term(∞)=1`。

### 2.2 Cameron 递推(在 w 域)
`cameron_recursive` (`generalized_chebyshev_helpers.rs:142`)。初值(第 1 个零点 `z₁`):
```
U₀ = 1 − w/z₁                       V₀ = √(1 − 1/z₁²)
```
对后续每个零点 `zₙ`(含 ∞;`reciprocal = 1/zₙ`(∞→0),`temp = √(1−1/zₙ²)`)迭代:
```
Uₙ = (1 − w/zₙ)·Uₙ₋₁ + temp·(1 − w²)·Vₙ₋₁
Vₙ = (1 − w/zₙ)·Vₙ₋₁ + temp·Uₙ₋₁
```
(代码 `generalized_chebyshev_helpers.rs:179-199` 是就地展开的卷积:U 更新 = `U·(1−w/z)`(长度 +1)
加 `temp·V·(1−w²)`(同下标加、下标+2 减);V 更新 = `V·(1−w/z)` 加 `temp·U`
—— 注意 V 的第二项是 `temp·U`,**不带 w 因子**。)结果 `U` 转成 `F(s)`
(降幂 w → 升幂 s,`descending_w_to_ascending_s`,`generalized_ops.rs:33`)。
**F(s) 的根即反射零点(在 ±j·z 处)**。

### 2.3 P(s)(传输零点)
`find_p_polynomial` (`generalized_chebyshev_helpers.rs:214`):`P(s) = ∏(s² + z_k²)`(有限零点),
若 `(N − finite_count)` 为偶再乘 `j`。P 承载**传输零点**。

### 2.4 eps 与 eps_r(纹波参数)
`find_eps` (`generalized_chebyshev_helpers.rs:247`):
```
rl_linear = 10^(RL/10) − 1
eps   = |P(j)| / |F(j)| / √rl_linear
eps_r = 1                                   （有无穷零点,即 finite_count < N）
      = eps / √(eps² − 1)                   （全有限零点）
```
> **全有限零点的有效域**:`eps_r = eps/√(eps²−1)` 仅当 `eps > 1` 才有有限实值。
近带边零点(`|z| ≈ 1`)或高 RL 会使 `eps ≤ 1`,该配置对全有限公式不可行 ——
`find_eps` 在入口以带域信息的 typed error(`PreconditionViolation`,含 eps 值与
替代建议:改用 mixed 零点,即至少保留一个无穷零点,eps_r ≡ 1)拒绝,而非让 NaN
流到 E 多项式阶段后被通用有限性检查拦下(2026-10-12 修复,review P1-1)。

### 2.5 辅助 A 多项式
`build_a_polynomial_stage` (`generalized_chebyshev_helpers.rs:328`):`P(s)→P(w)`,对每个有限零点求
`Ψ` 求和(`Ψₙ = rₙ·∏_{m≠n}(w²−z_m²)`,`rₙ = zₙ·safe_sqrt_term(zₙ)`),
`A(w) = (N − finite_count)·P(w) + Ψ_sum`,求根后旋回 s 域。全极点时 A = None。

### 2.6 E 多项式(关键:在 w 域做,强制上半平面,再做 monic)
`build_e_polynomial_stage` (`generalized_chebyshev_helpers.rs:407`):
1. `F(s),P(s) → F(w),P(w)`,各自 **normalize to monic**(s→w 旋转会让首项变
   `j^N` 而非 1,不归一会破坏 F/eps_r 与 P/eps 的相对权重 —— 这正是奇数阶曾产生
   虚假通带零点的根因,`generalized_chebyshev_helpers.rs:417-424`);
2. **Cameron 公式 `E(w) = F(w)/eps_r + j·P(w)/eps`**(`generalized_chebyshev_helpers.rs:426-433`)。
   那个 `j` 至关重要 —— 缺了它奇数阶的 P 贡献变实部,E 的根结构全错;
3. 求根(AdaptiveRootSolver),`reflect_to_upper_half_plane` 把根映到稳定上半平面,
   重建、旋回 s 域、**再 monic**(`generalized_chebyshev_helpers.rs:443-447`,w→s 首项带 `(−j)^N`)。
E 的根 = 网络的极点(左半 s 平面)。

### 2.7 求根器链(所有求根共用)
`AdaptiveRootSolver` (`complex_poly.rs:805-823`):**faer 复伴随矩阵特征值优先**
(首一归一 + Parlett–Reinsch 相似度平衡 + 后向残差门限 ≤1e-3 + 非有限检查),
失败 fall-through → **Durand–Kerner** → **Aberth–Ehrlich**。
> 注意:不是"DK 优先"。nalgebra 的 `CompanionMatrixRootSolver` 被**故意弃用**
> (其复系数 2n×2n 实嵌入分支会返回共轭翻转的伪根,有专门回归测试钉死)。

---

## 3. Stage 2 — 导纳多项式 (`synthesis/residues/admittance.rs`)

入口 `synthesize_admittance_polynomials` (`admittance.rs:22`)。把 E、F/eps_r、
P 拼成二端口短路导纳矩阵 `[Y] = (Y11, Y12, Y22)`,共享一个分母。

记 paraconjugate `X~ = X(−s)* = alternating_conjugate(X) · (−1)^N`
(全局 `(−1)^N` 因子是 Cameron 公式与朴素 `alternating_conjugate` 的差别,
`admittance.rs:58-65`;奇数阶缺它会丢阶、极点离开虚轴、得到全反射矩阵)。

**首项抵消守卫**(`admittance.rs:40-54`):防御性检查 —— 仅当 E(s) 首项系数为**负**
且与 `F/eps_r` 首项相消(`|e_lead + f_lead|² < 1e-20` 且 `e_lead.re < −1e-15`)时**取反 E**。
库内部路径上 E 在 `build_e_polynomial_stage` 末尾已 monic(首项恒 `+1+0j`,实测
order 4/6/8),`F/eps_r` 首项 `1/eps_r ∈ (0,1]`,两首项相加恒在 `(1,2]`,守卫
**不触发**(no-op);它只保护经公共 `pub` 字段直接构造非 monic `e_s` 的误用路径
(用户手搓 `GeneralizedChebyshevData` 时 E 首项可为负,守卫此时生效)。
注:旧版注释曾声称 "N mod 4 == 2 时 E_lead = −1 导致相消" —— 那是 re-monic
(`:447`)引入**之前**的状态,现已失效,勿据此"修复"守卫。

四元构造:
```
denominator = E_n + F/eps_r + (F/eps_r)~ + E~
y11 = E_n − F/eps_r + (F/eps_r)~ − E~
y22 = E_n + F/eps_r − (F/eps_r)~ − E~
y12 = −2P/eps                      （承载传输零点;M_SL 从这里或闭式兜底来）
```
归一化 `R_S = R_L = 1`。**降级兜底**(`admittance.rs:84-96`):若四元构造因数值
噪声丢阶(`denominator.degree() < order`),退回简化的两元 Cameron 式
`denominator = E_n + F/eps_r` 等 —— 有意的紧急路径,不是主算法。

---

## 4. Stage 3 — 残差展开 (`synthesis/residues/expansion.rs`)

入口 `synthesize_residue_expansions` (`expansion.rs:170`):对 y11/y12/y22 各自做
部分分式展开(共享同一分母根集)。

1. **求根**:分母全部根(AdaptiveRootSolver)。
2. **重根聚类**(E2 Laurent):求解器把真 m 重根拆成 ~1e-5 内的 m 个根;
   **单链聚类**(根离任一成员 ≤1e-5 即并入,`:234-241`)。极点=簇内均值,
   **重数=簇大小**(不从 shift 后系数向量重推 —— 病态,`:247-252`)。
3. **Laurent 主部**(`principal_part`,`:318-360`):把分子分母平移到极点 p
   (`P(x)=num(p+x)`,`Q(x)=den(p+x)`),对 `R(x)=q_m+q_{m+1}x+…` 做幂级数除法
   `a_k = (p_k − Σ r_i·a_{k−i})/q_m`,得 `Σ b_k/(s−p)^(m−k)`。
   只有 m≥2 才置 `laurent = Some`。
4. **假分数处理**(`:189-213`):`deg(num)>deg(den)` 做长除,**只保留商的常数项**
   作 `constant_term`(它成为横向矩阵的 M_SL);deg 相等则常数 = 首项之比。

### 残差分类 `classify_residues` (`:80`, tolerance = 1e-6)
- 极点 `|Im(p)| ≤ 1e-10` ⇒ **恒 Real**(那里的虚部是求根噪声);
- 否则残差 `|Im| > tol` 判为复,按 `|Im(p_a)+Im(p_b)|<1e-8 且 |Re 差|<1e-8`
  配对共轭;
- **未配对的复极点降级为 Real**(不报错)—— 对 alternating_conjugate 构造出的
  复系数分母是有意的兜底(`:150-160`);
- `residue_real_component` (`:63`):优先取 `.re`,仅当 `|re|≤1e-10` 才取 `.im`
  —— 吸收奇数阶全极点残差纯虚的 `j^N` 相位(直接取 `.re` 会把所有耦合清零)。

> **每个 Y 独立求 Laurent 长度**(review):同一共享重极点在 y11/y12/y22 经过
> 分子分母相消后有效重数可以不同,下游 `pole_expansion.rs:183-199` 因此对每个 Y
> 用**自己的** Laurent 长度求和,而非统一重数。

---

## 5. Stage 4 — 横向矩阵组装 (`synthesis/residues/transversal.rs`)

入口 `build_transversal_from_residues` (`transversal.rs:12`)。横向矩阵 = N 条一阶
支路并联在源/负载之间:**无谐振器间耦合**,只有 `M_Sk`、`M_kL`、对角 `M_kk`、
可选 `M_SL`。

**守卫**:
- 三个展开各必须恰有 N 个残差(`:16-24`);
- **任一残差 `laurent.is_some()` ⇒ `NotImplemented`**(`:33-43`):横向链每个谐振器
  只表示一个简单极点,`1/(s−p)^m (m≥2)` 不可表示 —— 显式报错而非静默构造错矩阵。
- 极点实部超自适应容差 `max(1e-6,1e-6·(N/5)²)` ⇒ `NumericalFailure`(`:53,65-72`)。

**每个极点(Real 或 ComplexPair 中逐个)填一行**(`:57-168`):
```
按 |r11| ≥ |r22| 选支路:
  若用 r11:  M_Sk = √|r11| ,        M_kL = r12 / √|r11|
  否则:       M_Sk = r12 / √|r22| ,  M_kL = √|r22|
对角:         M[k,k] = −pole.im      （= 实 ω 域失谐 B_k）
对称写入:     set(k,k,−pole.im); set_symmetric(0,k,M_Sk); set_symmetric(k,N+1,M_kL)
```
源耦合恒为所用自残差的**正平方根**;负载耦合携带 y12 的符号。
`nonzero_sqrt_abs` 拒绝 `|残差| ≤ 1e-12`。

**直接源-负载耦合 M_SL**(`:169-184`):
- 若 `y12.constant_term` 存在:其实部必须 `≤ 1e-6`(否则报错 —— 该路径下直接项
  必须纯虚),`M_SL = constant.im`;
- 否则,当**所有零点有限**(`len == order`,即无无穷零点)时用闭式
  `M_SL = eps·(eps_r − 1)/eps_r`。

> **偏差提示**:`combine_conjugate_pair` (`:190-246`, `#[allow(dead_code)]`) 描述的是
> "2·Re(r)" 合并捷径,但**现行代码不加倍** —— ComplexPair 的两个极点各按上面的
> sqrt/ratio 逻辑独立处理,等效贡献是两条独立行之和。不要把这个 dead helper 当
> 成已实现的路径来描述。

**不变量(测试钉死)**:order-4、±1.5 两零点、20 dB → 2 个 ComplexPair;所有源
耦合为正;`M_Sk·M_kL` 符号匹配 `2·Re(r12)`;全矩阵对称到 1e-10。

---

## 6. Stage 5 — 拓扑变换 (`transform/` + `matrix/core.rs`)

`transform_matrix(m, TopologyKind)` facade 只做两件事:调 `m.transform_topology`
再附加**模式校验**(默认不附响应校验,`transform/mod.rs:64-115`)。
`transform_topology` (`core.rs:218-230`) 是个纯 match:
`Folded → to_folded()`,`Arrow → to_arrow()`,`Wheel → to_arrow()+改标签`
(即 **Wheel == Arrow**,core.rs:224-227),`Transversal → clone`。

### 6.1 相似度符号翻转 `flip_sign_inplace`(`matrix/mod.rs:42`,共享原语)
对谐振器 `index` 做 `R·M·Rᵀ`(`R[index][index]=−1`):**整行取负 + 整列取负**。
对角元被取负两次 ⇒ **净不变**。绝不能取负第三次(历史 bug:偶数阶非对称滤波器
丢传输零点,见 `docs/folded-fidelity-flip-sign-fix.md`)。
- `core.rs::flip_sign`(`&self→Self`)= clone + 调共享原语(不可变包装,`to_folded` 用);
- `sections.rs::ensure_positive_main_couplings`(CT/CQ 收尾用)直接调共享原语。
回归测试 `flip_sign_negates_row_and_column_only_and_preserves_diagonal` 钉死
对角不变 / 仅行列变号 / 两种形式逐位一致 / 双翻为恒等。

### 6.2 `to_folded`(`core.rs:372-437`)— 三段 Givens + 符号规范
`end_oper = N/2`(偶)或 `(N−1)/2`(奇)。
1. **列消去(源耦合)**:`row_oper_num = 0..end_oper`,`col` 从 `N−row_oper_num`
   降序到 `row_oper_num+2`;若 `|M[row,col]| > 1e-7`,在平面 `(col−1, col)` 做
   列操作,`θ = safe_angle(−M[row,col], M[row,col−1])`,消去 `M[row,col]`。
   净效果:保留 `M[0,1]`,消掉 `M[0,k], k≥2`。
2. **行消去(负载耦合)**:镜像列 `col_oper_num = side−1−row_oper_num`,
   `row_in_col_oper` 升序;平面 `(row, row+1)`,`θ = safe_angle(M[r,c],M[r+1,c])`。
3. **正主链符号规范**:`index = 0..side−1`,若 `M[index,index+1] < 0` ⇒
   `flip_sign(index+1)`。
4. `clean_small_values()`(≤1e-10 清零)+ `topology = Folded`。
`safe_angle(y,x)`(`rotations.rs:30-38`):`|x|<1e-10` 时按 `y` 符号返回 ±π/2,
否则 `atan(y/x)` —— 避免除以极小邻耦合。

### 6.3 `to_arrow`(`core.rs:439-467`)
对 `resonator = 1..N`、`target = resonator+1..=N`,平面 `(resonator, target)`,
`θ = safe_angle(−M[res−1, target], M[res−1, resonator])`,把每个谐振器到更高号
谐振器的耦合**聚到最后一个谐振器列**(箭头形)。**无条件旋转**(已 ≈0 时 θ≈0
为 no-op),**无正主链符号阶段**。收尾 `clean_small_values()` + `topology = Arrow`。

### 6.4 级联三节 CT / 四节 CQ(`matrix/cascade.rs` + `sections.rs`)
两条路径都**先转到 Arrow**,再抽取。
- **CT** `to_cascaded_trisections(zeros)`:守卫 `zeros` 非空、`order ≥ zeros.len()+2`
  (即 max 零点 = N−2,写法避免 order=1 的 usize 下溢)。然后**单次**
  `arrow.extract_triplet(zeros[0], 1)` —— 条件化负载尾 + 逐左移 peel 到源端。
  `zeros[0]` 是主零点(源侧、Fig.21 锚定);其余零点已编码在 arrow 里,只用于计数校验。
  > 历史的"逐零点循环"(每趟重新条件化尾)被**删除**(review P2-4):它会切断主链、
  > 产出稠密矩阵。
  `extract_triplet`(`sections.rs:11-56`):`θ = safe_angle(M[t−1,t], z + M[t,t])`
  消掉 (t−1,t) 主链(在零点频率把尾三节与尾端解耦),然后逐行旋转把三节移向源端。
- **CQ** `to_cascaded_quadruplets(zero_pairs)`:守卫 `order ≥ 4·pairs`。**源侧不重叠
  平铺**:第 i 个 quadruplet(0-based)span 谐振器 `4i−3..=4i`,每对推两个 z_loc
  `[s, s+2]`(span 2 = 伴生三节)和 `[s, s+3]`(span 3 = quadruplet 对角),`s=4i+1`。
  `extract_cascade_quadruplets`(`sections.rs:226-344`)按 span 分类,每个 quadruplet
  找伴生三节(**源侧平铺下 `common_resonator` 恒 = 1**,`=4` 负载端公式不可达),
  所有操作**按位置从尾到顶交错**执行(尾先),`extract_quadruplet` = 两个相邻三节 +
  一次合并旋转(消掉两节间的 inline 耦合,留下 2-skip 交叉作 quadruplet 签名),
  收尾 `ensure_positive_main_couplings()`。
  验证布局:order-6 1 对 → `[[1,3],[1,4]]`;order-8 2 对 → `[[1,3],[1,4],[5,7],[5,8]]`。

CT/CQ 都是相似度变换 ⇒ S21 幅度在 21 点网格上保持到 <1e-6(有测试)。

---

## 7. Stage 6 — S 参数响应 (`response/`)

### 7.1 LU 路径(权威,所有拓扑通用)`backend.rs`
每个频率点:
```
A(ω) = M + ω·I_r − jR_S·e₀e₀ᵀ − jR_L·eₙeₙᵀ      （I_r 只在谐振器对角为 1）
有耗(Q 有限):谐振器对角再 − j/Qu  （backend.rs:130-138,把极点移出虚轴）
S11 = 1 + 2jR_S·[A⁻¹]₀₀
S21 = −2j√(R_S R_L)·[A⁻¹]_{N+1,0}
群延迟 τ = Im( Σ_k [A⁻¹]_{N+1,k}·[A⁻¹]_{k,0} / [A⁻¹]_{N+1,0} )   （= −Im(d/dω ln S21)）
```
**选择性求解**(`backend.rs:172-194`):只解 `A x = e₀`(第一列)和 `A y = eₙ`(最后一列;
因 A 对称,最后一行 = 最后一列),**一次 LU 分解解两个右端**,省掉完整求逆的 N−2 次回代。
> 偏差:LU 群延迟利用对称取 `最后一行=最后一列`(`backend.rs:164-169`),但**消费端不
> 验证对称** —— 若经 `new_with_topology_asymmetric` / `allow_asymmetric(true)` 喂入非对称
> 矩阵,S11/S21 幅度仍对、群延迟会**静默算错**(review 0.6.1 P1-3)。构造端现在有对称
> 守卫兜底,但这是已知不变量。
`ResponseSolver::evaluate_normalized(_with_settings)` (`response/mod.rs:169-185`)。

### 7.2 极点展开快路径(O(N)/点,仅无耗横向)`pole_expansion.rs`
`PoleExpansionData::from_matrix` (`pole_expansion.rs:100`) 从横向矩阵**直接读极点/残差**(不求根):

**结构性前置条件**(`assert_transversal_structure`,`:164`,round-3 修复):
该路径只能表达**横向**矩阵(谐振器块对角 + 每谐振器只耦合 Source/Load)。
折叠/箭头等拓扑携带的谐振器间耦合 `M[k,k']` 无法表达 —— 若直接读取会被**静默丢弃**、
给出物理上错误(但静默算出)的响应。`from_matrix` 现在先扫描所有谐振器块非对角项与
**端口对角** `M[0,0]`/`M[N+1,N+1]`,任何超过 `1e-9·max(1, max|M|)` 者(远高于 f64
舍入 ~1e-16,远低于任何 O(1) 耦合)即以 typed `PreconditionViolation` 拒绝 ——
查的是**元素**而非拓扑标签,故标着 `Transversal` 却仍含谐振器间耦合的矩阵也拦得住。
拒绝后主流程落回 LU(任意对称矩阵都正确,只是慢)。5 个回归测试钉住该契约。

读取本身(不求根):
```
极点  pole_k  = −B_k            （实 ω 域,= 对角相反数,P2-12 约定）
残差  y11: M_Sk² ,  y12: M_Sk·M_kL ,  y22: M_kL²
常项  y12 += M_SL
```
`evaluate_at(ω)` (`:233`) 在 Schur 2×2 型下解析求:
```
S11 = 1 + 2j·[Σ r11,k/(ω−pole_k) + const11] / D(ω)
S21 = −2j·[Σ r12,k/(ω−pole_k) + M_SL] / D(ω)
D(ω) = 1 + Σ r22,k/(ω−pole_k) + …（网络分母）
群延迟 = −Im(d/dω ln S21) = 对同一单极点项求导累加
        d/dω [r/(ω−p)] = −r/(ω−p)²        （简单极点）
        d/dω [b/(ω−p)^(m−k)] = −(m−k)·b/(ω−p)^(m−k+1)   （Laurent,每 Y 用自己长度）
```
单遍 O(N)。**在传输零点(S21→0)群延迟返回 `None`**(`inv_10.norm_sqr() ≤ 1e-18`,
`:330`),与 LU 的退化守卫一致 —— 显式未定义,不是假 0。

### 7.3 快路径的信任策略(`design.rs` `fast_path_gate`,`:72`)
`response_normalized` / `response_bandpass` 先建极点展开,在 **5 点检查网格**上与 LU
比对 `fast_path_gate`:比较**复数** S21(幅度**与**相位)与 S11,任一点
`max(|ΔS21|/|S21_LU|, |ΔS11|/|S11_LU|) < 1e-3`,且**至少比对了 1 个**
`|S21| > 1e-10`(深阻带噪声地板)的点、所有值有限,才信任;否则整格落回 LU。
纯幅度门会放过"幅度不变但相位旋转/仅 S11 畸变"的错误,复数门(含 S11)堵住了
(round-3 修复,5 个测试钉住,含"等幅旋转相位"回归用例)。信任后单点若近奇异仍逐点落回 LU。

### 7.4 群延迟单位约定(`response/mod.rs:15-35`)
`group_delay` 是**归一化原型单位** `dτ/dω`(ω = 归一化坐标,`normalized_omega`),
两条求解路径同单位。**换算成秒**:对带通映射 `ω = (1/FB)(f/fc − fc/f)`,
`dω/df = (1/FB)(1/fc + fc/f²)`,故 `秒 = group_delay·(1/FB)(1/fc+fc/f²)/(2π)`。
> **陷阱(2026-10-11 记录)**:数值 `−d(phase)/df` 是"每赫兹",比"每秒"大 2π 倍。
> 直接和秒值比会读出 `(2π−1)·gd ≈ 85 ns` 的假偏差。判据:该偏移**随 |S21| 增大**、
> **细化网格不缩小** —— 指向单位/约定错配而非算法 bug。对比时先乘 `dω/df` 或除以 2π。

---

## 8. 高级功能

### 8.1 Macchiarella–Tamiazzo 全谱迭代(`synthesis/machiarella.rs`)
双工器(N=2,与历史路径**逐位一致**)与多工器(N≥2)共享的固定点迭代,
把 N 个通道耦合到一个公共 junction。

**公共带通域**:单一变量 `s = j·Ω(f)`,`Ω(f) = (f/f₀ − f₀/f)/b`,
`f₀ = √(f_L·f_H)`、`b = (f_H−f_L)/f₀` 取**最外侧两通道**的边,使两个极端边落在 `s=±j`。
每通道原型 `P_p(ω)` 经仿射 `P_c(Ω) = P_p(AᵢΩ + Bᵢ)` 合成到公共域,
`Aᵢ = (Ωᵢ,hi − Ωᵢ,lo)/2`,`Bᵢ = (Ωᵢ,hi + Ωᵢ,lo)/2`。

**固定分子 N(Ω)**(`:277-301`,循环外建一次):所有通道**反射零点**(各 F 的根,
在原型虚轴)映到公共域 `j·(bᵢ + aᵢ·r.im)`;Type-II junction 额外加实零点 `1/c0`。

**初始 S 状态**:`Sᵢ = (Eᵢ + pad(Fᵢ, deg Eᵢ))/2`(单带综合的左半平面特征多项式)。

**junction 支持**:理想 Type-I(n=1, b0=0)与 Type-II(c0>0)对任意 N≥2;
**非理想 Type-I 且 N>2 被拒绝**(`PreconditionViolation`,无论文依据不臆造数学)。

**每迭代**(`run_full_iteration`,`:540-695`):
1. 原型根 → 公共 Sᵢ(Ω)、公共 Pᵢ;
2. **传输多项式** `Ptᵢ = Pᵢ·∏_{j≠i} Sⱼ`;
3. **p0² 解**(`:347-435`):在每通道**最外侧边** eᵢ 上列 RL 条件
   `Σⱼ |Ptⱼ(eᵢ)|²·p0ⱼ² = |N(eᵢ)|²·(1−εᵢ)/εᵢ`,`εᵢ = 10^(−rlᵢ/10)`。
   对角主导 ⇒ 对角占优。**N=2 用 Cramer 规则逐位一致**(det |·|<1e-30 ⇒ 均分);
   **N≥3 用 nalgebra LU 部分主元**。`p0² < −1e-12` ⇒ typed 拒绝(通道强耦合/重叠);
   `(−1e-12, 0]` 当舍入、地板到 1e-15 保证 sqrt 良定义。
4. **R**(`:576-585`)= `N·Ñ + Σᵢ p0ᵢ²·Ptᵢ·Pt̃ᵢ`(~= alternating_conjugate)。
   在虚轴 `R(jω) ≥ 0` ⇒ **para-identity `D·D~ = R`**,D 为 LHP 谱因子。
5. **D** = `spectral_factor(R, deg N)`;
6. **∏Sᵢ**(`:438-474`):Type-I `(a·N + b·D)/2`(理想 a=b=1 ⇒ `(N+D)/2`,junction
   关系 `N + D = 2·∏S`);Type-II `(c0/2)(D_monic − N)`(D 的领先项与 N 抵消,
   留下 `(2/c0)·∏S`)。**D 先 re-monic**(spectral_factor 只修相位不修幅值),
   `|lead| < 1e-15` ⇒ `NumericalFailure`(避免 200 十年的缩放灾难);
7. **根分配** `assign_blocks`(`:476-505`):∏Sᵢ 的根按 `(imag, real)` 升序,取
   **连续块**(块大小 nᵢ,按通道中心升序);映射回原型(`to_proto` 丢弃非有限
   与 `|r|≥50`);丢根 ⇒ typed `NumericalFailure` 指名越界根;
8. **收敛**(`:637-690`):`max_delta` = 各通道对齐根位移的最大值;
   接受当 `max_delta < tol`(严格)**或**(高阶停滞)迭代 ≥5 且连续 3 次不改善
   且 `best_delta < 1e-3`。
   > **一次性捕获滞后(§C.5)**:返回的 `common.{d,pt,p0_sq}` 描述第 k 次迭代,
   > 而 `s_roots` 是 k+1;滞后被报告的 convergence_error 界定,有意不重算。

**响应**:`junction_response` 扫 `S_ant = N/D`、`S_i = √p0ᵢ·Ptᵢ/D`(在 `s=jΩ(f)`)。

### 8.2 谱分解(`spectral.rs`)
问题:给定 para-Hermitian `R(s)`(`R(jω)≥0`),求度 = target、全左半平面的
`D(s)` 使 `R = D·D~`。`spectral_factor`(`:391`)调度:**先 cepstrum**,
`para_product_verified` 不过则退回伴随;伴随结果**必须**过校验,否则 typed 报错
(不静默出错因子)。两者输出过 `monic()`(相位归一,`|lead|` 保留)。

- **Cepstrum**(`:33-113`,无求根):FFT 大小 `next_pow2(max(8·deg,128))`;
  在 `ω∈[−ω_max,ω_max]` 直采 `|R(jω)|`,`log`(地板 −690=ln 1e-300);
  IFFT → 实 cepstrum;**最小相位窗**(causal 半保留,Nyquist 与 cep[0] 各半);
  前向 FFT → log D(jω) → exp;**系数恢复**用最小二乘(Vandermonde 法方程 LU,
  在良好条件子集上,不是精确 DFT 反卷积)。
- **平衡伴随**(`:205-239`):AdaptiveRootSolver(同 §2.7)求 R 全部根;
  根成 `(z, −z*)` 对,**LHP 选择**(`select_lhp_roots` `:242-345`):按实部升序,
  取 `re<−1e-12` 并标记其镜像已用;不足时把纯虚根 `re≤1e-12` 微 nud 到
  `−1e-14` 入 LHP;最后兜底只补 `re<−1e-9` 的真 LHP 根(防止静默拉进 RHP/纯虚
  根 → 非最小相位)。`s=0` 增益校正 `scale = √(R(0)/(D(0)D̃(0)))`。
- **校验** `para_product_verified`(`:367-385`):虚轴 8 点,
  相对误差 `(|D|²−R|)/R`(R 小时用绝对)< 0.01。
  > **局限(有测试钉死)**:这是**幅值**校验 —— 对实 R,最小相位与最大相位反因子
  > 共享 `|D(jω)|²` 都能过;最小相位由 `select_lhp_roots` 的 LHP 偏置强制,不由该校验。

### 8.3 多频带 route-α 双通带(`approx/multi_band.rs`)
`multi_band_chebyshev_polynomials`(`:338`)。特征函数 `K(w) = F(w)/P(w)`
(实根、对称)在**单一正轴区间**上强制等纹波。反射零点 `F(w)=∏(w²−r_k²)`,
`n = order/2` 个正零点(阶数必须偶)。一个正区间 (a,b) 上 |K| 在 n+1 处取峰
(n 零点之间 + 两端),等纹波 = n+1 方程 / n+1 未知(含共享电平)⇒ **只对 M=1
个正轴区间是平方(可解)**。

**LM 求零点** `solve_rz`(`:115-203`):20000 点网格找 n+1 峰,残差 = 各峰偏离均值;
前向 FD Jacobian `δ_j = 1e-7·max(1,|rz_j|)`;正规方程 `(JᵀJ + λ·diag(JᵀJ))·dp = −Jᵀr`
(高斯消元部分主元);可行域约束(零点须严格在 (lo,hi) 内);接受当 trial cost 下降
(λ×0.5)否则 λ×2,λ>1e12 中断;末门 `√cost > 1e-6` ⇒ `NumericalFailure`。

**K ≥ 3 的 DOF 墙**(`:383-392`):若正轴通带区间数 `pos.len() ≠ 1` ⇒
**`NotImplemented`**,明说"实根特征只能在单一正轴区间等纹波;K≥3(跨原点中频带或
2+ 正轴带)过定 M−1 个,需复根(Zolotarev)或多函数特征"。约束计数:等纹波施加
`N/2 + M − 1` 个独立约束于 `N/2` 个零点,仅 M=1 平方。

**Feldtkeller 装配**:解出正零点后 `F_s = ∏(s² + r_k²)`,P_s 由传输零点
(±jz 对,原点零点单根;去重避免对称列表 {0,±z} 双计),`find_eps` +
`build_e_polynomial_stage` 复用 §2.6,输出 `GeneralizedChebyshevData`
(A=None, route-α 不走 A 阶段),`PolynomialSet` 直接喂 `MatrixSynthesisEngine` 不变。
锚定: Cameron 2011 Fig.25 双通带(阶 10, RL 22, bands [−1,−0.35]∪[0.35,1])
解出零点与 oracle 差 ≤5e-4(`:455-475`)。

### 8.4 优化 / 提取(`optimize.rs`)
`optimize_matrix(initial, target, grid, config)`(`:337`);
`tune_matrix`(目标响应匹配)/`extract_matrix`(从测量 S2P 提取,模板=粗 inline 链)。
- **参数(冻结非零模式)**:`extract_parameters`(`:61-108`)按固定序取
  `|val|>1e-15` 的源/负载/交叉/对角/直接耦合;`rebuild_matrix` 把参数对称写回模板
  的 clone,≤1e-15 的冻结为模板值。⇒ 恢复矩阵**共享模板的稀疏模式**,
  `extract_matrix` 只能拟合横向族。
- **残差**:`4·|grid|` 个(S11 re/im、S21 re/im 全复数,非幅度),
  当前响应经 `evaluate_normalized_with_settings`。目标 `Σr²` 对单分量符号翻转
  不敏感(相位符号歧义不可钉)。
- **Jacobian**:前向 FD `δ_j = 1e-7·max(1,|param_j|)`;扰动求值失败(远离解盆奇异)
  时**步长减半重试 ≤4 次**。
- **LM 步**:解 `(JᵀJ + λ·diag(JᵀJ))·Δp = −Jᵀr`(Marquardt 对角阻尼),高斯消元
  部分主元(枢轴 <1e-20 跳过)。
- **接受/收敛**(`:373-427`):接受当 `trial_cost < cost`(λ×0.5)否则 λ×2,
  λ>1e12 中断;**收敛判据**只在"上一步被接受"且 `|prev_cost − cost| < tol`
  时成立(对运行最优比是错的 —— 接受后次步 cost==best 会误判;对拒绝步判也错 ——
  那是停滞不是收敛)。

---

## 9. 已知限制与不变量(与代码一致)

1. **通带传输零点(`|z|<1`)不支持** —— 单函数实根广义 Chebyshev 原型无法实现
   带内零点(Feldtkeller 破坏,E1 2026-10-05);需多函数特征(新研究)。
2. **横向综合的重复极点** —— 残差展开能处理(Laurent),快路径精确求和,但横向
   链每谐振器只表示一个简单极点 ⇒ typed `NotImplemented`。
3. **Wheel 拓扑** 已弃用,当前是 Arrow 别名。
4. **LU 群延迟依赖对称** —— 非对称输入(经 opt-out 构造)会静默算错群延迟,
   幅度不受影响(review 0.6.1 P1-3)。
5. **多工器简化路径保真度有限**(相邻通道),全谱路径 `full_iteration(true)`
   已收敛并 gate 锚定(Type-I 阶 2–8 + Type-II + N=3/4)。
6. **多频带仅双通带(K=2)**,K≥3 是 DOF 墙(typed 拒绝)。

---

*本文档基于 `7dd617a`(0.6.3-unreleased)源码核对;行号锚点以该版本为准,
重构可能使其漂移,但算法描述不受影响。*
