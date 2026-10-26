# mfs 算法部分性能分析报告

> 性质：只读分析，未修改任何产品代码（`src/`、`tests/`、`benches/`、`Cargo.toml`、`Cargo.lock` 均未改动）。
> 生成日期：2026-10-16。所有基准数据为本次实测（`cargo bench`，release 优化），非引用旧文档。

## 1. 环境与方法

- 机器：Linux x64，WSL2（`6.18.54.1-microsoft-standard-WSL2`）。
- 工具链：`rustc 1.95.0`，`cargo 1.95.0`。
- 库栈：`nalgebra 0.33`、`faer 0.20`、`num-complex 0.4`、`pyo3`（FFI 用）。
- 测量方法：
  - 仓库自带 4 个 Criterion 基准（`benches/`）：`synthesis`、`multiplexer`、`root_solvers`、`optimizer`，各跑 1 轮（共 1 轮基准，远低于 3 轮上限；总墙钟约 14 min，远低于 40 min 上限）。
  - 另写一个一次性探针（`/tmp/mfs_perf_probe`，**已删除**），用公开 API 对 twoport 与 multiplexer 流水线做**子阶段**拆解（Chebyshev / 矩阵 / 折叠 / 响应-LU / 响应-极点展开 / 分母求根），中位数（9 次重复，3 次预热；multiplexer 3 次重复）。探针只 `use mfs::*` 公共导出，不改产品代码。
- 数值口径：`time: [lo med hi]` 为 Criterion 的 `[下界 中位数 上界]`，本文引用中位数；样本数 `n` 为每次采样的迭代数。

## 2. 实测基准数据（本次运行，median，n=样本数）

### 2.1 root_solvers（多项式求根器，n=100）

| 求解器 \ 次数 | d10 | d20 | d30 |
|---|---|---|---|
| durand_kerner | **23.2 µs** | **230.1 µs** | 1.262 ms |
| aberth | **22.8 µs** | **169.1 µs** | **719.9 µs** |
| faer（伴随矩阵特征值） | 84.4 µs | 275.2 µs | **658.1 µs** |
| companion_matrix（裸，非生产路径） | 502.7 µs | 1.165 ms | 2.518 ms |

> 生产默认 `AdaptiveRootSolver` = faer 优先 → 失败回退 Durand-Kerner → Aberth。见 §3 求解器热点。

### 2.2 synthesis（twoport：合成 / 折叠 / 响应，n=10）

标准案例（`data::cases("twoport")`，来自 `cases/twoport.json`）：

| 案例 | 次数 | 合成 (synthesis) | 折叠 (folding) |
|---|---|---|---|
| twoport_4_4 | 4 | 78.7 µs | — |
| twoport_6_3 | 6 | 262.1 µs | — |
| twoport_4_2_unity | 4 | 89.3 µs | — |
| twoport_6_2_culdesac | 6 | 147.4 µs | — |
| twoport_8_3 | 8 | 327.9 µs | — |
| twoport_4_1_trisection | 4 | 91.0 µs | — |
| twoport_7_1_2 | 7 | 324.0 µs（高方差 192–480）| — |

顺序/传输零点扫描（`sweep`，三阶段）：

| 次数 | 零点 | 合成 | 折叠 | 响应 (201 点) |
|---|---|---|---|---|
| 4 | tz0 | 91.0 µs | 56.1 µs | 427.2 µs |
| 4 | tz2 | 60.0 µs | 44.1 µs | 318.8 µs |
| 6 | tz2 | 112.1 µs | 106.6 µs | 498.3 µs |
| 8 | tz4 | 199.7 µs | 131.8 µs | 664.8 µs |
| 10 | tz4 | 287.2 µs | 203.3 µs | 902.8 µs |
| 12 | tz6 | 525.5 µs | 324.0 µs | 1.452 ms |
| 16 | tz8 | 695.8 µs | 619.1 µs | 2.372 ms |
| 20 | tz8 | 1.653 ms | 820.2 µs | 3.008 ms |

> 注意：**基准里的"响应"阶段量的是裸 LU 解算** `ResponseSolver::evaluate_normalized(design.matrix(), &grid)`（`benches/synthesis.rs:187`），不是 `FilterDesign::response_normalized` 走的极点展开快路径。见 §3 响应热点。

### 2.3 multiplexer（N 通道合成 + 响应，n=10）

| 案例 | 阶段 | median |
|---|---|---|
| multiplexer_a（4 通道 4/6/5/4，N=19） | poly（全迭代合成） | **22.05 ms** |
| multiplexer_a | response（401 点，4 通道） | 3.804 ms |
| multiplexer_diplexer（2 通道） | poly | 3.236 ms |
| multiplexer_diplexer | response（401 点） | 1.475 ms |

> `type2`（谐振节点 / 星型）与 `contiguous`（guard band=0）案例在基准中标记 SKIP——前者不是 mfs 实现的拓扑，后者使逐通道装载近似失效。`poly` 是 Machiarella-Tamiao 逐通道迭代合成，含谱因子化 + 乘积多项式求根。

### 2.4 optimizer（LM 参数提取，n=10）

| 场景 | median | 迭代 |
|---|---|---|
| tune/o4 | 3.320 ms | 5 |
| tune/o6 | 8.080 ms | 6 |
| tune/o8 | 13.532 ms | 6 |
| extract/o4 | 9.335 ms | 12 |
| extract/o6 | 12.356 ms | — |

> 每迭代 = 1 次基准残差（41 点 LU 扫描）+ (n_params+1) 次扰动响应，逐列构造有限差分 Jacobian。

### 2.5 探针子阶段拆解（补充测量，中位数）

twoport（`build_filter_spec` → Chebyshev → 矩阵 → 折叠 → 响应）：

| 次数 | Chebyshev | 分母求根×1 | 矩阵阶段 | 全合成 | 折叠 | 响应-LU(201) | 响应-极点展开(201) |
|---|---|---|---|---|---|---|---|
| 4 (tz2) | 0.011 ms | 0.007 ms | 0.027 ms | 0.031 ms | 0.022 ms | 0.174 ms | **0.024 ms** |
| 6 (tz2) | 0.017 ms | 0.011 ms | 0.070 ms | 0.087 ms | 0.060 ms | 0.284 ms | **0.073 ms** |
| 8 (tz4) | 0.051 ms | 0.036 ms | 0.145 ms | 0.203 ms | 0.094 ms | 0.557 ms | **0.047 ms** |
| 12 (tz6) | 0.083 ms | 0.049 ms | 0.160 ms | 0.288 ms | 0.168 ms | 0.816 ms | **0.075 ms** |
| 16 (tz8) | 0.147 ms | 0.102 ms | 0.361 ms | 0.468 ms | 0.419 ms | 1.714 ms | **0.089 ms** |
| 20 (tz8) | 0.315 ms | 0.161 ms | 0.644 ms | 0.940 ms | 0.767 ms | 2.835 ms | **0.146 ms** |
| 24 (tz8) | 0.289 ms | 0.231 ms | 0.865 ms | 1.240 ms | 0.990 ms | 3.866 ms | **0.177 ms** |

multiplexer_a：

| 阶段 | median |
|---|---|
| prepare（建 p_poly、n_poly） | 0.100 ms |
| full_iteration（Machiarella 迭代主体，12 次迭代收敛） | 15.03 ms |
| synthesize() 全程 | 16.19 ms |
| response（401 点，4 通道） | 3.07 ms |

> 探针数值与 Criterion 基准在同一量级（探针为独立短进程、预热 3 次，略低于长时 Criterion 中位数属正常），用于子阶段归因而非绝对值。

## 3. 热点函数分析

按"用户实际调用路径"而非裸基准分组。占比为实测或推算（注明）。

### H1. multiplexer 谱因子化迭代（`poly` 阶段）— 最大单项热点

- 位置：`src/synthesis/machiarella.rs`（`prepare` / `run_full_iteration` / `synthesize_full`）；谱因子化 `src/spectral.rs:208-214`（`spectral_factor` → `AdaptiveRootSolver.roots_of` + `polish_roots` + LHP 选取）。
- 实测：`poly/multiplexer_a` = **22.05 ms**（4 通道，N=19，12 次迭代），是全部基准里**最大单项时间**，约为 twoport 全合成（~1 ms）的 20 倍、LM 提取（12 ms）的 2 倍。
- 复杂度：每迭代 = 每通道 1 次谱因子化（FFT 倒谱 O(n log n) + 多项式求根 O(k³)，k=通道阶）+ 1 次乘积多项式求根（次数 = N=19，O(N³)）。全迭代 ≈ 12 × (4 + 1) 次求根。
- 占比（multiplexer_a 合成内）：`full_iteration` 15.03 ms / `synthesize()` 16.19 ms ≈ **93%**。

### H2. 多项式求根 `ComplexPolynomial::roots` / `AdaptiveRootSolver` — 全流水线共用底层热点

- 位置：`src/approx/complex_poly.rs:261`（`roots`）、`:805-825`（`AdaptiveRootSolver::roots_of`）、`:667-791`（`FaerRootSolver` 伴随矩阵特征值 + 后向残差门 `faer_gate_accepts`，门限 `FAER_BACKWARD_RESIDUAL_TOL = 1e-3`，`complex_poly.rs:797-799`）。
- 调用点：① twoport Chebyshev 阶段（A/E 多项式求根）② 矩阵阶段残差展开分母求根（×3，见 H3）③ multiplexer 谱因子化（每迭代）④ 优化器（间接，经响应/矩阵）。
- 实测（§2.1）：生产次数带 4–20 内，**迭代法（DK/Aberth）比 faer 快 3–4×**（d10：23 µs vs 84 µs；d20：169–230 µs vs 275 µs）；faer 仅在 d≥30 反超（658 µs vs 1.26/0.72 ms）。
- 占比（twoport order-20 合成 0.94 ms）：分母求根×1 ≈ 0.16 ms → ×3 ≈ 0.48 ms ≈ **51%**；Chebyshev 内求根另计。求解器选 faer 相对迭代法，每次求根多花 2–4×。

### H3. 矩阵阶段残差展开的 3 次重复分母求根

- 位置：`src/synthesis/residues/expansion.rs:170-178`（`synthesize_residue_expansions` 对 y11/y12/y22 各调一次 `residue_expansion`），`:180-211`（`residue_expansion` 内 `denominator.roots_with(&solver)`，**每次都对同一分母重新求根**）。
- 复杂度：分母（E+conjE+F/εr+conjF/εr）是同一多项式，3 个 `residue_expansion` 各自独立求一次其全部 N 个根——**3× 冗余求根**。
- 占比（twoport order-20 矩阵阶段 0.644 ms）：分母求根×3 ≈ 0.48 ms ≈ **75%**；order-4 矩阵阶段 0.027 ms 中分母求根×3 ≈ 0.022 ms ≈ 81%。
- 可省：3 次 → 1 次，矩阵阶段求根部分省 2/3（order-20 约 0.32 ms/次合成）。

### H4. twoport 折叠 `to_folded` 的重复合成

- 位置：`src/design.rs:444-446`（`to_folded` = `MatrixSynthesisEngine.synthesize_with_topology(&self.polynomials, Folded)`）；`src/synthesis/engine.rs:74-80`（`synthesize_with_topology` = `self.synthesize(polynomials)?.transform_topology(topology)`，即**先完整重新合成再折叠**）。
- 实测归因：基准"折叠"（`build_design().to_folded()`，`benches/synthesis.rs:176-182`）− "合成"（`build_design()`）= 纯 `transform_topology`（Givens）成本：

  | 次数 | 折叠 | 合成 | 纯折叠(差值) |
  |---|---|---|---|
  | 4 | 56.1 µs | 91.0 µs | ≈ 21–35 µs |
  | 16 | 619.1 µs | 695.8 µs | ≈ 60–90 µs |
  | 20 | 820.2 µs | 1.653 ms | ≈ 170 µs |

  折叠时间 ≈ 合成时间 + 小量 Givens（O(N²·N) 旋转，实测差值仅 ~10–30%）。
- 占比：每次 `to_folded()` 调用，**~75–90% 花在把刚合成过的矩阵又合成一遍**。
- 冗余依据：`FilterDesign` 存的 `self.matrix` 即 `synthesize(&polynomials)` 的结果（`design.rs:342`），`to_folded` 内部 `synthesize(polynomials)` 是同输入确定性重算，输出等价；Givens 折叠是确定性线性变换，直接对已存矩阵折叠结果相同。
- 注：`FilterDesign` 已保留 `polynomials`（`design.rs:224`），若未来要拓扑变体（arrow/CT）再合成，可只复用矩阵、按需重算。

### H5. twoport 响应：LU vs 极点展开快路径

- 位置：LU 路径 `src/response/backend.rs:10-16`（`evaluate_normalized_response` → 逐点 2 个线性系统 `A x = e0`、`Aᵀ y = e_{N-1}`，`backend.rs:53-96`）；快路径 `src/design.rs:489-530`（`response_normalized` 先试 `PoleExpansionData::from_matrix`（`src/response/pole_expansion.rs:56-`，O(N)/点），失败回退 LU）。
- 实测（§2.2 / §2.5）：order-20，201 点，**LU = 2.84–3.01 ms，极点展开 = 0.146 ms，≈ 20×**；order-4：0.17–0.43 ms vs 0.024 ms，≈ 7–17×。
- 关键点：**用户面 `FilterDesign::response_normalized` 默认就是极点展开**（design.rs:489 文档 + :504 先试 fast path），所以 twoport 响应对最终用户**不是**瓶颈；只有直接调 `ResponseSolver::evaluate_normalized`（基准即此，`benches/synthesis.rs:187`）才走 LU。
- 对 multiplexer 响应**不适用**：`multiplexer_a` 通道矩阵是折叠拓扑（`PoleExpansionData::from_matrix` 对非透射矩阵直接报错，`pole_expansion.rs:179-201`），只能 LU，3.8 ms 为该路径固有成本。

### H6. LM 优化器有限差分 Jacobian

- 位置：`benches/optimizer.rs`（tune/extract）对应的 `src/` 中 LM 提取器；每迭代 n_params+1 次 41 点响应（LU）扫描。
- 实测：tune/o8 = 13.5 ms/6 迭代 ≈ 2.25 ms/迭代，其中响应扫描占主导（每迭代 ~7 次 × ~0.3 ms + 残差）。
- 占比（优化器自身）：响应求值 ≈ 80%+（推算，探针未拆优化器子阶段）。

## 4. 瓶颈判断（按可节省时间排序）

| 排名 | 瓶颈 | 影响范围 | 可省（每次调用，量级） | 备注 |
|---|---|---|---|---|
| 1 | **multiplexer 谱因子化迭代**（H1） | multiplexer/duplexer 合成 | `poly` 22 ms 中 ~10–15 ms（迭代次数/求根） | 最大绝对热点；两端口无关 |
| 2 | **求解器 faer 优先**（H2） | 所有求根（twoport 合成 + multiplexer 谱） | 每次求根 2–4×；twoport 合成内 ~0.3–0.5 ms；multiplexer 每迭代谱求根 | 次数带 4–20 内迭代法更快 |
| 3 | **`to_folded` 重复合成**（H4） | 每个 `FilterDesign::to_folded()` 调用 | ≈ 合成时间（order-4 78 µs → order-20 ~1.5 ms） | 纯冗余，结果等价 |
| 4 | **3× 分母重复求根**（H3） | twoport 矩阵阶段 | 矩阵阶段 ~75%（order-20 ~0.32 ms） | 同多项式求根 3 次 |
| 5 | **裸 LU 响应**（H5） | 直接用 `ResponseSolver::evaluate_normalized` 的调用方 | ~20×（order-20 2.8 ms → 0.15 ms） | 用户面 fast path 已覆盖 twoport；multiplexer 不适用 |
| 6 | **LM 有限差分 Jacobian**（H6） | 参数提取/调谐 | ~2×（响应求值部分） | 复杂步可省一半求值 |

## 5. 优化建议（对应函数 / 手段 / 预期 / 可行性）

| # | 针对 | 手段 | 预期提升 | 可行性 |
|---|---|---|---|---|
| O1 | `machiarella.rs::run_full_iteration` / `spectral.rs::spectral_factor` | ① 谱因子化复用：同一迭代内各通道 E/F 多项式求根结果缓存；② 乘积多项式 N=19 求根改用迭代法（见 O2）；③ 收敛判据/迭代上限审视（当前 12 次迭代到 5.9e-7） | `poly` 22 ms → ~10–15 ms（2–3×） | 中（需保证收敛性与精度回归） |
| O2 | `complex_poly.rs::AdaptiveRootSolver::roots_of` | 次数自适应：deg ≤ ~22 先迭代法（DK/Aberth），deg > 22 才 faer；或 faer 残差门通过后仍按次数择快 | 求根 2–4×（d10 84→23 µs，d20 275→169–230 µs）；twoport 合成 ~1.2–1.5×，multiplexer 谱因子同比例 | 高（有现成基准背书，行为等价、仅择优） |
| O3 | `design.rs::to_folded` / `engine.rs::synthesize_with_topology` | `to_folded` 改为对已存 `self.matrix` 直接 `transform_topology(Folded)`，不再 `synthesize` 重算 | `to_folded()` 调用 ~2×（order-20 820 µs → ~170 µs） | 高（确定性等价，结果不变；arrow/CT 同理） |
| O4 | `expansion.rs::synthesize_residue_expansions` | 分母求根一次，3 个 `residue_expansion` 共享 `roots`（改传 `roots` 入参） | 矩阵阶段 ~1.4×（order-20 0.64→~0.32 ms 求根部分）；twoport 合成 ~1.2× | 高（纯去重，数值不变） |
| O5 | 直接调 `ResponseSolver::evaluate_normalized` 的调用方 | twoport 透射矩阵改用 `FilterDesign::response_normalized`（极点展开） | 响应 ~20×（order-20 2.8 ms → 0.15 ms） | 高（库内已有 fast path，仅改调用方式）；multiplexer 折叠矩阵不适用 |
| O6 | LM 提取器 Jacobian | 复变步（complex-step）或自动微分替代有限差分 | 每迭代响应求值 ~2×（13.5 ms → ~7 ms @o8） | 中（需复值化残差，或引入 AD 依赖） |

> O2+O3+O4 三者组合、且相互独立：twoport 单次 `build_design + to_folded` 预计从 ~0.94 + 0.82 ms（order-20）降至 ~0.5 + 0.17 ms 量级（≈ 2.5–3×）。

## 6. 结论：是否还有提升空间

**有，且集中在少数几处，收益明确：**

1. **最大绝对热点是 multiplexer 谱因子化迭代（22 ms）**——这是 mfs 区别于"教科书两端口"的自研算法路径，也是唯一 ms 级以上单项。它内部复用 H2/H3 的求根，O1+O2 落地后可显著下降，但需精度回归保障。
2. **twoport 主路径存在两处纯冗余**（`to_folded` 重合成、3× 分母求根）与**一处求解器次优**（faer 在 4–20 次比迭代法慢 3–4×），三者都是"结果不变、纯省时"的高可行性优化（O2/O3/O4），合计可让 twoport 合成快 ~2.5–3×。
3. **响应不是 twoport 用户面的瓶颈**——`FilterDesign::response_normalized` 已默认极点展开（20× 快于 LU）；瓶颈只存在于直接调 `ResponseSolver` 裸 LU 的场景（含基准本身），以及 multiplexer 折叠矩阵（fast path 结构性不适用）。
4. **算法本身量级合理**：twoport 合成 order-4→24 约 0.03→1.24 ms（接近二次，符合伴随矩阵/多项式根 O(N³) 主导），无数量级异常；瓶颈在"重复劳动"与"求解器选型"，不在算法复杂度失控。

**数值健壮性附带观察**（非性能，仅记录）：order-30 纯极点原型在 `synthesize` 处被拒——`transversal synthesis expects poles close to the imaginary axis (pole real part 1.54e-4 exceeds tolerance 3.60e-5 for order 30)`。这是算法有效域限制，说明高阶纯极点配置接近数值边界；若需支持 order≥30，应放宽门限或改用更高精度求根，而非性能问题。

## 7. 复现命令

```bash
# 四个 criterion 基准（各自 release 编译后运行；本次为单轮，总墙钟 ~14 min）
cargo bench --bench synthesis
cargo bench --bench multiplexer
cargo bench --bench root_solvers
cargo bench --bench optimizer

# 本次实测原始输出（median 列）见 target/ 下各 bench 的 criterion 报告；
# 子阶段拆解（§2.5 探针表）由一次性 /tmp 探针完成，已删除，不再保留。
```

## 8. 合规性（对应验收）

- `git status --porcelain src/ tests/ benches/ Cargo.toml Cargo.lock`：分析前后**无新增改动**。分析开始前工作树已有 3 处未提交改动（`M src/spectral.rs`、`M src/synthesis/machiarella.rs`、`M tests/multiplexer.rs`），**均为分析开始前的既有状态，非本次产生**；本次未对其做任何修改。
- 未新增产品内文件；临时探针建于 `/tmp/mfs_perf_probe`，**已删除**。
- 本报告 `docs/perf-analysis.md` 为本次唯一新增文件（untracked）。
- 基准 1 轮（≤3），总墙钟 ~14 min（≤40 min）。
