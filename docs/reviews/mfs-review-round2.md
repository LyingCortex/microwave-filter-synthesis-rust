# MFS 第二轮深度 Review 报告(复核修订版)

- 审查日期:2026-10-12(复核修订)
- 审查对象:HEAD `7dd617a`(master),工作树干净,rustc/cargo 1.95.0
- 审查方式:只审查,不修改源码/文档/依赖/发布记录;反例在仓库外 `/tmp/mfs_probe` 独立构建
- **本版为复核修订版**:对初稿逐项二次复核(调用上下文/已有保护/数学约定/测试覆盖),修正了初稿的三处误判,并新增一条初稿遗漏的发现(N1)。修订处见各条 "复核" 标注与 "总评"。

## 一、审查范围与项目理解

对象:mfs — Rust 微波耦合谐振滤波器综合库,实现 Cameron 广义 Chebyshev 耦合矩阵综合。
审查维度:A 数学约定 / B 综合流水线 / C 数值稳定性 / D 拓扑与响应 / E 高级功能 / F 测试独立性 / G 文档与发布。

### 基线核对(用户基线 vs 实际)

| 项目 | 用户基线 | 实际 | 一致性 |
|---|---|---|---|
| HEAD | 7dd617a | 7dd617adfa8a… | ✅ 一致 |
| 工作树 | 干净 | `git status --porcelain` 空 | ✅ 一致 |
| 工具链 | — | rustc/cargo 1.95.0 | ✅ |
| 版本声明 | "0.6.2-unreleased" | Cargo.toml `0.6.2`(无 -unreleased) | ⚠️ 见 P3-4 |
| cargo test | 435 通过 | 默认 **419** 通过;`--features cli,python` **435** 通过 | ⚠️ 435 是 feature 档 |
| Clippy 两档 | 0 警告 | 两档均 0 警告(`-D warnings`) | ✅ 一致 |

理解要点:核心数学链 `Cameron 递推 → F/P → find_eps(ε,ε_r) → A → E → 导纳四元 [Y] → 极点展开/M_SL → 折叠/CT/CQ → 响应`。

### 核心数学结论(已验证为正确,非缺陷)

- ✅ Cameron 递推 `U_new[i] = U[i+1] + temp·(1−w²)·V[i-1]`,V 更新项 `temp·U[i]`,与 1999 文献一致。
- ✅ E 公式 `E(w)=F(w)/eps_r + j·P(w)/eps`(`generalized_chebyshev_helpers.rs:408`,`j` 因子在位,奇数阶关键)。
- ✅ 导纳 paraconjugate 含全局 `(−1)^N`(`admittance.rs:61-69`),奇数阶防丢阶。
- ✅ 响应正确性:mixed 设计(N mod 4 = 0 与 2,order 4/6/8/10/12)在 ω=±1 **精确命中 RL spec(0.000 dB 偏差)**;all-finite 可行设计同样命中;功率守恒 |S11|²+|S21|²−1 ≤ 1.8e-15;order-30 极点展开快路径正常。
- ✅ **(复核新增)** all-finite 全规范成功路径实测正确:`twoport_4_4`(4-4, 4 有限零点)与 `6-6`(6-6)在 ω=1 命中 spec(|S21|=0.9968, −0.027 dB),功率守恒 ≤ 9e-16;order-1 全有限(z=2, RL=1)经 M_SL 闭式分支合成后 w=2 处 |S21|→0(零点正确)、功率守恒 ≤ 1e-15。

## 二、已确认缺陷(按 P0/P1/P2/P3)

无 P0。

### P1-1 — all-finite 不可行 spec 返回误导性硬错误,且该域限制完全未文档化

- 文件/位置:`src/approx/generalized_chebyshev_helpers.rs:276`(`eps_r = eps/√(eps²−1)`)、`src/approx/complex_poly.rs:95`(错误实际抛出点)、`docs/algorithms.md:119-120`(未标注有效域)。
  - **复核更正**:初稿把错误抛出点写成 `complex_poly.rs:96`/`85`;实际有限性检查在 `:91-97`,`PreconditionViolation` 抛出在 `:95`。
- 触发条件/最小反例:
  ```rust
  FilterDesign::prototype(4, 20.0).zeros([1.05,-1.05,1.1,-1.1]).synthesize()
  // → Err(MfsError::PreconditionViolation("complex polynomial coefficients must be finite"))
  FilterDesign::prototype(4, 34.0).zeros([1.5,-1.5,2.0,-2.0]).synthesize() // 同样触发
  FilterDesign::prototype(4, 20.0).zeros([1.5, 2.0,-1.5,-2.0]).synthesize() // 对照组:OK(可行)
  ```
  即 all-finite 零点(finite_count==order)且 ε<1:近带边零点(|z|≈1)或高 RL 均落入此域。
- 数学/代码证据:`find_eps` 中 `eps = |P(j)|/|F(j)|/√rl`;all-finite 时 `eps_r = eps/√(eps²−1)`(line 276),ε<1 时 `√(eps²−1)=NaN` → `1/eps_r=NaN` → `build_e_polynomial_stage` 的 `scale(1/eps_r)`(line 410)产生 NaN 系数 → 被 `complex_poly.rs:95` 的通用有限性检查拦截,抛出与根因无关的 "coefficients must be finite"。实测 ε 边界:order-4 RL20,`zeros=[1.5,2.0]`→ε=4.80(OK)、`[1.2,1.6]`→ε=1.42(OK)、`[1.05,1.1]`→ε<1(INFEASIBLE);RL 扫描 RL≤33 OK(eps_r 1.02→2.83),RL≥34 FAIL。
- 预期 vs 实际:预期返回带域信息的 typed error(如 `InfeasibleSpec("all-finite: eps<1, eps_r 无实值,建议改用 mixed 零点")`);实际是误导性硬错误,无指引告知可行替代(mixed 零点、eps_r 强制=1,同样近边零点可合成,已验证)。
- 影响/严重度:all-finite 是合法交付路径;可行时(ε≥1)数学完全正确(实测命中 spec),所以是错误信息质量+可发现性缺陷,非算法错误。但 all-finite 是文档支持的配置,用户撞墙后无指引 → P1(高影响、低概率,无文档/错误告知)。
- 最小修复方向:`find_eps`(line 273-277) all-finite 分支加 `if eps < 1.0 { return Err(含 eps 值与替代建议的 typed error) }`;`docs/algorithms.md:120` 与 `docs/architecture.md` Known Limitations 标注 all-finite 有效域 ε≥1。
- 回归验证建议:integration test 断言 (a) all-finite 可行 spec(远零点低 RL)成功且命中 RL spec;(b) all-finite 不可行 spec 返回特定 typed error 而非泛化 "coefficients must be finite"。
- 置信度:高(多次独立探针复现,根因链完整)。

### P2-1 — `PoleExpansionData::from_matrix` 对非横向(折叠)矩阵静默返回错误 S 参数(新发现,初稿遗漏)

- 文件/位置:`src/response/pole_expansion.rs:87-127`(`from_matrix`),读取循环 `:96-108`;文档 `:84-86`。
- 触发条件/最小反例:
  ```rust
  let d = FilterDesign::prototype(4, 22.0).zeros([1.3217, 1.8082]).synthesize()?;
  let folded = transform_matrix(d.matrix(), TopologyKind::Folded)?.matrix;   // 折叠矩阵(含谐振器间耦合 M[1,2]=0.960, M[2,3]=0.286)
  let ped  = PoleExpansionData::from_matrix(&folded)?;                      // 无任何错误/警告,直接成功
  let fast = ped.evaluate_normalized(&grid)?;
  let lu   = ResponseSolver{}.evaluate_normalized(&folded, &grid)?;         // LU 基准(正确)
  // w=0:   LU |S21|=0.999949   fast |S21|=0.000000 (-300 dB)   relErr=1.00
  // w=1:   LU |S21|=0.996840   fast |S21|=0.000000 (-300 dB)   relErr=1.00
  // 全网格 max |ΔS21| ≈ 1.00(100% 相对误差)
  ```
- 数学/代码证据:`from_matrix` 的模型假设谐振器块 D 是**对角阵**(Schur 分块 `S[0,1]=M_SL − Σ M_Sk·M_kL/(jω+B_k)`,`:157-161`),只读取 `M[0,k]`(源耦合)、`M[k,k]`(对角)、`M[k,N+1]`(载耦合)、`M[0,N+1]`(`:97-99,111`),**完全忽略谐振器间耦合 `M[k,k±1]`**。折叠矩阵恰有大量非零 `M[k,k±1]`(实测 0.960/0.286),因此部分分式模型对折叠网络不成立。`from_matrix` 不检查 `matrix.topology()`,也不检查 `M[k,k±1]==0`,对折叠/箭头矩阵一律静默接受。
- 预期 vs 实际:文档(`:84-86`)明言 "a non-transversal (folded/arrow) matrix … **must be transformed to transversal form first**";预期对非横向输入返回 typed error 或至少拒绝;实际是静默计算出一个**物理上错误的网络**(忽略所有谐振器间耦合)的响应,S21 在中心频率处从正确的 ≈1 变为 0(−300 dB),无任何错误、无警告。
- 影响/严重度:
  - **主路径受保护**:`FilterDesign::response_normalized`(`design.rs:419-449`)与 `response_bandpass`(`:527-558`)在采用快路径前,于 5 个采样点上对 LU 做 `max_diff < 1e-3` 交叉校验,折叠矩阵的 100% 偏差必然触发回退到 LU。故**经 `FilterDesign` 主 API 不会返回错误结果**。
  - **直接调用不受保护**:`PoleExpansionData` 是 `pub`(lib.rs 导出 `response`),`from_matrix`/`evaluate_normalized` 均 `pub`;折叠矩阵是库自身 `to_folded`/`to_topology` 的正常产物。直接调用者绕过 `FilterDesign` 的保护(例如已持有折叠矩阵、追求 O(N) 快路径)会得到**静默 100% 错误**的 S21。
  - 综合:P2(公开 API 在文档声明的输入类上静默返回错误数值;主路径有护栏,故非 P1)。
- 最小修复方向:在 `from_matrix` 入口校验横向性 —— 扫描 `|M[k,k±1]| > tol`(或 `matrix.topology() != MatrixTopology::Transversal`),命中则 `return Err(MfsError::PreconditionViolation("PoleExpansionData requires a transversal matrix; transform folded/arrow matrices first"))`。把静默错误转为 typed error,与文档声明一致。
- 回归验证建议:integration test 断言对折叠矩阵调用 `PoleExpansionData::from_matrix` 返回 typed error(或,若保留静默语义,则断言返回的响应与 LU 偏差 >1e-3 时不采用);现有 `folded_response_invariance` 已覆盖主路径回退,补一条直接调用路径的用例即可。
- 置信度:高(外部探针实测 100% 偏差;主路径护栏经源码逐行确认)。

### P3-1 — `admittance.rs` 首项抵消守卫为内部死代码,注释描述失效机制(初稿 P2-1,复核降级为 P3)

- 文件/位置:`src/synthesis/residues/admittance.rs:40-54`(守卫本体+注释 40-44);`docs/algorithms.md:157-159`(描述为活跃机制)。
- 触发条件:所有库内部调用路径上**不触发**;经公共 `pub` 字段可触发(见下)。
- **复核更正**(初稿两处算术/措辞错误):
  1. 初稿写 `f_leading = 1/eps_r ≥ 1` —— 错。`eps_r ≥ 1`,故 `f_leading = 1/eps_r ≤ 1`(混合零点时 =1,全有限时 <1)。
  2. 初稿写 "永不触发" —— 过强。该守卫经**公共 API 可触发**:`GeneralizedChebyshevData.e_s` 是 `pub` 字段(`ComplexPolynomial.coefficients` 亦 `pub`),用户可构造非 monic 的 `e_s`(leading=−0.978)喂给 `synthesize_admittance_polynomials`,此时 `e_leading.re<−1e-15` 且 `sum_leading≈0`,守卫**会触发**并翻转 E。故它是"针对未文档误用路径的防御",非纯死代码。
- 代码/数学证据:`build_e_polynomial_stage` `:425` 对 E 做 `normalize_to_monic`,`:468` 存入 `generalized.e_s` → 内部路径消费的 `e_leading` 恒为 +1+0j(实测 order 4/6/8 均 1+0i);`f_over_eps_r` 由 monic `f_s`(`:452`)缩放 `1/eps_r` → `f_leading=1/eps_r ∈ (0,1]`。内部 `e_leading+f_leading = 1+1/eps_r ∈ [2,3]`,`norm_sqr ≥ 4 ≫ 1e-20`,且 `e_leading.re=1 ≮ −1e-15`,两条件在内部路径均不成立 → 守卫在库自身综合流程中是 no-op。注释 `:41-44` 声称 "N mod 4==2 时 E_lead = −1 … exact cancellation" —— 这是 0.6.1 之前(E 未 re-monic)的状态;`:425` re-monic 之后 E 恒 +1,与 F/eps_r(+1/eps_r>0)相加而非相消。正确性无损:实测 order-6/10 mixed 设计 ω=±1 精确命中 spec。
- 影响/严重度:内部无行为影响(死 no-op)+ 误导性注释/文档,维护者可能被引导去"修复"不存在的抵消;但守卫本身是对公共字段误用的有效防御,可保留 → **P3(文档/维护),非 P2**。
- 最小修复方向:保留守卫;重写 `:41-44` 注释,说明其实际角色是"防御经公共 `pub` 字段构造的非 monic `e_s`",并删除 "N mod 4==2 → E_lead=−1" 这一已失效陈述;同步更正 `algorithms.md:157-159`(勿描述为活跃的抵消消解机制)。
- 回归验证建议:守卫保留则无需行为回归;若删除,跑 `tests/`(`synthesis_regression`、`folded_response_invariance`、order 2/6/10/14 用例)确认 bit-exact 不变。
- 置信度:高(静态推导 + 实测 E 首项 + 公共字段可达性三重印证)。

### P3-2 — `docs/internal/todo.md` 头部陈旧(HEAD 漂移)

- 位置:`docs/internal/todo.md:11`("master head `e0c2699`")与 `:1-3`("0.6.2-unreleased")。
- 证据:实际 HEAD = `7dd617a`(e0c2699 之后;`git merge-base` 确认 e0c2699 是祖先)。todo.md 的 "435 passed" 对应 `--features cli,python`,而头部语境是普通 `cargo test`(实为 419)。
- 影响/严重度:文档时效性,误导对"当前已验证 head"的判断。P3。
- 最小修复方向:头部 head 更新为 `7dd617a`,测试计数注明 feature 档。
- 置信度:高。

### P3-3 — `algorithms.md` 文件引用指向不存在的 `helpers.rs`

- 位置:`docs/algorithms.md`(11 处 `helpers.rs:NNN`)。
- 证据:`src/approx/` 下无 `helpers.rs`;真实实现文件是 `generalized_chebyshev_helpers.rs`,`helpers` 只是 `mod.rs:35` 的 re-export 模块。行号(214/247/385/409-411)实际对应 `generalized_chebyshev_helpers.rs`,但文件名引用错误。
- 影响/严重度:文档可导航性。P3。
- 最小修复方向:批量替换 `helpers.rs` → `generalized_chebyshev_helpers.rs`(或统一用模块路径 `approx::helpers`)。
- 置信度:高。

### P3-4 — 版本/发布状态不一致

- 位置:`Cargo.toml:3`(`0.6.2`)、`pyproject.toml:7`(`0.6.2`)、`CHANGELOG.md:89`(已定版 `[0.6.2] — 2026-10-11`)与 `:3`(非空 `[Unreleased]`)、git tag(无 `v0.6.2`,最后 tag `v0.6.1`)。
- 证据:CHANGELOG 同时存在已定版的 `[0.6.2]`(10-11 发布)与非空 `[Unreleased]`(flip_sign 去重等 post-release 改动)。Cargo.toml 写 `0.6.2`(无 -unreleased)却无对应 tag,而 todo.md 又称 "0.6.2-unreleased"。
- 影响/严重度:发布状态语义混乱;CI `version-sync` 只校验 Cargo↔pyproject(两者都 0.6.2,通过),无法捕获此漂移。P3。
- 最小修复方向:`[Unreleased]` 非空且 0.6.2 已发布 → bump 至 `0.6.3(-unreleased)`,或后移 0.6.2 定版边界。
- 置信度:高。

## 三、待验证风险(与已确认缺陷分开)

- **R1**(保留):all-finite ε 接近 1 的数值条件数。ε→1⁺ 时 `√(eps²−1)→0`、`eps_r→∞`、`1/eps_r→0`,求根/monic 归一条件数恶化;本轮未做 ε→1⁺ 极限扫描。待验证,非确认缺陷。
- **R2**(升级:根因已确认):`multi_band.rs:416` 复用 `find_eps`,其 all-finite 分支受同一 ε<1 域限制。本轮确认多频带 all-finite(5 频带,4 有限零点)在 RL=20 **可行成功**,与单带共享同一 `find_eps:276` 根因;可行域可用,不可行域(近带边零点)理论同 P1-1,未逐频带实测。根因确认,待补 typed error。
- **R3**(解决):奇数阶 all-finite 本轮已实测 order-1(z=2)与 order-3(3 有限零点)均合成成功(RL=10/40/42.4 命中 spec、零点正确、功率守恒)。奇数阶 all-finite 工作正常,原 "未覆盖" 风险解除。

## 四、测试与独立验证缺口

1. **(复核更正)** 初稿 "all-finite 成功路径零测试" **有误**。`tests/folded_response_invariance.rs:99`(twoport_4_4, 4-4 全有限)与 `:105`(6-6 全有限)是真实 all-finite 成功用例,且通过(命中 spec、折叠保持响应)。
   - 真实缺口是:**all-finite 不可行路径的 typed error 无回归测试** —— 全仓无任何测试断言 P1-1 那种 ε<1 spec 返回特定错误(当前会掉进泛化 "coefficients must be finite");order-1 不可行角落亦无覆盖。
2. 测试计数口径:"435" 仅 `--features cli,python` 成立;默认 `cargo test` 为 419。CI `test` job 跑默认(419),`cli` job 跑 cli_integration。基线陈述易被误读为默认档。
3. **N1 直接调用路径无测试**:`folded_response_invariance` 经 `FilterDesign::response_normalized`(有 LU 回退护栏)覆盖折叠响应,但无任何测试断言 `PoleExpansionData::from_matrix` 对非横向矩阵的拒绝/行为。
4. order-1 有限零点可行域(RL<~7 dB)与不可行域(RL≥~8)的边界无专项测试;P1-1 的 order-1 角落未覆盖。

## 五、文档及发布状态冲突

P3-2(todo.md HEAD 漂移)、P3-3(algorithms.md 文件名)、P3-4(版本/发布状态),以及 P1-1 / P3-1 的文档部分(algorithms.md:119-120 未标 all-finite 域;:157-159 描述 P3-1 守卫为活跃抵消机制)。

## 六、实际执行的检查及结果

| 检查 | 命令 | 结果 |
|---|---|---|
| 基线 | `git log -1`、`git status --porcelain`、`rustc/cargo --version` | HEAD 7dd617a,干净,1.95.0 |
| 默认测试 | `cargo test --all-targets` | 419 通过,0 失败 |
| Feature 测试 | `cargo test --all-targets --features cli,python` | 435 通过,0 失败 |
| Clippy 两档 | `cargo clippy --all-targets -- -D warnings`(default + cli,python) | 各档 0 警告 |
| 响应正确性探针 | 外部 `mfs_probe` crate(`/tmp/mfs_probe`) | mixed order4/6/8/10/12 在 ω=±1 精确命中 RL spec(0.000 dB);功率守恒 ≤1.8e-15 |
| all-finite 成功探针 | 同上 | twoport_4_4(4-4)/6-6(6-6)/order-1(z=2,RL1)命中 spec,功率守恒 ≤1e-15 |
| all-finite 不可行探针 | 同上 | ε 边界、ε<1→NaN→"coeffs must be finite" 复现;order-1 边界 RL*≈7 dB(z=2) |
| N1 折叠矩阵探针 | 同上 | `from_matrix`(折叠)静默成功;fast vs LU 全网格 maxΔS21≈1.00(100%);主路径 LU 回退经源码确认 |
| M_SL 分支探针 | 同上 | order-1 all-finite 矩阵 M[0,2]=0.12934 == `eps·(eps_r−1)/eps_r`;闭式分支可达且正确 |
| E 首项探针 | 同上 | order 4/6/8 `E.leading=1+0i`(内部守卫为 no-op) |
| 版本/tag | `git tag`、`grep version` | 无 v0.6.2 tag,最后 v0.6.1;Cargo/pyproject 均 0.6.2 |

执行限制:WSL2/NTFS 文件系统,全量 `cargo test` 单轮 ~800s;探针在仓库外独立构建,未触碰被审查工作树。

## 七、未覆盖范围

- 未运行 `cargo bench`(仅 CI 编译门禁);未构建/运行 Python wheel 与 Python 端集成测试(仅跑 Rust `--features cli,python` 测试)。
- 未跨 OS(CI matrix windows/macos)实测,仅 local WSL2。
- 未深入 `multiplexer.rs`/`duplexer.rs` 全谱路径的数学复核(0.6.1 已修 P1,本轮仅抽查响应正确性)。
- 未做 ε→1⁺ 极限(R1);多频带逐频带不可行域(R2)未逐一实测;order-1 不可行角落未覆盖。
- 未审计 `optimize.rs` LM 收敛的数值细节(仅确认结构与 scale-aware FD Jacobian 在位)。

## 八、最优先处理的三项

1. **P2-1(新发现,最高)**:`PoleExpansionData::from_matrix` 入口加横向性校验,对折叠/箭头矩阵返回 typed error。这是本轮唯一"公开 API 在文档声明的输入类上静默返回 100% 错误数值"的缺陷;主路径虽有 LU 回退护栏,但直接调用(已持有折叠矩阵、追求 O(N) 快路径)无任何保护。修复最小(入口扫描 `M[k,k±1]`/`topology()`),回归补一条直接调用拒绝用例。
2. **P1-1**:为 all-finite 不可行 spec 提供带域信息的 typed error(`find_eps` line 273-277 加 `eps<1` 前置检查),并文档化 all-finite 有效域 ε≥1 + 补不可行路径 typed-error 回归测试。唯一让合法用户撞墙且无指引的缺陷。
3. **P3 文档/发布批次**:P3-1(改写 `admittance.rs` 死守卫注释 + 更正 algorithms.md:157-159)、P3-3(algorithms.md 文件名)、P3-4(版本/发布状态理顺)合并处理,避免维护者被失效注释误导。

## 总评(复核修订说明)

本轮核心数学未发现错误——Cameron 递推、E 公式 `j` 因子、paraconjugate `(−1)^N`、mixed/all-finite 响应正确性均经独立探针验证(命中 spec 至 0.000 dB、功率守恒 ~1e-15)。

**对初稿的三处复核更正**:
1. 初稿 §四.1 "all-finite 成功路径零测试" 有误——`twoport_4_4`/`6-6` 是真实 all-finite 成功用例且在测;真实缺口是**不可行路径 typed error 无测试**。
2. 初稿 P2-1 把 `admittance.rs` 守卫定性为"死代码 + 永不触发" 过强——守卫经公共 `pub` 字段可触发(防御非 monic `e_s`),故降级为 P3(文档/维护);并更正初稿 `1/eps_r ≥ 1` 的算术错误(应为 ≤1)。
3. 初稿 §四.3 "order-1 all-finite 必 ε<1" 有误——order-1 all-finite 在低 RL(RL<~7 dB)可行,经 M_SL 闭式分支正确合成(功率守恒、零点正确);不可行仅在 RL≥~8 dB。

**新发现**:P2-1 `PoleExpansionData::from_matrix` 对非横向(折叠)矩阵静默返回 100% 错误 S21(主路径有 LU 回退护栏,直接调用无保护)。

已确认缺陷集中于:非横向矩阵快路径的输入校验缺失(P2-1)、all-finite 不可行域的误导性错误 + 未文档化(P1-1)、一处内部死守卫的失效注释(P3-1)与文档/发布状态漂移(P3-2/3/4)。验证边界明确:未证明整个库正确,仅覆盖上述已验证路径;R1-R3 列为待验证/已解决,非缺陷。
