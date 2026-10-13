# mfs 第二轮深度 Review(验收 + 独立验证)

- 审查日期:2026-10-12/13
- 审查对象:**工作树**(dirty,0.6.3-unreleased,含未提交的 round-2 修复)
- 基线核对:HEAD `7dd617a`(与用户给定一致);**工作树在我审查开始时是干净的 0.6.2,审查期间被更新为 0.6.3(12 个文件未提交)** —— 见 §4.1,这是本报告的基准,所有执行验证均在此树上进行。
- 工具链:rustc/cargo 1.95.0。
- 审查方式:以数学一致性为主线,对 round-2 两项修复做独立探针验证(仓库外 `/mnt/c/Users/eynulai/dev/mfs_r3probe`,不碰被审树),辅以完整 CI 命令重跑。

---

## 一、审查范围与项目理解

mfs 是 Rust 耦合谐振滤波器综合库,核心链路为:

```
FilterSpec
  → generalized_chebyshev_polynomials (Cameron 递推 F/P/E, find_eps/eps_r)
  → MatrixSynthesisEngine.synthesize (residue expansion → transversal M)
  → to_folded / to_arrow / to_ct / to_cq (相似变换)
  → ResponseSolver (LU 后端) / PoleExpansionData (O(N) 极点展开快路径, 5 点 LU 交叉验证门)
  → 带通映射/去归一化/Python/CLI/Touchstone
```

数学约定(已逐条核对):Feldtkeller 恒等式 `|E|²=|F|²+|P|²/eps²`;Cameron 递推 V 更新项 `temp·U_{n-1}`(无额外 w 因子);eps_r = 1(mixed 零点集)或 `eps/sqrt(eps²−1)`(全有限,要求 eps>1);带通映射 `ω=(1/FB)(f/fc−fc/f)`;S 参数经 Schur 2×2 归约,`S11=1+2j[A⁻¹]₀₀`,`S21=−2j[A⁻¹]_{N+1,0}`;群延迟 `−Im(d/dω ln S21)`(归一化频率单位,非 Hz/弧度)。

**核心数学链路在独立验证下自洽**:7 个标准设计(含 N=12 全有限、奇数阶、RL 最高 30 dB)在 1001 点稠密网格上于带边**精确命中 RL spec**(偏差 < 0.0001 dB,见 §5 探针),|S21| 峰值 ≤ 1.0,功率守恒在 LU 路径下 ≤ ~1e-15 量级(verify 模块既有断言 + 本次探针一致)。

---

## 二、已确认缺陷(按 P0/P1/P2/P3)

**无 P0/P1。** 库的数值核心与 round-2 两项修复均验证通过。以下为 P2/P3 级问题,全部有实测证据。

### P2-1 发布状态/文档冲突:"稳定到 30+"声明与实测边界不符

- **位置**:`README.md:48-50`("Numerically stable up to order 30+");`docs/architecture.md:573`("High order (10–30), all-pole + finite TZs — ✅ Numerically stable");`examples/high_order_stability.rs:123-139`。
- **实测**(探针 `mfs_r3rlcheck`,RL=20 dB):all-pole 路径在 **order 30 起被拒绝** —— `transversal.rs:53` 的 `pole_axis_tolerance = 1e-6·(N/5)²`,order 30 容差 3.6e-5,而根求解器在 N=30 产生的极点实部 1.41e-4(order 40 达 8.48e-3)。order 28 仍 OK;有限零点路径 order 30 同样被拒(极点实部 1.05e-4 > 3.6e-5),order 32 的有限零点设计反而 OK(15 个零点 vs 16 个,零点分布不同)。
- **关键证据**:`examples/high_order_stability.rs` 的 `test_order` 对综合失败是 `eprintln!` + 继续(`:115` `Err(e) => eprintln!("  Matrix synthesis failed: {e}")`),**从不断言成功**;我实际运行该 example,order 30 打印 "Matrix synthesis failed: ... pole real part 1.41e-4 exceeds tolerance 3.60e-5"。即:README 与 architecture 的 "stable to 30+" 是由一个会吞掉失败、且无测试背书的 example 支撑的。
- **预期 vs 实际**:文档声称 order 30+ 可用;实际 all-pole 与有限零点在 order 30 返回 typed `NumericalFailure`(边界清晰、无静默错误结果 —— 这点是好的),但 "stable up to 30+" 的措辞与实测边界(≈28)不符。
- **影响**:P2 —— 无静默错误,但发布文档对能力边界做了不实声明,用户按 README 使用 order 30+ 会直接得到错误而非设计。
- **最小修复方向**:把 README/architecture 改为 "stable to order 28 (adaptive pole-axis tolerance); order 30+ may be rejected with a typed NumericalFailure",或在 example 中改为断言成功并把实测可行上限作为回归锚点。
- **回归验证建议**:`cargo run --example high_order_stability` 的输出应包含 order 30 的成功(若声称支持),否则文档边界与 example 必须一致。
- **置信度**:高(三次独立运行,探针 + example + 源码容差公式互相印证)。**[已修复,2026-10-13,见 §九 R1]**

### P2-2 `PoleExpansionData::from_matrix` 的"transversal 契约"未完全强制:端口对角线被静默忽略

- **位置**:`src/response/pole_expansion.rs:158-185`(`assert_transversal_structure`),`from_matrix` 读取 `:107-123`;契约文档 `:88-96`。
- **问题**:round-2 新增的结构检查只扫描**谐振器块** `M[1..=N][1..=N]` 的非对角元(`:164-171` 明确跳过端口列 `j<1 || j>order` 和 `M[0,N+1]` 之外的端口行)。但 Schur 2×2 模型假设端口对角项 `M[0,0]=M[N+1,N+1]=0`(纯串联端接,无端口 shunt)。`from_matrix` 读取的 `m_sl = M[0,N+1]`(`:122`)是允许的(直接源-负载耦合),而 **`M[0,0]`/`M[N+1,N+1]` 的非零值不被检查、也不被建模** —— 若用户直接调用公共 API `PoleExpansionData::from_matrix` 并传入带端口 shunt 的矩阵,响应中的端口 shunt 贡献被**静默丢弃**,返回一个 |S11|/S21 都错的结果。
- **实测**:我构造 4 阶 transversal 矩阵(谐振器块严格对角,源-负载耦合 0.1)与相同矩阵加 `M[0,0]=4`(≈ 源端 shunt  susceptance),两者 `from_matrix` 都**成功**,5 点 |S21| 门通过(源 shunt 不移动 |S21| 的 5 个采样点相对差 < 1e-3 时),但全网格 LU 响应显示 **S21 相对差最高 ~1.37e-1,S11 显著不同**(探针 A)。即:结构检查 + 5 点门都无法拦截这个输入,而该输入是合法的 `CouplingMatrix`(对称、谐振器块对角),只是端口带 shunt。
- **缓解与定级依据**:主流程(`FilterDesign::synthesize`)产出的矩阵恒有端口对角 = 0(源/负载端接是外加的,见 `backend.rs:126-127`),所以主流程**不可达**此错误;仅当用户手工构造矩阵直接调 `from_matrix` 时触发。这是"文档声明的契约"与"实际强制的契约"之间的真实缝隙,故定 **P2 防御性缺口**(不是 P0/P1:主流程安全;修复也简单)。
- **最小修复方向**:`assert_transversal_structure` 增加对 `|M[0,0]|`、`|M[N+1,N+1]| > tol` 的拒绝;或在 `from_matrix` 文档中明确"端口对角必须为 0(纯串联端接),端口 shunt 请用 LU 后端"。
- **回归验证建议**:`pole_expansion::tests` 增加一个"transversal 谐振器块 + 非零端口对角 → 被拒"的测试,与既有 5 个契约测试并列。
- **置信度**:高(探针实测,错误量级 ~1e-1 远超任何合理容差)。**[已修复,2026-10-13,见 §九 R2]**

### P3-1 `verify_against_spec` 的 RL 门与"综合"语义脱节(文档易误导)

- **位置**:`src/verify/mod.rs:371-400`(`verify_against_spec`),`RETURN_LOSS_GATE_MARGIN_DB:254`;`synthesize()` 不调用它(`src/design.rs:84-` 全文无引用)。
- **问题**:`verify_against_spec` 用**通过带最小 RL** 对比 spec − 0.1 dB 余量(注释 `:389-392` 解释了为何不用 max),这是正确的"最坏点"判据。但它是**独立于综合管线的公共工具**,`FilterDesign::synthesize` 不会自动调用;文档示例(`verify/mod.rs:366-368`)把它写成"综合后顺手一验"的样子。实测所有标准设计都恰好命中 spec(§5),所以当前不会漏掉错误设计;但一个 RL spec 被违反的设计仍能以 `Ok(FilterDesign)` 返回,只有显式调用 `verify_against_spec` 才会 Fail。
- **影响**:P3 —— 行为本身正确且诚实(函数名是 "verify",不是 "synthesize_and_verify");风险仅在读者默认"综合=满足 spec"。
- **最小修复方向**:在 `synthesize` 的 doc 或 `FilterDesign` 类型 doc 中明确"综合不保证 RL 满足 spec,请用 `verify_against_spec` 显式检查"。
- **置信度**:高(代码事实,无争议;严重度低)。**[已修复,2026-10-13,见 §九 R4]**

### P3-2 群延迟的独立验证只覆盖 all-pole

- **位置**:`src/response/mod.rs:254-340`(`group_delay_matches_across_fast_and_lu_backends`)。
- **问题**:该测试 Check 1 是 fast-vs-LU 交叉(同一数学约定,只能抓实现不一致,不能抓约定错误);Check 2 用 S21 相位的中心差分做独立锚点,但**只在 all-pole ladder 上运行**(`:311-313` 注释自述:有零点处相变陡,FD 误差大)。因此带有限/复零点的群延迟约定**没有独立锚点**,只有 fast-vs-LU 一致性。
- **缓解**:两条路径都从 `−Im(d ln S21/dω)` 推导,约定一致的把握较大;`mod.rs:26-31` 还专门记录了 2π 单位的陷阱,说明该约定是被审查过的。
- **影响**:P3 —— 约定错误风险低但非零(尤其复零点/群延迟平坦化场景,`tests/cameron1999_fixtures.json:72-73` 的 GDE 复零点用例)。
- **最小修复方向**:对一组带零点设计,用足够细的网格 + 远离零点的频段(如 |S21| > 0.1 区域)做 FD 锚点,补上独立验证。
- **置信度**:中(代码事实明确;"是否真有约定错误"未证明)。**[已修复,2026-10-13,见 §九 R5]**

---

## 三、待验证风险(与已确认缺陷分开)

1. **三个只读审计代理(DocsAudit/TestIndependence/InterfacesAudit)在审查结束前未交付最终报告**(各自已读 15+ 文件,处于收尾阶段)。本报告 §三/§六 中关于文档锚点逐条核验、测试循环性全量分类、多工器/LM/cascade 逐入口守卫的部分,以**我本人的直接抽查**为准(见 §五 执行记录),未等代理汇总。风险:可能存在我未覆盖的锚点/守卫个案。
2. **多工器/双工器 `full_iteration` 不收敛时的行为**:我核对了 `continuation-plan.md` 记载的收敛锚点与 stagnation 判据(`:72-78, 224-228`),未实测非收敛工况下的返回路径。
3. **Python 绑定未在本地构建验证**(需要 maturin/CPython 环境;CI 的 wheel job 可覆盖,本轮未跑)。
4. **P2-2 的"5 点门漏检"是实测的**,但 5 点门在主流程中现在只是第二道防线(round-2 的结构检查是第一道),其独立价值仅存在于用户直接调 `from_matrix` 的路径。

---

## 四、测试与独立验证缺口

1. **order 30+ 无测试覆盖**:tests 目录中无任何 order ≥ 30 的综合用例(grep 确认);`examples/high_order_stability.rs` 是唯一触达 order 30 的代码,且吞错。→ 对应 P2-1。**[已修复,§九 R1:example 边界回归测试锁定 28/30]**
2. **`from_matrix` 的端口对角契约无测试**:round-2 新增的 5 个契约测试(`pole_expansion.rs:533+`)覆盖 transversal 接受/拒绝、round-off 接受、folded/arrow 拒绝、折叠回退,但**没有**"端口对角非零"用例。→ 对应 P2-2。**[已修复,§九 R2:新增 3 个端口对角拒绝/边界用例]**
3. **带零点群延迟无独立锚点**。→ 对应 P3-2。**[已修复,§九 R5:新增 with-zeros FD 锚点测试]**
4. **`verify_against_spec` 只有 good/tightened 两个用例**(`verify/mod.rs:966,1009`),无"RL 不满足 spec 但综合成功"的端到端用例(因为管线不强制,这类用例需要显式构造)。
5. 基准 `bench-data/golden/twoport_anchors.json` 的来源是 **Cameron 2018 印刷矩阵的转录**(文件 `:4` 自述,原出 fuqtgui/tests/cameron_golden.py)—— 这是**外部独立数据**,`benchmark_standard.rs` 用它做 transversal-vs-folded 响应不变性锚点(1e-10 容差),**非循环验证**。这是测试独立性的强项,记录在案。

---

## 五、文档及发布状态冲突

### 5.1 基线差异(重要)

- 用户给定基线:HEAD 7dd617a、0.6.2-unreleased、435 通过、工作树干净。
- 实测:**HEAD 一致(7dd617a)**;但**工作树在审查期间从干净 0.6.2 变为 dirty 0.6.3**(12 文件未提交:`pole_expansion.rs` +206、`generalized_chebyshev_helpers.rs` +152、`admittance.rs`、`design.rs`、`synthesis_regression.rs` +61、CHANGELOG +64、版本号 0.6.2→0.6.3 等)。committed HEAD 的 Cargo.toml 仍是 `0.6.2`;工作树为 `0.6.3`。
- CHANGELOG:committed 顶部是 `## [Unreleased]`(0.6.2 的内容全在未发布段);工作树顶部新增 `## [0.6.2] — 2026-10-11` 已发布段,0.6.3 的改动进了新的 `[Unreleased]`。**"0.6.2 已于 10-11 发布"与工作树一致,与 committed 状态不一致**(committed 里 0.6.2 仍标 Unreleased,且**没有 v0.6.2 git tag**,只有 v0.6.0/v0.6.1)。
- 用户报告的 435 通过对应的是审查前的干净 0.6.2 树;**我在 0.6.3 工作树上实测 431 通过、0 失败**(round-2 新增测试 12 个:5 个 pole_expansion 契约 + 5 个 find_eps + 2 个 integration,净增与既有测试结构调整并存,总数 431)。两个数字都对,但属于不同的树 —— 报告以 0.6.3 工作树为准。

### 5.2 文档与代码的不一致

| # | 文档 | 实际 | 依据 |
|---|------|------|------|
| 1 | README.md:48 "stable up to order 30+" | order 30+ 被 typed 拒绝,实测上限 28 | P2-1 |
| 2 | architecture.md:573 "High order (10–30) ✅" | 同上 | P2-1 |
| 3 | `pole_expansion.rs:88-96` "Input contract (enforced)" | 端口对角部分未强制 | P2-2 **[已修复,§九 R2:端口对角纳入 `assert_transversal_structure` 并写入契约文档]** |
| 4 | `CHANGELOG.md` [Unreleased] 声称 round-2 修复"已完成" | 属实(我独立验证通过,§5.3),但**未提交** | §5.1 |

### 5.3 round-2 两项修复的独立验收(通过)

1. **fast path 结构前置检查**:探针 A —— 合法 transversal 矩阵被接受(fast=LU 逐点 ≤ 2e-14);注入 O(1) 谐振器间耦合 → `PreconditionViolation`(typed,信息指明 M[i][j] 与容差);folded/arrow 矩阵被拒;`response_normalized` 在折叠设计上正确回退 LU(与 LU 参考 1e-6 一致)。容差 `1e-9·max(1,max|M|)` 尺度感知正确:低于容差(1e-10/1e-12)的 round-off 耦合被接受,响应偏差 ≤ ~1e-6(物理可忽略);主流程合成矩阵实测端口对角恒为 0。
2. **all-finite eps≤1 拒绝**:`find_eps` 守卫(`generalized_chebyshev_helpers.rs:284-297`)同时覆盖单带(`:476`)与多带(`multi_band.rs:416`)两个入口;数学正确 —— order-2 全有限 F=0.25s²+1、P=4s²+1、RL20 时 eps=0.75/3/√99≈0.251<1,守卫触发,错误信息点名 eps 值、域与可行替代(mixed 零点集)。mixed 零点(eps_r=1)不受影响(既有测试 + 新测试 `find_eps_mixed_zeros_always_feasible_with_eps_r_one`)。可行域 eps>1 行为不变(N12 RL25 全有限设计仍合成且命中 spec)。

---

## 六、实际执行的检查及结果

| 检查 | 命令/方法 | 结果 |
|------|-----------|------|
| 全量测试(0.6.3 工作树) | `cargo test --all-targets` | **431 通过,0 失败**(lib 299 含新契约测试;集成测试全绿;~3.5 min) |
| Clippy 默认档 | `cargo clippy --all-targets` | 0 warning,exit 0 |
| Clippy cli+python 档 | `cargo clippy --all-targets --features cli,python` | 0 warning,exit 0 |
| RL spec 符合性探针 | `/mnt/c/Users/eynulai/dev/mfs_r3rlcheck`(1001 点,-1..1) | 7 设计全部 MEETS SPEC:min-RL 于带边 = spec(偏差 <0.0001 dB);\|S21\|max ≤ 1.0;N12 全有限 RL25 命中(此前用 1e-9 过严判据误报,改用库的 0.1 dB 门余量后通过) |
| 高阶稳定性扫描 | 同探针,RL20 | all-pole OK 至 N=28;N=30/32/40 被 typed 拒绝(极点实部 1.41e-4/3.19e-4/8.48e-3 > 容差);finite-tz N=30 被拒、N=32 OK |
| `high_order_stability` example | `cargo run --example high_order_stability` | order 30 两例打印 "Matrix synthesis failed"(吞错继续)—— 证实文档声明无 example 背书 |
| from_matrix 契约探针 A | 仓库外 crate:合法 transversal / 注入 M[1,2]=1e-9..1e-3 / 端口 M[0,0]=4 / folded / arrow | 见 §5.3;端口对角 4 被接受但响应错(S21 差 ~1.4e-1)→ P2-2 |
| 5 点门盲区 | 探针 D:fast 与 LU 的 5 点 \|S21\| 对比 | 源 shunt 不移动 5 点 \|S21\| 相对差(<1e-3 门)时通过,全网格差 ~1e-1 —— 门只比较 \|S21\| 不比较相位/S11(`design.rs:544-550`),但主流程已不可达(结构检查前置) |
| git 基线 | `git log/status/tag` | HEAD 7dd617a 一致;工作树 dirty(12 文件);tags 仅 v0.6.0/v0.6.1 |

---

## 七、未覆盖范围

- 三个审计代理(DocsAudit/TestIndependence/InterfacesAudit)的最终汇总未在本报告截止前交付;其逐条文档锚点核验、测试循环性全表、多工器/LM/cascade/Python-CLI 逐入口守卫审查以我的抽查替代,未做穷举。
- 多工器/双工器 `full_iteration(true)` 的非收敛路径、macOS/Windows CI 矩阵(仅本地 Linux 验证)。
- Python wheel 构建与 CLI 二进制行为(未构建)。
- 复零点群延迟(§P3-2)未实测,仅代码审查。
- 性能基准未跑(CI 只编译不执行,本轮保持一致)。

---

## 八、最优先处理的三项

1. **修正 "stable to 30+" 声明**(P2-1):README.md:48 与 architecture.md:573 改为实测边界(order ≤ 28 保证;30+ 可能 typed 拒绝),或让 `high_order_stability` example 断言成功并把可行上限固化为回归锚。依据:发布文档当前做不实能力声明,且该 example 会吞掉综合失败。**[已修复,§九 R1]**
2. **补全 `from_matrix` 的端口对角契约**(P2-2):`assert_transversal_structure` 增加 `M[0,0]`/`M[N+1,N+1]` 检查(或文档明确"端口 shunt 请用 LU"),并加 1 个拒绝用例。依据:实测存在合法 `CouplingMatrix` 输入被静默给出 ~1e-1 量级错误响应,且现有 5 个契约测试 + 5 点门都拦不住。**[已修复,§九 R2]**
3. **提交并固化 0.6.3 状态**(§5.1):round-2 修复已验证通过,但版本号 bump、CHANGELOG 重排、12 个文件改动全部未提交,且没有 v0.6.2 tag。当前仓库处于"committed 0.6.2-unreleased / 工作树 0.6.3-unreleased"的分裂状态,发布追溯性受损;建议按 CHANGELOG 现状提交(或打 v0.6.2 tag 后提交 0.6.3 段)。**[待用户决定:是否提交/打 tag 属发布动作,本修复未代执行]**

---

## 九、修复记录(2026-10-13,round-3 缺陷全部关闭)

上述 P2-1 / P2-2 / P3-1 / P3-2 已在本工作树全部修复并通过完整 CI 命令验证(本地,WSL2 / Rust 1.95.0)。逐项如下:

**R1(P2-1,高阶边界)** — README 与 `architecture.md` 的 "stable to 30+" 改为实测边界表述(all-pole ≤ order 28 保证;finite-TZ ≤ order 24 @ RL 20–22 dB;order 30+ 可能 typed `NumericalFailure` 拒绝)。`examples/high_order_stability.rs` 拆为"必须成功的已验证集"(任一回归即 exit 非零)与"显式标注的探索扫描"(order 30+,拒绝为预期);新增 `#[test]` 锁定边界(order 28 成功、order 30 拒绝)。

**R2(P2-2,端口对角契约)** — `assert_transversal_structure`(`src/response/pole_expansion.rs`)现同时检查 `M[0,0]` 与 `M[N+1,N+1]`:超过尺度感知容差 `1e-9·max(1,max|M|)` 的端口对角以 `PreconditionViolation` 拒绝(与谐振器块同一容差),`from_matrix` 契约文档写明"纯串联端接,端口 shunt 请用 LU 后端"。新增 3 个 `pole_expansion` 测试:`source_port_diagonal_is_rejected`、`load_port_diagonal_is_rejected`、`port_diagonal_tolerance_boundary`(1e-12 接受 / 1e-8 拒绝)。既有两个使用 `M[3][3]=0.05` 占位值的 fixture 已改为 0.0(纯串联),使拒绝归因唯一、夹具落在快路径合法域内。主流程合成矩阵端口对角恒为 0,故无主路径回归。

**R3(P2-2 连带,5 点门)** — 快路径 5 点交叉验证门从"仅比较 |S21|"升级为共享 `fast_path_gate` 助手:比较**复数** S21(幅值 + 相位)与 S11,任一非有限值即自动回退 LU。两处调用点(`response_normalized` / `response_bandpass`)改用该助手。新增 5 个 `design` 测试,含回归用例(等幅、相位旋转——旧仅幅值门会放行)。

**R4(P3-1,文档)** — `synthesize` 文档明确"综合成功不保证 RL 满足 spec,需显式调用 `verify_against_spec`";`verify_against_spec` 新增失败处理 doctest(17 个 doctest 全绿)。

**R5(P3-2,群延迟锚点)** — 新增 `response` 测试 `group_delay_fd_anchor_on_with_zeros_design`:对带有限传输零点(±2, ±2.5)的 6 阶设计,在相位平滑区(|S21| > 0.2)用 S21 相位中心差分做**独立**值锚点(LU 与快路径双双对照),附 step 收敛守护(n=401 vs n=801,差 < 1e-3)与零点区检查(传输零点处 `group_delay` 必须为 `None`,两个后端均查)。

**R6(§八.3,发布状态)— 未执行,待用户决定**:提交/打 tag 属发布动作。round-2 + round-3 全部改动仍在工作树未提交;按 CHANGELOG 现状提交或打 v0.6.2 tag 后提交 0.6.3 段均可,需用户定夺。

### 验证矩阵(修复后,全部实测)

| 检查 | 命令 | 结果 |
|------|------|------|
| 默认档(库 + 集成 + doctest) | `cargo test` | **457 通过,0 失败**(= lib 308 + 集成 132 + doctest 17) |
| cli,python 档(库 + 集成 + doctest) | `cargo test --features cli,python` | **473 通过,0 失败**(= 默认档 + cli 集成 16) |
| 全目标档(库 + 集成 + bench + example 单测,**不含 doctest**) | `cargo test --all-targets` | **441 通过,0 失败**(= lib 308 + 集成 132 + example `#[test]` 1) |
| Clippy 默认档 | `cargo clippy --all-targets` | 0 warning |
| Clippy cli,python 档 | `cargo clippy --features cli,python --all-targets` | 0 warning |

**计数口径(各命令是不同 target 集合,不可跨集合相加):**
- 默认档 `cargo test` = lib 308 + 集成 132 + doctest 17 = **457**(不编译 example/bench 作为测试目标)。
- cli,python 档 = 默认档 + cli 集成 16(`cli_integration` / `cli_save_resume_e2e` 等)= **473**。
- 全目标档 `cargo test --all-targets` = lib 308 + 集成 132 + example `#[test]` 1 = **441**;它把 `examples/high_order_stability.rs` 内的 `#[test]` 边界回归测试编入运行,但**不跑 doctest**(doctest 由 `cargo test` 或 `cargo test --doc` 触发)。
- **441 + 17 = 458 ≠ 457 的差异原因**:默认档 457 含 17 个 doctest、不含那 1 个 example `#[test]`;全目标档 441 含那 1 个 example `#[test]`、不含 17 个 doctest。即 `441 = 457 + 1(example 测试) − 17(doctest)`。两命令跑的是互不相同的 target 集合,把 441 与 17 直接相加是错误口径,多出的 1 正是 `examples/high_order_stability.rs` 的边界回归测试。

---

## 附:验证边界声明

本轮验证覆盖了核心数学链路(Cameron 递推→residue→transversal→变换→响应→带通映射)的自洽性与 round-2 修复的正确性,并有外部数据锚点(Cameron 2018 印刷矩阵)与独立 FD 锚点(all-pole 群延迟)支撑;本轮发现的 round-3 缺陷(§二)已在本报告 §九 中修复并复核通过。**本报告不宣称整个库已证明正确**;修复范围外的未覆盖范围见 §七,待验证风险见 §三。所有缺陷均经两次复核(调用上下文、既有保护、数学约定、测试覆盖),未采用"测试数量"或"Clippy 干净"作为正确性论据。
