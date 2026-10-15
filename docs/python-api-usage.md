# MFS Python API 使用手册

## 安装

```bash
cd /path/to/mfs
maturin develop --features python --release
# 或构建wheel:
maturin build --features python --release -o dist/
pip install dist/mfs-*.whl
```

## 快速开始

```python
import mfs

# 设计一个6阶带通滤波器
d = mfs.bandpass(order=6, rl=23.0, center=2e9, bw=100e6,
                 zeros=[2090e6, 1910e6, 2150e6])

# 获取耦合矩阵
print(d.folded())   # 折叠形式
print(d.arrow())    # 箭头形式

# 计算S参数
freqs, s21_db, s11_db = d.response(1.8e9, 2.2e9, 201)
```

---

## 滤波器创建

### bandpass() — 从中心频率+带宽

```python
d = mfs.bandpass(
    order=6,           # 阶数
    rl=23.0,           # 回波损耗 [dB]
    center=2e9,        # 中心频率 [Hz] (几何均值)
    bw=100e6,          # 带宽 [Hz]
    zeros=[2.09e9, 1.91e9],  # 传输零点 [Hz] (可选)
    q=1000.0,          # 无载Q值 (可选, 默认无损)
)
```

> **注意**: `center`应为passband edges的**几何均值** `√(f_low × f_high)`。
> 使用`bandpass_edges()`可自动计算。

### bandpass_edges() — 从通带边缘频率

```python
d = mfs.bandpass_edges(
    order=6,
    rl=23.0,
    f_start=1.95e9,    # 通带下边缘 [Hz]
    f_stop=2.05e9,     # 通带上边缘 [Hz]
    zeros=[2.09e9, 1.91e9],
    q=1000.0,
)
```

### prototype() — 归一化原型

```python
d = mfs.prototype(
    order=4,
    rl=20.0,
    zeros=[1.8, -1.8],  # 归一化Omega值
    q=500.0,
)
```

---

## 耦合矩阵获取

所有方法返回2D Python列表 (N+2)×(N+2)。

### 基本拓扑

```python
d.matrix()          # 横向矩阵 (Transversal)
d.folded()          # 折叠矩阵 (Folded) — Cameron标准形式
d.arrow()           # 箭头矩阵 (Arrow)
d.scaled_matrix()   # 带通缩放矩阵 (物理频率单位)
```

### CT — 级联三元组 (Cascaded Trisection)

每个零点对应一个trisection交叉耦合。

```python
# z_loc: [from, to] 表示交叉耦合位置, span=2
# center = (from+to)//2

# 6阶3零点CT
d = mfs.bandpass(6, 25.0, 2e9, 100e6, zeros=[2080e6, 1920e6, 2120e6])
cm_ct = d.ct([[1, 3], [2, 4], [4, 6]])
# 结果: M[1,3]≠0, M[2,4]≠0, M[4,6]≠0
```

**z_loc规则**:
- `[k, k+2]`: trisection跨resonator k到k+2, 中心=k+1
- center必须 ≥2 且 ≤N-1 (不旋转Source/Load直连的谐振器)

### CQ — 级联四元组 (Cascaded Quadruplet)

每对零点对应一个quadruplet对角交叉耦合。

```python
# z_loc混合: span=2(trisection部分) + span=3(quadruplet对角)
# 程序自动配对和推导

# 6阶2零点CQ, M[1,4]对角耦合, 保留M[1,3]
d = mfs.bandpass(6, 25.0, 2e9, 100e6, zeros=[2080e6, 1920e6])
cm_cq = d.cq([[1, 3], [1, 4]])

# 8阶4零点, 两个CQ
d = mfs.bandpass(8, 25.0, 2e9, 100e6, zeros=[2080e6, 1920e6, 2150e6, 1850e6])
cm_cq = d.cq([[1, 3], [1, 4], [5, 7], [5, 8]])
```

**z_loc规则**:
- `[k, k+2]` (span=2): trisection交叉耦合 → 保留的那条
- `[k, k+3]` (span=3): quadruplet对角 → 定义4谐振器块范围
- 每个quadruplet需要一个span=2和一个span=3配对
- 程序自动推导`common_res`: span=2起始于quadruplet start→1, 否则→4

**common_res效果** (非对称零点):
| z_loc | common_res | 保留 | 消除 |
|-------|-----------|------|------|
| `[[1,3],[1,4]]` | 1 | M13 | M24=0 |
| `[[2,4],[1,4]]` | 4 | M24 | M13=0 |

### BOX — 方阵拓扑

每个零点对应一个4谐振器BOX块。

```python
# spans: [start, end] with end-start=3

# 4阶1零点BOX
d = mfs.bandpass(4, 25.0, 2e9, 100e6, zeros=[2090e6])
cm_box = d.box_topology([[1, 4]])

# 8阶2零点, 两个BOX
d = mfs.bandpass(8, 25.0, 2e9, 100e6, zeros=[2090e6, 1910e6])
cm_box = d.box_topology([[1, 4], [5, 8]])
```

**BOX结构**: M(k,k+2)≠0, M(k+1,k+3)≠0, M(k+1,k+2)=0

**与CQ区别**:
- BOX: 1个零点/块, 消除inline耦合(Jacobi旋转)
- CQ: 2个零点/块, 保留inline耦合(cross-rotation合并)

---

## S参数计算

### response() — 幅度响应 (dB)

```python
freqs, s21_db, s11_db = d.response(
    start=1.8e9,    # 起始频率 [Hz]
    stop=2.2e9,     # 终止频率 [Hz]
    points=201,     # 频率点数
)
```

### response_lossy() — 有损响应

```python
freqs, s21_db, s11_db = d.response_lossy(
    start=1.8e9,
    stop=2.2e9,
    points=201,
    q=1000.0,       # 无载Q值
)
```

### response_full() — 完整响应 (幅度+相位+群延迟)

```python
data = d.response_full(start=1.8e9, stop=2.2e9, points=201)
# 返回 list of dict:
# [{'freq': Hz, 's21_db': dB, 's11_db': dB,
#   's21_mag': linear, 's11_mag': linear,
#   's21_phase': deg, 's11_phase': deg,
#   'group_delay': float or None}, ...]
# group_delay 单位是归一化原型单位 dτ/dω(不是秒);LU 与极点展开快路径都会计算,
# 仅在数值退化处(传输零点 S21→0)为 None —— 显式标记,不是假 0.0
```

---

## Touchstone导出

```python
# 导出为字符串
s2p_text = d.to_touchstone(freq_unit="GHz", format="RI", impedance=50.0)

# 直接保存文件
d.save_touchstone("filter.s2p", freq_unit="MHz", format="DB")
```

---

## 耦合矩阵优化/调谐

```python
# 从目标S参数调谐耦合矩阵
result = d.tune(
    freqs=[...],              # 频率点 [Hz]
    target_s21_re=[...],      # 目标S21实部
    target_s21_im=[...],      # 目标S21虚部
    target_s11_re=[...],      # 目标S11实部
    target_s11_im=[...],      # 目标S11虚部
    max_iter=200,             # 最大迭代次数
)
# result: {'matrix': [[...]], 'cost': float, 'iterations': int, 'converged': bool}
```

---

## Duplexer设计

```python
dpx = mfs.duplexer(
    tx_order=4, tx_rl=23.0, tx_center=1.96e9, tx_bw=60e6,
    rx_order=4, rx_rl=23.0, rx_center=2.14e9, rx_bw=60e6,
    tx_zeros=[1.92e9],
    rx_zeros=[2.18e9],
)

resp = dpx.response(1.8e9, 2.3e9, 401)
# resp: dict with 'frequencies', 'tx_s21_db', 'rx_s21_db',
#       'tx_s11_db', 'rx_s11_db', 'antenna_s11_db', 'isolation_db'

print(dpx.tx_matrix())  # Tx滤波器耦合矩阵
print(dpx.rx_matrix())  # Rx滤波器耦合矩阵
```

---

## 设计属性

```python
d.order         # 滤波器阶数
d.center        # 中心频率 [Hz] (None for prototype)
d.bw            # 带宽 [Hz] (None for prototype)
d.f_start       # 通带下边缘 [Hz]
d.f_stop        # 通带上边缘 [Hz]
```

---

## 频率约定

| API参数 | 单位 | 说明 |
|---------|------|------|
| `center`, `bw` | Hz | 几何均值中心频率 |
| `f_start`, `f_stop` | Hz | 通带边缘 |
| `zeros` (bandpass) | Hz | 物理频率零点 |
| `zeros` (prototype) | 无量纲 | 归一化Omega |
| `response` freqs | Hz | 输入输出都是Hz |

**频率变换**: `Ω = (1/FBW) × (f/f0 - f0/f)`
- `f0 = center` (几何均值)
- `FBW = bw / center`
- `Ω(f_start) = -1`, `Ω(f_stop) = +1` (精确)
