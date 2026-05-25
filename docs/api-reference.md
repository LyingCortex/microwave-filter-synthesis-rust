# MFS API 参考手册

## 概览

```rust
use mfs::prelude::*;
```

一行导入，获得所有需要的类型。

---

## 1. 创建滤波器设计

### 1.1 带通滤波器（最常用）

```rust
let design = FilterDesign::bandpass(order, return_loss_db, center_hz, bandwidth_hz)
    .zeros_hz([...])      // 可选：物理频率零点
    .unloaded_q(value)    // 可选：无载Q值
    .synthesize()?;
```

**参数说明：**

| 参数 | 类型 | 说明 | 示例 |
|------|------|------|------|
| `order` | usize | 谐振器数量（滤波器阶数） | `6` |
| `return_loss_db` | f64 | 通带回波损耗 (dB) | `23.0` |
| `center_hz` | f64 | 中心频率 (Hz) | `6.75e9` |
| `bandwidth_hz` | f64 | 3dB 带宽 (Hz) | `300e6` |

**完整示例：**

```rust
// 6阶，23dB回损，中心6.75GHz，带宽300MHz，3个传输零点，Q=3000
let design = FilterDesign::bandpass(6, 23.0, 6.75e9, 300e6)
    .zeros_hz([6.4e9, 6.5e9, 7.0e9])
    .unloaded_q(3000.0)
    .synthesize()?;
```

### 1.2 归一化原型（学术/调试用）

```rust
let design = FilterDesign::prototype(order, return_loss_db)
    .zeros([...])         // 可选：归一化零点坐标
    .unloaded_q(value)    // 可选
    .synthesize()?;
```

**示例：**

```rust
// 4阶全极点，20dB回损
let design = FilterDesign::prototype(4, 20.0).synthesize()?;

// 4阶带2个对称零点
let design = FilterDesign::prototype(4, 20.0)
    .zeros([-1.5, 1.5])
    .synthesize()?;
```

---

## 2. 获取耦合矩阵

```rust
// 横向矩阵（transversal，综合直接输出）
let m = design.matrix();

// 折叠矩阵（folded）
let m = design.to_folded()?;

// 箭头矩阵（arrow）
let m = design.to_arrow()?;
```

**矩阵操作：**

```rust
let m = design.to_folded()?;

m.order()              // 阶数 → usize
m.side()               // 矩阵维度 (order+2) → usize
m.at(row, col)         // 读取元素 → Option<f64>
m.source_coupling()    // 源耦合 |M[0,1]|
m.load_coupling()      // 负载耦合 |M[N,N+1]|
m.chain_couplings()    // 相邻耦合列表 → Vec<f64>
m.as_slice()           // 底层行主序数据 → &[f64]
```

---

## 3. 计算 S 参数频率响应

### 3.1 带通设计（自动使用存储的中心频率/带宽）

```rust
let design = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;

// 只需给起止频率和点数
let response = design.response(6.5e9, 7.0e9, 201)?;
```

### 3.2 指定带通参数

```rust
let response = design.response_bandpass(
    center_hz,      // 6.75e9
    bandwidth_hz,   // 300e6
    start_hz,       // 6.5e9
    stop_hz,        // 7.0e9
    points,         // 201
)?;
```

### 3.3 归一化频率响应

```rust
let response = design.response_normalized(
    start,    // -3.0（归一化频率）
    stop,     //  3.0
    points,   //  201
)?;
```

### 3.4 对变换后的矩阵求响应

```rust
let folded = design.to_folded()?;
let response = design.eval(&folded, -3.0, 3.0, 201)?;
```

---

## 4. 读取 S 参数数据

```rust
let response = design.response(6.5e9, 7.0e9, 201)?;

for sample in &response.samples {
    // 频率
    sample.frequency_hz       // 物理频率 (Hz)
    sample.normalized_omega   // 归一化频率

    // S21（传输）
    sample.s21_mag()          // 幅度（线性）
    sample.s21_db()           // 幅度 (dB)
    sample.s21_phase_deg()    // 相位（度）
    sample.s21_re             // 实部
    sample.s21_im             // 虚部

    // S11（反射）
    sample.s11_mag()          // 幅度（线性）
    sample.s11_db()           // 幅度 (dB)
    sample.s11_phase_deg()    // 相位（度）
    sample.s11_re             // 实部
    sample.s11_im             // 虚部

    // 群延迟
    sample.group_delay        // 群延迟（归一化单位）
}
```

---

## 5. 带通缩放

```rust
// 带通设计：自动使用存储的参数
let design = FilterDesign::bandpass(4, 20.0, 6.75e9, 300e6).synthesize()?;
let scaled = design.scale()?;

// 手动指定参数
let scaled = design.scale_to(6.75e9, 300e6)?;

// 缩放任意矩阵
let folded = design.to_folded()?;
let scaled_folded = design.scale_matrix(&folded, 6.75e9, 300e6)?;
```

---

## 6. 访问内部数据

```rust
design.order()          // 阶数
design.center_hz()      // 中心频率 → Option<f64>
design.bandwidth_hz()   // 带宽 → Option<f64>
design.spec()           // 完整规格 → &FilterSpec
design.polynomials()    // 原型多项式 (E, F, P) → &PolynomialSet
design.matrix()         // 横向耦合矩阵 → &CouplingMatrix
```

---

## 7. 完整工作流示例

### 7.1 典型带通滤波器设计

```rust
use mfs::prelude::*;

fn main() -> Result<()> {
    // 设计
    let design = FilterDesign::bandpass(6, 23.0, 6.75e9, 300e6)
        .zeros_hz([6.4e9, 6.5e9, 7.0e9])
        .unloaded_q(3000.0)
        .synthesize()?;

    // 矩阵
    let folded = design.to_folded()?;
    println!("源耦合: {:.4}", folded.source_coupling());
    println!("负载耦合: {:.4}", folded.load_coupling());
    println!("链耦合: {:?}", folded.chain_couplings());

    // 响应
    let response = design.response(6.0e9, 7.5e9, 501)?;
    for s in &response.samples {
        println!("{:.4} GHz | S21={:.2} dB | S11={:.2} dB",
            s.frequency_hz / 1e9, s.s21_db(), s.s11_db());
    }

    Ok(())
}
```

### 7.2 对比不同拓扑的响应

```rust
use mfs::prelude::*;

fn main() -> Result<()> {
    let design = FilterDesign::prototype(4, 20.0)
        .zeros([-1.5, 1.5])
        .synthesize()?;

    let transversal = design.response_normalized(-3.0, 3.0, 101)?;
    let folded = design.to_folded()?;
    let folded_resp = design.eval(&folded, -3.0, 3.0, 101)?;

    // 两者应该完全一致（相似变换不改变响应）
    for (a, b) in transversal.samples.iter().zip(folded_resp.samples.iter()) {
        assert!((a.s21_db() - b.s21_db()).abs() < 1e-6);
    }

    Ok(())
}
```

---

## 8. 错误处理

所有可能失败的操作返回 `Result<T>`，错误类型为 `MfsError`：

```rust
match FilterDesign::bandpass(0, 20.0, 6.75e9, 300e6).synthesize() {
    Ok(design) => { /* ... */ }
    Err(MfsError::InvalidOrder { order }) => {
        eprintln!("阶数无效: {order}");
    }
    Err(e) => {
        eprintln!("其他错误: {e}");
    }
}
```

常见错误：
- `InvalidOrder` — 阶数为 0
- `InvalidReturnLoss` — 回损 ≤ 0 或非有限值
- `InvalidFrequency` — 频率参数无效
- `NumericalFailure` — 数值计算未收敛（极高阶时可能出现）
