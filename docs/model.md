# 模型推导：梯形驱动 + 正弦扰动的频点 SNR

> 这是一份**方法框架**文档。所有量都用输入参数表达，不绑定任何厂商芯片
> 参数或实测数据；把模型换算到真实量纲需要一个**标定增益 G**（见 §4），
> 未标定时取 G = 1，dB 曲线只反映**相对形状**。

## 1. 驱动波形

激励由"激励段"重复构成，一个激励段包含 `na` 个脉冲：

- 每个脉冲宽度为一个半周期 `ta`；脉冲内是**梯形波**：从低电平 `vl` 以
  有限斜率 `k` 上升到高电平 `vh`，剩余时间为平台；
- 相邻脉冲**极性交替**（相当于交流激励，用于抵消直流偏置）；
- 每 `na` 个脉冲之后留一段静默时间 `td`，再开始下一激励段。

单段总时长

$$\boxed{T_{seg} = na \cdot ta + td}$$

上升时间与平台时间：

$$\boxed{tr = \frac{vh - vl}{k}, \qquad t_{plateau} = ta - tr \ \ (>0)}$$

## 2. 梯形积分（无扰动时的"信号"）

单个半周期内梯形波的面积（伏·微秒）：

$$\boxed{\Phi_{trap} = \frac{vl + vh}{2}\cdot tr + vh \cdot t_{plateau}}$$

## 3. 叠加正弦扰动

假设存在幅值 $A_{sin}$、角频率 $\omega$ 的外部扰动 $A_{sin}\sin(\omega t)$，
它叠加在每一个脉冲的采样窗口上。第 $i$ 段第 $j$ 个脉冲的起始时刻为

$$\boxed{t_{0}(i,j) = i\cdot T_{seg} + j\cdot ta}$$

扰动在窗口 $[t_{0},\, t_{0}+ta]$ 上的积分为闭式解：

$$\int_{t_{0}}^{t_{0}+ta} A_{sin}\sin(\omega t)\,dt
= \frac{A_{sin}}{\omega}\Big[\cos(\omega t_{0}) - \cos\big(\omega(t_{0}+ta)\big)\Big]$$

考虑脉冲极性交替（$j$ 偶取正、$j$ 奇取负），该脉冲的总积分：

$$\Phi_{j}(f)
= \Phi_{trap} + (-1)^{j}\cdot \frac{A_{sin}}{\omega}
\Big[\cos(\omega t_{0}) - \cos\big(\omega(t_{0}+ta)\big)\Big]$$

## 4. Block 平均与 SNR

将解调后的信号按 **Block** 组织：每 `nb` 个相邻激励段划入一个 Block，
`n_blocks` 个 Block 覆盖全部 $N_{seg} = nb \cdot n_{blocks}$ 段。

- 段求和 $S_{i}(f) = \sum_{j=0}^{na-1}\Phi_{j}(f)$
- Block 均值 $B_{k}(f) = \tfrac{1}{nb}\sum_{m} S_{k\cdot nb + m}(f)$

于是对每个频点 $f$：

$$\mu(f) = \frac{1}{n_{blocks}}\sum_k B_k(f),
\qquad \sigma(f) = \operatorname{std}_k\big(B_k(f)\big)$$

把 Block 间波动视为等效噪声、Block 均值视为信号，定义频点 SNR：

$$\boxed{SNR_{dB}(f) = 10\log_{10}\frac{\mu(f)\cdot G}{\sigma(f)\cdot f}}$$

说明：

- $\mu(f)\le 0$ 或 $\sigma(f)=0$ 处 SNR 无定义（置为 NaN）；
- $f$ 归一化保证不同频点可比（噪声与积分窗口数相关）；
- **$G$ 是量纲 / 接收链路标定增益**：它把"伏·微秒量级的数值"换算到
  真实信噪比。$G$ 需要针对具体硬件标定，本框架不假设其数值。

## 5. 关键数值示例（演示值）

取 $ta=6.4\,\mu s,\ na=16,\ nb=1,\ vh=4.0\,\text{V},\ vl=0.5\,\text{V},\ k=1.8\,\text{V/}\mu s$：

- $tr = (4.0-0.5)/1.8 = 1.94\,\mu s$
- $T_{seg} = 16\times 6.4 + td$
- 对 $f$ 从 10 kHz 扫描到 500 kHz，`snr_model.sweep_snr` 给出完整曲线，
  可见在某些频点（周期结构与扰动频率满足"整数/半整数锁定"关系处）SNR 显著塌陷。

这正是反向应用的价值：**给定需要屏蔽的目标干扰频点，反过来枚举驱动参数，
使目标处 $\min_f SNR$ 最大或覆盖低于阈值的频段最少**（见 `demo.py` 的
maximin 小演示）。

## 6. 实现与向量化

`snr_model.py` 的核心是一次三维广播（避免按频率写 Python 循环）：

```
频率轴 ω [n_freq, 1, 1]
  × 脉冲时刻 t0 [1, n_seg, na]
  → 正弦积分 sin_int [n_freq, n_seg, na]
  → 叠加梯形面积 + 极性 → 段求和 [n_freq, n_seg]
  → Block 均值 [n_freq, n_blocks] → μ / σ → SNR_dB [n_freq]
```

计算量为 $O(n_{freq}\times n_{seg}\times na)$，全部向量化，参数寻优时对每组
候选只重算目标频点即可（开销极小）。
