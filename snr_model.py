"""通用方法框架：梯形脉冲驱动 + 正弦扰动 → 各频点 SNR 解析模型。

本文件不含任何厂商芯片参数或实测数据，仅实现一套可复现的建模方法：

问题背景
--------
不少模拟 / 传感前端的激励是"分段梯形脉冲"：
每个激励段由若干脉冲组成，脉冲在高低电平间以有限斜率翻转，
相邻脉冲极性交替以形成交流激励；相邻激励段之间留一段静默时间。
当存在频率为 f 的外部正弦扰动时，扰动按确定性的采样关系进入解调结果，
使某些频点上的信噪比明显塌陷。我们希望解析地预测 SNR(f)，
从而在给定目标干扰频点时反推驱动参数，或评估某组参数的覆盖能力。

模型要点
--------
1. 单个半周期(宽度 ta)内的梯形积分是闭式可求的：上升沿 + 平台两段；
2. 叠加正弦扰动后，每个脉冲的积分多出一项解析可积的正弦项；
3. 解调按"段平均"进行：把若干个相邻段划成一个 Block 取平均，
   信号取所有 Block 的均值、噪声取 Block 间标准差，两者之比即 SNR。

所有电平 / 时间 / 幅值 / 段数都是**输入参数**。gain 是把模型输出
换算到具体量纲的标定系数——未标定时取 1，此时 dB 曲线只反映相对形状。
"""

from __future__ import annotations

import numpy as np

# 频率范围与分段常数不再内置：全部由调用方传入。
# 以下仅保留推导出的物理量命名约定。

# ---------------------------------------------------------------------------
# 核心：给定一组驱动参数，求指定频点上的 SNR
# ---------------------------------------------------------------------------


def snr_for_freqs(
    freqs_khz: np.ndarray,
    *,
    vh: float = 4.0,       # 驱动高电平 [V]
    vl: float = 0.5,       # 驱动低电平 [V]
    k: float = 1.8,        # 上升/下降斜率 [V/us]
    ta: float = 6.4,       # 单脉冲(半周期)宽度 [us]
    td: float = 1.5,       # 激励段间静默(死区) [us]
    na: int = 16,          # 每激励段脉冲数
    nb: int = 1,           # 每 Block 的段数（段平均窗口）
    n_blocks: int = 200,   # Block 总数
    a_sin: float = 0.1,    # 外部正弦扰动幅值 [V]
    gain: float = 1.0,     # 量纲/链路标定增益（未标定时为 1）
) -> np.ndarray:
    """返回与 freqs_khz 等长的 SNR_dB 数组。

    params ：
        freqs_khz : 待求频点 [kHz]，任意形状，自动展平处理。
    其余均为关键字参数，见上面的默认值说明。
    """
    tr = (vh - vl) / k
    if tr >= ta:
        raise ValueError(
            f"上升时间 tr={tr:.3f}us >= 半周期 ta={ta}us，波形不存在平台段"
        )

    # 梯形半周期面积（伏·微秒）
    phi_trap = (vl + vh) / 2.0 * tr + vh * (ta - tr)

    # 脉冲序号 j = 0..na-1，极性交替（偶正奇负）
    j = np.arange(na)
    signs = np.where(j % 2 == 0, 1.0, -1.0)  # (na,)

    # 段骨架：第 i 段的第 j 个脉冲起始时刻
    t_seg = na * ta + td
    n_seg = nb * n_blocks
    i_seg = np.arange(n_seg)
    t0 = (i_seg * t_seg)[:, None] + j[None, :]  # (n_seg, na)，与频率无关

    # ---- 频率轴向量化：扰动正弦的解析积分 ----
    freq = np.atleast_1d(np.asarray(freqs_khz, dtype=np.float64)).ravel()
    om = 2.0 * np.pi * freq * 1e-3  # 角频率 [rad/us]
    om = om[:, None, None]          # (n_freq, 1, 1)

    t0_b = t0[None, :, :]           # (1, n_seg, na)
    # 正弦扰动在半周期 [t0, t0+ta] 上的积分
    sin_int = (a_sin / om) * (
        np.cos(om * t0_b) - np.cos(om * (t0_b + ta))
    )                               # (n_freq, n_seg, na)

    # 每个脉冲总积分 = 梯形面积 + 交替符号的正弦项
    phi = phi_trap + signs[None, None, :] * sin_int
    seg_sum = phi.sum(axis=2)       # (n_freq, n_seg)：各激励段积分

    # Block 平均：每 nb 个相邻段取一次平均
    blocks = seg_sum.reshape(-1, n_blocks, nb).mean(axis=2)  # (n_freq, n_blocks)

    mu = blocks.mean(axis=1)        # 信号：Block 均值
    sd = blocks.std(axis=1)         # 噪声：Block 间标准差

    # SNR(线性) = mu * gain / (sd * f)；mu<=0 或 sd==0 处无定义
    snr_lin = np.full_like(mu, np.nan)
    ok = (sd > 0) & (mu > 0)
    snr_lin[ok] = mu[ok] * gain / (sd[ok] * freq[ok])
    snr_db = np.where(ok, 10.0 * np.log10(snr_lin), np.nan)

    return snr_db.reshape(np.shape(freqs_khz))


# ---------------------------------------------------------------------------
# 便捷：整段扫频
# ---------------------------------------------------------------------------


def sweep_snr(
    f_start: float = 10.0,
    f_step: float = 1.0,
    f_end: float = 500.0,
    **kwargs,
):
    """在 [f_start, f_end] kHz 上以 f_step 步长扫频，返回 (freq_khz, snr_db)。"""
    freq = np.arange(f_start, f_end + 0.5 * f_step, f_step)
    return freq, snr_for_freqs(freq, **kwargs)


# ---------------------------------------------------------------------------
# 覆盖统计：低于阈值的连续频段
# ---------------------------------------------------------------------------


def uncovered_spans(freq_khz, snr_db, threshold_db: float):
    """找出 SNR 连续低于阈值的频段，返回 [(f_start, f_end), ...]。"""
    if len(snr_db) == 0:
        return []
    below = np.where(np.isfinite(snr_db) & (snr_db < threshold_db), 1, 0)
    spans = []
    start = None
    for idx, flag in enumerate(below):
        if flag and start is None:
            start = idx
        elif not flag and start is not None:
            spans.append((freq_khz[start], freq_khz[idx - 1]))
            start = None
    if start is not None:
        spans.append((freq_khz[start], freq_khz[-1]))
    return spans
