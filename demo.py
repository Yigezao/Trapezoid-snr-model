"""可运行的演示：扫频看 SNR 形状 → 找低于阈值的未覆盖频段 → 小网格寻参。

用法：  python demo.py

会打印：
  1. 三组示例驱动的 SNR 扫频对比、低于 3dB 的未覆盖频段；
  2. 以 [55, 88, 120] kHz 为目标的极小规模 maximin 寻参结果。

若本机装有 matplotlib，还会把对比曲线存成 docs/demo_snr_curves.png。
所有数字均为演示用假设值，不代表任何具体产品参数。
"""

from __future__ import annotations

import numpy as np

from snr_model import snr_for_freqs, sweep_snr, uncovered_spans

# 三组示例配置（纯演示值）
CONFIGS = [
    {"name": "cfg-A", "ta": 6.4, "na": 16, "nb": 1, "vh": 4.0, "vl": 0.5, "k": 1.8},
    {"name": "cfg-B", "ta": 8.0, "na": 8, "nb": 1, "vh": 4.0, "vl": 0.5, "k": 1.8},
    {"name": "cfg-C", "ta": 5.0, "na": 12, "nb": 2, "vh": 3.6, "vl": 0.4, "k": 2.0},
]

THRESHOLD_DB = 3.0      # 演示用灵敏度阈值
TARGETS_KHZ = [55.0, 88.0, 120.0]  # 目标干扰频点


def _params(cfg):
    return {k: v for k, v in cfg.items() if k != "name"}


def report_coverage():
    print("=" * 64)
    print("1) 扫频 SNR 对比（10–500 kHz, 步长 1 kHz, 阈值 %g dB）" % THRESHOLD_DB)
    print("=" * 64)
    for cfg in CONFIGS:
        freq, snr = sweep_snr(10.0, 1.0, 500.0, **_params(cfg), a_sin=0.1)
        spans = uncovered_spans(freq, snr, THRESHOLD_DB)
        n_uncovered = int(round(sum(f2 - f1 + 1 for f1, f2 in spans)))
        line = f"{cfg['name']:<6s} ta={cfg['ta']:<4} na={cfg['na']:<3} "
        line += f"| 未覆盖频段数 {len(spans):<2d} | 累计宽度 {n_uncovered:>3d} kHz"
        print(line)
        for f1, f2 in spans[:6]:
            print(f"           ├─ {f1:7.1f} ~ {f2:7.1f} kHz")


def run_optimize():
    print("\n" + "=" * 64)
    print("2) 极小规模 maximin 寻参（目标 %s kHz）" % TARGETS_KHZ)
    print("=" * 64)
    targets = np.array(TARGETS_KHZ, dtype=float)
    best = None
    for ta in (5.0, 6.4, 8.0):      # 演示只枚举 3x3 = 9 组
        for na in (8, 12, 16):
            snr = snr_for_freqs(targets, ta=ta, na=na, nb=1,
                                vh=4.0, vl=0.5, k=1.8, a_sin=0.1)
            worst = float(np.nanmin(snr))
            if best is None or worst > best[0]:
                best = (worst, ta, na, snr)
    worst, ta, na, snr = best
    print("最优: ta=%g us, na=%d → 目标处 SNR(dB) = %s"
          % (ta, na, np.round(snr, 2).tolist()))
    print("其中最弱频点 %g dB（maximin 目标值）。" % worst)


def plot_curves(save_path: str = "docs/demo_snr_curves.png"):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
        plt.rcParams["axes.unicode_minus"] = False
    except ImportError:
        print("\n[提示] 未安装 matplotlib，跳过绘图。")
        return

    fig, ax = plt.subplots(figsize=(8.4, 4.8), dpi=130)
    for cfg in CONFIGS:
        freq, snr = sweep_snr(10.0, 1.0, 500.0, **_params(cfg), a_sin=0.1)
        ax.plot(freq, snr, lw=1.3, label=cfg["name"])
    ax.axhline(THRESHOLD_DB, color="tab:red", ls="--", lw=1,
               label="灵敏度阈值 %g dB" % THRESHOLD_DB)
    for t in TARGETS_KHZ:
        ax.axvline(t, color="tab:gray", ls=":", lw=0.8, alpha=0.7)
    ax.set_xlabel("频率 [kHz]")
    ax.set_ylabel("SNR [dB]")
    ax.set_title("示例：梯形驱动 + 正弦扰动的频点 SNR（演示值，非真实数据）")
    ax.set_ylim(-12, 40)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc="upper right")
    fig.tight_layout()
    fig.savefig(save_path)
    print(f"\n对比曲线已保存: {save_path}")
    plt.close(fig)


if __name__ == "__main__":
    report_coverage()
    run_optimize()
    plot_curves()
