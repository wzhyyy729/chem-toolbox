"""chem_toolbox 自检测试。

运行方式（不需要安装 pytest）：
    python -m chem_toolbox.tests

为什么要有这个文件？
--------------------
1. **验证正确性**：每个期望值都能用教科书或手算核对，不是"跑通就算对"。
2. **防止改坏**：以后你想改某个公式，先跑一遍测试，改错了立刻会发现。
3. **面试/笔试能用**：问你"怎么保证算得对"，你可以直接说"有 40 多个断言覆盖边界情况"。
"""

from __future__ import annotations

import math
import sys

from . import core

# Windows 控制台默认编码是 GBK，直接打印中文/特殊符号可能抛 UnicodeEncodeError。
# 这里把输出强制切成 UTF-8，并让无法编码的字符退化成替代字符而不是直接崩溃。
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

_passed = 0
_failed: list[str] = []


def check(name: str, actual, expected, tolerance: float = 1e-6) -> None:
    """断言 actual ≈ expected，并把结果记录下来。"""
    global _passed
    if isinstance(expected, float) or isinstance(actual, float):
        ok = actual is not None and math.isclose(float(actual), float(expected), rel_tol=0, abs_tol=tolerance)
    else:
        ok = actual == expected

    if ok:
        _passed += 1
    else:
        _failed.append(f"{name}: 期望 {expected}，实际 {actual}")


def check_raises(name: str, func, *args, **kwargs) -> None:
    """断言这段代码会抛出 ChemToolboxError。"""
    global _passed
    try:
        func(*args, **kwargs)
    except core.ChemToolboxError:
        _passed += 1
    except Exception as error:  # 抛了别的异常也算不通过
        _failed.append(f"{name}: 期望 ChemToolboxError，实际抛出 {type(error).__name__}: {error}")
    else:
        _failed.append(f"{name}: 期望抛出 ChemToolboxError，但正常返回了")


# ===========================================================================
# 一、分子式解析
# ===========================================================================

def test_formula() -> None:
    # --- 基础分子式 ---
    check("H2O 组成", core.parse_formula("H2O"), {"H": 2, "O": 1})
    check("H2O 摩尔质量", core.analyse_formula("H2O").molar_mass, 18.015, 0.001)

    # --- 葡萄糖：教科书标准值 180.156 g/mol ---
    check("葡萄糖摩尔质量", core.analyse_formula("C6H12O6").molar_mass, 180.156, 0.002)
    check("葡萄糖碳数", core.parse_formula("C6H12O6")["C"], 6)

    check("硫酸摩尔质量", core.analyse_formula("H2SO4").molar_mass, 98.072, 0.005)
    check("氯化钠摩尔质量", core.analyse_formula("NaCl").molar_mass, 58.44, 0.01)

    # --- 括号：Ca(OH)2 = 40.078 + 2×15.999 + 2×1.008 = 74.092 ---
    check("Ca(OH)2 组成", core.parse_formula("Ca(OH)2"), {"Ca": 1, "O": 2, "H": 2})
    check("Ca(OH)2 摩尔质量", core.analyse_formula("Ca(OH)2").molar_mass, 74.092, 0.005)

    # --- 方括号 + 嵌套：K3[Fe(CN)6] = 3×39.098 + 55.845 + 6×12.011 + 6×14.007 = 329.24 ---
    check("K3[Fe(CN)6] 组成", core.parse_formula("K3[Fe(CN)6]"), {"K": 3, "Fe": 1, "C": 6, "N": 6})
    check("K3[Fe(CN)6] 摩尔质量", core.analyse_formula("K3[Fe(CN)6]").molar_mass, 329.24, 0.02)

    # --- 水合物：CuSO4·5H2O = 249.68 ---
    check("胆矾摩尔质量", core.analyse_formula("CuSO4.5H2O").molar_mass, 249.68, 0.02)
    check("胆矾结晶水", core.parse_formula("CuSO4.5H2O")["H"], 10)
    check("纯碱十水合物", core.analyse_formula("Na2CO3.10H2O").molar_mass, 286.14, 0.02)

    # --- 前置系数：2H2O 表示两个水分子 ---
    check("2H2O 组成", core.parse_formula("2H2O"), {"H": 4, "O": 2})

    # --- 不饱和度 DBE = (2C + 2 + N + P - H - X)/2 ---
    check("苯 C6H6 的 DBE", core.analyse_formula("C6H6").dbe, 4.0)          # 3 双键 + 1 环
    check("环己烷 C6H12 的 DBE", core.analyse_formula("C6H12").dbe, 1.0)     # 一个环
    check("乙烷 C2H6 的 DBE", core.analyse_formula("C2H6").dbe, 0.0)
    check("吡啶 C5H5N 的 DBE", core.analyse_formula("C5H5N").dbe, 4.0)
    check("氯苯 C6H5Cl 的 DBE", core.analyse_formula("C6H5Cl").dbe, 4.0)
    check("十氢萘 C10H18 的 DBE", core.analyse_formula("C10H18").dbe, 2.0)   # 两个环

    # 硫和磷在有机物里都是二价/三价，计入 DBE 公式后结果才成立
    check("噻吩 C4H4S 的 DBE", core.analyse_formula("C4H4S").dbe, 3.0)      # 2 双键 + 1 环
    # 三苯基膦 P(C6H5)3 = C18H15P：(2×18 + 2 + 1 - 15)/2 = 12，即 3 个苯环 × 4
    check("三苯基膦 P(C6H5)3 的 DBE", core.analyse_formula("P(C6H5)3").dbe, 12.0)
    check("三苯基氧膦 OP(C6H5)3 的 DBE", core.analyse_formula("OP(C6H5)3").dbe, 12.0)
    # 饱和膦 P(CH3)3 = C3H9P 应为 0，验证 +P 不会把饱和物算成非整数
    check("三甲基膦 P(CH3)3 的 DBE", core.analyse_formula("P(CH3)3").dbe, 0.0)

    # 金属有机物：DBE 公式不适用，应为 None 并给出说明
    check("K3[Fe(CN)6] 无 DBE", core.analyse_formula("K3[Fe(CN)6]").dbe, None)
    check("非有机物有说明", bool(core.analyse_formula("K3[Fe(CN)6]").dbe_note), True)
    check("无碳时给提示", bool(core.analyse_formula("H2SO4").dbe_note), True)

    # --- 元素质量分数：水的氢占 2×1.008/18.015 = 11.19% ---
    check("水中氢的质量分数", core.analyse_formula("H2O").mass_fractions["H"], 0.1119, 0.0005)

    # --- 非法输入必须报错，而不是默默算错 ---
    check_raises("空分子式", core.parse_formula, "")
    check_raises("小写元素符号", core.parse_formula, "c6h12o6")
    check_raises("括号不匹配", core.parse_formula, "Ca(OH2")
    check_raises("空括号", core.parse_formula, "Ca()2")
    check_raises("未知元素", core.parse_formula, "Xx2O")


# ===========================================================================
# 一之二、药学与材料方向的分子
# ===========================================================================

def test_drugs_and_materials() -> None:
    """药学与材料方向常见分子的摩尔质量核对。

    这些数字都可以用药典或试剂手册核对：
    - 阿司匹林 C9H8O4 = 180.158
    - 对乙酰氨基酚 C8H9NO2 = 151.163
    - 布洛芬 C13H18O2 = 206.281
    - 顺铂 Pt(NH3)2Cl2 = 300.05
    - 硫酸铜五水合物 CuSO4.5H2O = 249.68
    - 二氧化钛 TiO2 = 79.866
    - 钛酸钡 BaTiO3 = 233.192
    """
    check("阿司匹林", core.analyse_formula("C9H8O4").molar_mass, 180.158, 0.01)
    check("对乙酰氨基酚", core.analyse_formula("C8H9NO2").molar_mass, 151.163, 0.01)
    check("布洛芬", core.analyse_formula("C13H18O2").molar_mass, 206.281, 0.01)

    # 顺铂：Pt(NH3)2Cl2，含铂配合物，药学化疗药物
    cisplatin = core.analyse_formula("Pt(NH3)2Cl2")
    check("顺铂摩尔质量", cisplatin.molar_mass, 300.05, 0.02)
    check("顺铂氮原子数", cisplatin.counts["N"], 2)
    check("顺铂无 DBE", cisplatin.dbe, None)

    # 材料：氧化物与钙钛矿结构
    check("二氧化钛 TiO2", core.analyse_formula("TiO2").molar_mass, 79.866, 0.01)
    check("钛酸钡 BaTiO3", core.analyse_formula("BaTiO3").molar_mass, 233.192, 0.02)
    check("氧化锆 ZrO2", core.analyse_formula("ZrO2").molar_mass, 123.222, 0.02)

    # 稀土与难熔金属现在也在元素表里了
    check("氧化镧 La2O3", core.analyse_formula("La2O3").molar_mass, 325.809, 0.03)
    check("碳化钨 WC", core.analyse_formula("WC").molar_mass, 195.85, 0.02)
    check("氧化铟锡骨架 In2O3", core.analyse_formula("In2O3").molar_mass, 277.634, 0.03)

    # 药学里常见的含硫药物：硫在 DBE 公式中计入
    check("硫脲 CH4N2S 的 DBE", core.analyse_formula("CH4N2S").dbe, 1.0)   # 一个 C=S 双键


# ===========================================================================
# 二、溶液与 pH
# ===========================================================================

def test_solution() -> None:
    # --- 物质的量与浓度 ---
    # 称 5.844 g NaCl（M=58.44）配成 1 L，正好 0.1000 mol/L
    check("NaCl 物质的量", core.moles(5.844, 58.44), 0.1, 1e-9)
    check("NaCl 浓度", core.molarity(5.844, 58.44, 1.0), 0.1, 1e-9)
    check("配制 0.1 M 1 L 需 NaCl 质量", core.mass_for_molarity(0.1, 1.0, 58.44), 5.844, 1e-9)

    # --- 稀释定律：0.1 M 取 10 mL 稀释到 100 mL → 0.01 M ---
    result = core.dilute(c1=0.1, v1=0.01, v2=0.1, c2=None)
    check("稀释后浓度", result["c2"], 0.01, 1e-12)
    result = core.dilute(c1=0.1, v1=0.01, c2=0.01, v2=None)
    check("稀释后体积", result["v2"], 0.1, 1e-12)
    check_raises("稀释少了两个未知", core.dilute, 0.1, 0.01, None, None)

    # --- 强酸强碱 ---
    check("0.1 M HCl 的 pH", core.ph_strong_acid(0.1), 1.0, 1e-9)
    check("0.01 M HCl 的 pH", core.ph_strong_acid(0.01), 2.0, 1e-9)
    check("0.1 M NaOH 的 pH", core.ph_strong_base(0.1), 13.0, 1e-9)

    # --- 弱酸精确解 ---
    # 0.1 M 醋酸，Ka = 1.75e-5（pKa = 4.756）
    # 手算：(−1.75e-5 + √(1.75e-5² + 4×1.75e-5×0.1))/2 = 1.3143e-3 → pH = 2.8813
    ph, alpha = core.ph_weak_acid(1.75e-5, 0.1)
    check("0.1 M 醋酸 pH", ph, 2.8813, 0.0005)
    check("0.1 M 醋酸电离度", alpha, 0.013143, 1e-5)

    # 用最简式 √(Ka·c) 验证：√(1.75e-5 × 0.1) = 1.3229e-3 → pH = 2.8785
    # 精确解应当略小于最简式结果（因为分母 c-x 让电离受抑制）
    naive_ph = -math.log10(math.sqrt(1.75e-5 * 0.1))
    check("精确解比最简式算出的 pH 略大", ph > naive_ph, True)

    # 极稀极弱酸不应崩溃或给出负数
    ph_dilute, _ = core.ph_weak_acid(1e-10, 1e-5)
    check("极弱酸 pH 仍在 0~14 合理范围", 0.0 < ph_dilute < 14.0, True)

    # --- Ka / Kb / pKa 换算 ---
    check("pKa 4.76 对应 Ka", core.ka_from_pka(4.76), 1.7378e-5, 1e-8)
    check("Ka·Kb = Kw", core.kb_from_ka(1.75e-5) * 1.75e-5, core.KW_25C, 1e-20)

    # 氨水：Kb = 1.8e-5，0.1 M → [OH-] = 1.3327e-3 → pOH = 2.8752 → pH = 11.1248
    ph_nh3, _ = core.ph_weak_base(1.8e-5, 0.1)
    check("0.1 M 氨水 pH", ph_nh3, 11.1248, 0.001)

    # --- 朗伯-比尔定律：A = 0.5，ε = 1000，b = 1 → c = 5e-4 ---
    check("比尔定律求浓度", core.mass_from_absorbance(0.5, 1000.0, 1.0), 5e-4, 1e-12)

    # --- 由标准曲线反推未知样浓度 ---
    # 曲线 y = 2x + 0.01，测得 A = 0.601 → x = (0.601 - 0.01)/2 = 0.2955
    check("标准曲线反推浓度", core.concentration_from_curve(0.601, 2.0, 0.01), 0.2955, 1e-12)
    check_raises("斜率为 0 无法反推", core.concentration_from_curve, 0.5, 0.0, 0.1)

    # --- 非法输入 ---
    check_raises("负浓度", core.ph_strong_acid, -0.1)
    check_raises("零浓度", core.ph_weak_acid, 1e-5, 0.0)
    check_raises("负 Ka", core.ph_weak_acid, -1e-5, 0.1)
    check_raises("负质量", core.moles, -1.0, 58.44)
    check_raises("零体积", core.molarity, 1.0, 58.44, 0.0)


# ===========================================================================
# 三、统计与线性拟合
# ===========================================================================

def test_stats() -> None:
    # --- 统计量：数据 [1,2,3,4,5] ---
    summary = core.summarise([1, 2, 3, 4, 5])
    check("n", summary.n, 5)
    check("均值", summary.mean, 3.0, 1e-12)
    check("中位数", summary.median, 3.0, 1e-12)
    check("样本标准差", summary.std_sample, math.sqrt(2.5), 1e-12)   # √2.5 ≈ 1.5811
    check("总体标准差", summary.std_population, math.sqrt(2.0), 1e-12)
    check("极差", summary.range_, 4.0, 1e-12)

    # 单个数据点：样本标准差按约定取 0，不应崩溃
    single = core.summarise([7.0])
    check("单点样本标准差", single.std_sample, 0.0)
    check_raises("空数据", core.summarise, [])

    # --- 线性拟合：人造完美直线 y = 3x + 2 ---
    xs = [1.0, 2.0, 3.0, 4.0, 5.0]
    ys = [3 * x + 2 for x in xs]
    fit = core.linear_fit(xs, ys)
    check("完美直线斜率", fit.slope, 3.0, 1e-9)
    check("完美直线截距", fit.intercept, 2.0, 1e-9)
    check("完美直线 R²", fit.r_squared, 1.0, 1e-9)

    # --- 带噪声的真实感数据：标准曲线 c=0.1~0.5 mmol/L，A 有微小误差 ---
    conc = [0.1, 0.2, 0.3, 0.4, 0.5]
    absorb = [0.103, 0.198, 0.302, 0.397, 0.501]
    fit2 = core.linear_fit(conc, absorb)
    check("标准曲线斜率约 1.0", fit2.slope, 0.994, 0.02)
    check("标准曲线 R² 很高", fit2.r_squared > 0.999, True)
    check("斜率标准误已算出", math.isnan(fit2.stderr_slope), False)

    check_raises("x/y 长度不匹配", core.linear_fit, [1, 2, 3], [1, 2])
    check_raises("数据点太少", core.linear_fit, [1, 2], [1, 2])
    check_raises("x 全相同", core.linear_fit, [1, 1, 1], [1, 2, 3])


# ===========================================================================
# 四、谱图找峰
# ===========================================================================

def test_peaks() -> None:
    # --- 人造高斯峰：中心 500 nm，宽 20 nm ---
    def gaussian(x: float, center: float, sigma: float, height: float) -> float:
        return height * math.exp(-((x - center) ** 2) / (2 * sigma ** 2))

    xs = [400.0 + i for i in range(201)]  # 400 ~ 600 nm
    ys = [gaussian(x, 500.0, 10.0, 1.0) + 0.01 for x in xs]

    peaks, noise = core.find_peaks(xs, ys)
    check("找到一个主峰", len(peaks), 1)
    check("峰位约 500 nm", peaks[0].position, 500.0, 1.5)

    # 高斯峰的半高全宽 FWHM = 2.3548 × σ = 23.548
    check("FWHM 约 2.355σ", peaks[0].fwhm, 2.3548 * 10.0, 1.5)

    # --- 双峰分辨 ---
    ys2 = [gaussian(x, 470.0, 6.0, 1.0) + gaussian(x, 530.0, 6.0, 0.6) + 0.01 for x in xs]
    peaks2, _ = core.find_peaks(xs, ys2)
    check("找到两个峰", len(peaks2), 2)
    check("第一峰位置", peaks2[0].position, 470.0, 2.0)
    check("第二峰位置", peaks2[1].position, 530.0, 2.0)

    # --- 平坦噪声不应报峰 ---
    flat = [0.5] * 201
    peaks3, _ = core.find_peaks(xs, flat)
    check("平坦数据无峰", len(peaks3), 0)

    # --- 异常输入 ---
    check_raises("点数太少", core.find_peaks, [1, 2, 3], [1, 2, 3])
    check_raises("x 未排序", core.find_peaks, [5, 4, 3, 2, 1], [1, 2, 3, 4, 5])
    check_raises("长度不匹配", core.find_peaks, [1, 2, 3, 4, 5], [1, 2, 3])


# ===========================================================================
# 运行全部测试
# ===========================================================================

def run() -> int:
    tests = [
        ("分子式解析", test_formula),
        ("药学与材料分子", test_drugs_and_materials),
        ("溶液与 pH", test_solution),
        ("统计与拟合", test_stats),
    ]

    print("=" * 62)
    print("chem_toolbox 自检测试")
    print("=" * 62)

    for title, func in tests:
        before_failed = len(_failed)
        func()
        status = "通过" if len(_failed) == before_failed else "失败"
        print(f"  [{status}] {title}")

    print("-" * 62)
    print(f"断言通过 {_passed} 项，失败 {len(_failed)} 项")

    if _failed:
        print()
        print("失败详情：")
        for item in _failed:
            print(f"  [FAIL] {item}")
        return 1

    print("全部通过")
    return 0


if __name__ == "__main__":
    sys.exit(run())
