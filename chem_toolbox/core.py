"""化学计算核心模块（chem_toolbox.core）

设计原则
--------
1. **只用 Python 标准库**：不需要 pip 安装任何东西，拷到任何装好 Python 的电脑上都能跑。
2. **每个函数只做一件事**：方便单独测试，也方便你以后自己改。
3. **中间步骤一起返回**：不只是给答案，还把过程给你看，方便对照手算结果。

命名约定
--------
- 摩尔质量单位一律 g/mol
- 浓度单位一律 mol/L
- 温度默认 25 ℃（此时水的离子积 Kw = 1.0e-14）

作者：wzhyyy729 · 化学专业大二在读 · 代码由 AI 协作产出并逐行理解
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

# 25 ℃ 时水的离子积常数
KW_25C = 1.0e-14

# 相对原子质量（IUPAC 常用值，保留到小数点后 3 位）
# 覆盖范围：基础无机 + 有机 + 药学（Pt/Sn/Sb/Te/As/Se）+ 材料（Zr/Hf/Nb/Ta/W/Mo/稀土/In/Ga/Ge）
ATOMIC_MASS: Dict[str, float] = {
    # --- 第一、二周期 ---
    "H": 1.008, "He": 4.003,
    "Li": 6.94, "Be": 9.012, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998, "Ne": 20.180,
    # --- 第三周期 ---
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948,
    # --- 第四周期 ---
    "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942,
    "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693,
    "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.630, "As": 74.922,
    "Se": 78.971, "Br": 79.904, "Kr": 83.798,
    # --- 第五周期 ---
    "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906,
    "Mo": 95.95, "Ru": 101.07, "Rh": 102.906, "Pd": 106.42, "Ag": 107.868,
    "Cd": 112.414, "In": 114.818, "Sn": 118.710, "Sb": 121.760, "Te": 127.60,
    "I": 126.904, "Xe": 131.293,
    # --- 第六周期 ---
    "Cs": 132.905, "Ba": 137.327, "La": 138.905, "Ce": 140.116, "Pr": 140.908,
    "Nd": 144.242, "Sm": 150.36, "Eu": 151.964, "Gd": 157.25, "Tb": 158.925,
    "Dy": 162.500, "Ho": 164.930, "Er": 167.259, "Tm": 168.934, "Yb": 173.045,
    "Lu": 174.967, "Hf": 178.486, "Ta": 180.948, "W": 183.84, "Re": 186.207,
    "Os": 190.23, "Ir": 192.217, "Pt": 195.084, "Au": 196.967, "Hg": 200.592,
    "Tl": 204.38, "Pb": 207.2, "Bi": 208.980,
    # --- 第七周期（天然存在）---
    "Th": 232.038, "U": 238.029,
}


class ChemToolboxError(ValueError):
    """本工具包统一使用的异常类型。

    遇到不合法的输入（比如写错的分子式、负浓度）时抛出，
    这样调用方只需要捕获这一种异常就够了。
    """


# ===========================================================================
# 一、分子式分析
# ===========================================================================

# 匹配"元素符号 + 可选下标"，例如 C、H2、Na、Cl2
_ELEMENT_RE = re.compile(r"([A-Z][a-z]?)(\d*)")

# 匹配"占位符 + 可选下标"，例如 @0@、@1@6
_TOKEN_RE = re.compile(r"@(\d+)@(\d*)")

# 匹配"元素 或 占位符 + 可选下标"，用于解析已展开的式子
_PART_RE = re.compile(r"@\d+@\d*|[A-Z][a-z]?\d*")

# 匹配最内层括号组（括号内不允许再有括号），以及紧跟其后的下标
_INNERMOST_GROUP_RE = re.compile(r"\(([^()]*)\)(\d*)")


def _normalise_brackets(text: str) -> str:
    """把方括号、花括号统一换成圆括号，方便统一处理。

    化学式里 [] 和 () 都起"把一个原子团括起来"的作用，例如
    K3[Fe(CN)6]、Cu(NH3)4SO4。统一后只需要处理一种括号，逻辑简单很多。
    遇到 ( ] 这种不配对的情况会直接报错，避免悄悄算错。
    """
    pairs = {"(": ")", "[": "]", "{": "}"}
    closers = {v: k for k, v in pairs.items()}
    stack: List[str] = []
    out: List[str] = []

    for character in text:
        if character in pairs:
            stack.append(character)
            out.append("(")
        elif character in closers:
            if not stack or stack.pop() != closers[character]:
                raise ChemToolboxError(
                    f"化学式 '{text}' 的括号不配对（{character} 没有对应的左括号）"
                )
            out.append(")")
        else:
            out.append(character)

    if stack:
        raise ChemToolboxError(f"化学式 '{text}' 的括号不配对（有左括号没有闭合）")
    return "".join(out)


def _resolve_tokens(text: str, tokens: Dict[str, str]) -> Dict[str, int]:
    """把式子里的占位符递归展开，返回 {元素: 个数}。

    为什么需要递归？
    ----------------
    嵌套括号是"从里往外"替换的，所以外层占位符的内容里可能还嵌着内层占位符。
    以 K3[Fe(CN)6] 为例，处理过程是：
        K3(Fe(CN)6)
          → 内层 (CN) 换成 @0@     得到 K3(Fe@0@6)
          → 外层 (Fe@0@6) 换成 @1@ 得到 K3@1@
    展开 @1@ 时拿到的是 "Fe@0@6"，里面还有 @0@，所以必须继续递归展开。
    """
    resolved: Dict[str, int] = {}
    position = 0

    for match in _PART_RE.finditer(text):
        if match.start() != position:
            raise ChemToolboxError(
                f"式子片段 '{text}' 在第 {position + 1} 个字符处无法解析，"
                f"请检查大小写（元素符号首字母大写，如 Na、Cl）"
            )
        position = match.end()
        part = match.group(0)

        token_match = _TOKEN_RE.fullmatch(part)
        if token_match:  # 占位符：取回内容继续展开，再乘上它自己的下标
            inner_text = tokens[f"@{token_match.group(1)}@"]
            counts = _resolve_tokens(inner_text, tokens)
            factor = int(token_match.group(2)) if token_match.group(2) else 1
            counts = {element: count * factor for element, count in counts.items()}
        else:  # 普通元素
            element_match = _ELEMENT_RE.fullmatch(part)
            if element_match is None:  # 理论上不会发生，留作保险
                raise ChemToolboxError(f"式子片段 '{text}' 中的 '{part}' 无法解析")
            element = element_match.group(1)
            if element not in ATOMIC_MASS:
                raise ChemToolboxError(
                    f"元素 '{element}' 不在相对原子质量表中。"
                    f"如需支持，请在 chem_toolbox/core.py 的 ATOMIC_MASS 里补充"
                )
            counts = {element: int(element_match.group(2)) if element_match.group(2) else 1}

        for element, count in counts.items():
            resolved[element] = resolved.get(element, 0) + count

    if position != len(text):
        raise ChemToolboxError(
            f"式子片段 '{text}' 在第 {position + 1} 个字符处无法解析，"
            f"请检查大小写（元素符号首字母大写，如 Na、Cl）"
        )
    return resolved


def parse_formula(formula: str) -> Dict[str, int]:
    """解析化学式，返回 {元素: 原子个数}。

    支持
    ----
    - 常见分子式       : "H2O"、"C6H12O6"、"H2SO4"
    - 括号与嵌套括号   : "Ca(OH)2"、"K3[Fe(CN)6]"、"Cu(NH3)4SO4"
    - 水合物（点号）   : "CuSO4.5H2O"、"Na2CO3.10H2O"
    - 前置系数         : "2H2O"（会把整段乘以 2）

    实现思路（从最内层括号往外剥，比写递归下降解析器短得多）
    -------------------------------------------------------
    第一遍：反复用正则找出**不含括号**的最小括号组 `(...)`，把它整体换成占位符
            `@n@`，并记录它的内容文本，直到式子里没有括号。
    第二遍：统一把占位符递归展开成元素计数（下标在展开时相乘）。
    这样天然支持任意层嵌套，例如 K3[Fe(CN)6] 这种方括号套圆括号的写法。
    """
    if not formula or not formula.strip():
        raise ChemToolboxError("分子式不能为空")

    text = _normalise_brackets(formula.strip().replace(" ", ""))

    # --- 处理水合物：按点号切开，逐段解析后相加 ---
    if "." in text:
        total: Dict[str, int] = {}
        for index, part in enumerate(text.split(".")):
            if not part:
                raise ChemToolboxError(f"分子式 '{formula}' 中存在空的片段（点号位置不对）")
            if index == 0:
                sub = parse_formula(part)
            else:
                # "." 后面的部分支持前置系数，例如 CuSO4.5H2O 中的 5H2O
                leading_match = re.fullmatch(r"(\d+)(.+)", part)
                if leading_match:
                    factor = int(leading_match.group(1))
                    sub = parse_formula(leading_match.group(2))
                    sub = {element: count * factor for element, count in sub.items()}
                else:
                    sub = parse_formula(part)
            for element, count in sub.items():
                total[element] = total.get(element, 0) + count
        return total

    # --- 处理整式前置系数，例如 "2H2O" ---
    leading_match = re.fullmatch(r"(\d+)(.+)", text)
    if leading_match:
        factor = int(leading_match.group(1))
        inner = parse_formula(leading_match.group(2))
        return {element: count * factor for element, count in inner.items()}

    # --- 第一遍：把每一层括号替换成占位符，记录"占位符 → 内容文本" ---
    tokens: Dict[str, str] = {}
    counter = 0

    while "(" in text or ")" in text:
        match = _INNERMOST_GROUP_RE.search(text)
        if match is None:
            raise ChemToolboxError(f"分子式 '{formula}' 的括号不匹配")

        inner_text = match.group(1)
        if not inner_text:
            raise ChemToolboxError(f"分子式 '{formula}' 中存在空括号")

        suffix = match.group(2)
        token = f"@{counter}@"
        counter += 1
        tokens[token] = inner_text
        # 下标暂时原样保留，等展开阶段统一相乘
        text = text[: match.start()] + token + suffix + text[match.end():]

    # --- 第二遍：展开所有占位符，得到元素计数 ---
    return _resolve_tokens(text, tokens)


@dataclass
class FormulaReport:
    """分子式分析结果。"""

    formula: str
    counts: Dict[str, int]
    molar_mass: float
    mass_fractions: Dict[str, float] = field(default_factory=dict)
    dbe: float | None = None
    dbe_note: str = ""


def analyse_formula(formula: str) -> FormulaReport:
    """分析分子式：摩尔质量、元素质量分数、不饱和度（DBE）。

    不饱和度公式：DBE = (2C + 2 + N + P - H - X) / 2
      C = 碳原子数, N = 氮原子数, P = 磷原子数, H = 氢原子数, X = 卤素原子数

    为什么硫和磷要算进去？—— 它们在有机物里都表现为二价：
    噻吩 C4H4S 的 DBE 应为 3（两个双键 + 一个环），用这个公式能得到 3；
    若把硫排除在外，公式就不成立了。所以适用元素集合是
    C、H、N、O、P、S、卤素；含金属或其他元素的式子不适用。
    """
    counts = parse_formula(formula)

    molar_mass = 0.0
    for element, count in counts.items():
        molar_mass += ATOMIC_MASS[element] * count
    molar_mass = round(molar_mass, 4)

    mass_fractions = {
        element: round(ATOMIC_MASS[element] * count / molar_mass, 6)
        for element, count in counts.items()
    }

    # 判断 DBE 公式是否适用：只看 C/H/N/O/P/S/卤素组成的有机物
    allowed = {"C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"}
    dbe: float | None = None
    note = ""
    unsupported = sorted(set(counts) - allowed)

    if unsupported:
        note = (
            f"分子含 {'、'.join(unsupported)}，不饱和度公式仅适用于 "
            f"C/H/N/O/P/S/卤素组成的有机物，此处不适用"
        )
    elif counts.get("C", 0) == 0:
        note = "分子中没有碳，不饱和度公式（针对有机物）不适用"
    else:
        carbon = counts.get("C", 0)
        hydrogen = counts.get("H", 0)
        nitrogen = counts.get("N", 0)
        phosphorus = counts.get("P", 0)
        halogen = sum(counts.get(x, 0) for x in ("F", "Cl", "Br", "I"))
        dbe = (2 * carbon + 2 + nitrogen + phosphorus - hydrogen - halogen) / 2

    return FormulaReport(
        formula=formula,
        counts=counts,
        molar_mass=molar_mass,
        mass_fractions=mass_fractions,
        dbe=dbe,
        dbe_note=note,
    )


# ===========================================================================
# 二、溶液配制、稀释与 pH
# ===========================================================================

def moles(mass_g: float, molar_mass: float) -> float:
    """由质量(g)和摩尔质量(g/mol)求物质的量(mol)。"""
    if molar_mass <= 0:
        raise ChemToolboxError("摩尔质量必须大于 0")
    if mass_g < 0:
        raise ChemToolboxError("质量不能为负数")
    return mass_g / molar_mass


def molarity(mass_g: float, molar_mass: float, volume_l: float) -> float:
    """由质量、摩尔质量、溶液体积求物质的量浓度(mol/L)。"""
    if volume_l <= 0:
        raise ChemToolboxError("体积必须大于 0（单位：L）")
    return moles(mass_g, molar_mass) / volume_l


def mass_for_molarity(target_molarity: float, volume_l: float, molar_mass: float) -> float:
    """要配制指定浓度和体积的溶液，需要称取多少克溶质。"""
    if volume_l <= 0:
        raise ChemToolboxError("体积必须大于 0（单位：L）")
    if target_molarity < 0:
        raise ChemToolboxError("浓度不能为负数")
    return target_molarity * volume_l * molar_mass


def dilute(
    c1: float | None,
    v1: float | None,
    c2: float | None = None,
    v2: float | None = None,
) -> Dict[str, float]:
    """稀释定律 c1·V1 = c2·V2，给出其中三个量求第四个。

    用法：恰好留一个参数为 None。
    """
    given = sum(x is None for x in (c1, v1, c2, v2))
    if given != 1:
        raise ChemToolboxError("c1、V1、c2、V2 中必须且只能留一个未知（传 None）")
    for name, value in (("c1", c1), ("c2", c2)):
        if value is not None and value < 0:
            raise ChemToolboxError(f"{name} 不能为负数")
    for name, value in (("V1", v1), ("V2", v2)):
        if value is not None and value <= 0:
            raise ChemToolboxError(f"{name} 必须大于 0")

    if c1 is None:
        return {"c1": c2 * v2 / v1}  # type: ignore[operator]
    if v1 is None:
        return {"v1": c2 * v2 / c1}  # type: ignore[operator]
    if c2 is None:
        return {"c2": c1 * v1 / v2}  # type: ignore[operator]
    return {"v2": c1 * v1 / c2}  # type: ignore[operator]


def ph_strong_acid(concentration: float) -> float:
    """强酸溶液 pH（忽略水的电离，适用于 c ≥ 1e-6 mol/L）。"""
    if concentration <= 0:
        raise ChemToolboxError("浓度必须大于 0")
    return -math.log10(concentration)


def ph_strong_base(concentration: float, kw: float = KW_25C) -> float:
    """强碱溶液 pH（忽略水的电离）。"""
    if concentration <= 0:
        raise ChemToolboxError("浓度必须大于 0")
    poh = -math.log10(concentration)
    return -math.log10(kw) - poh


def ph_weak_acid(ka: float, concentration: float) -> Tuple[float, float]:
    """一元弱酸溶液 pH，返回 (pH, 电离度 α)。

    推导（物料守恒 + 电荷守恒，只考虑一元弱酸 HA）
    ------------------------------------------------
    设电离出的 [H+] = x，则 [A-] = x，[HA] = c - x
        Ka = x² / (c - x)
    =>  x² + Ka·x - Ka·c = 0
    =>  x = (-Ka + √(Ka² + 4·Ka·c)) / 2      （取正根）

    注意：当 c·Ka 很小时（极弱酸），x 会趋近 √(Ka·c)，这正是最简式的来源；
    但本函数用的是**精确解**，所以极弱酸和大浓度下都比最简式更准。
    """
    if ka <= 0:
        raise ChemToolboxError("Ka 必须大于 0")
    if concentration <= 0:
        raise ChemToolboxError("浓度必须大于 0")

    x = (-ka + math.sqrt(ka * ka + 4.0 * ka * concentration)) / 2.0
    return -math.log10(x), x / concentration


def ph_weak_base(kb: float, concentration: float, kw: float = KW_25C) -> Tuple[float, float]:
    """一元弱碱溶液 pH，返回 (pH, 电离度 α)。

    与弱酸完全对称：设 [OH-] = x，则 Kb = x² / (c - x)。
    """
    if kb <= 0:
        raise ChemToolboxError("Kb 必须大于 0")
    if concentration <= 0:
        raise ChemToolboxError("浓度必须大于 0")

    x = (-kb + math.sqrt(kb * kb + 4.0 * kb * concentration)) / 2.0
    poh = -math.log10(x)
    return -math.log10(kw) - poh, x / concentration


def ka_from_pka(pka: float) -> float:
    """由 pKa 求 Ka。"""
    return 10.0 ** (-pka)


def kb_from_ka(ka: float, kw: float = KW_25C) -> float:
    """共轭酸碱对：Ka · Kb = Kw。"""
    if ka <= 0:
        raise ChemToolboxError("Ka 必须大于 0")
    return kw / ka


def mass_from_absorbance(
    absorbance: float,
    epsilon: float,
    path_length_cm: float,
) -> float:
    """朗伯-比尔定律 A = ε·b·c，已知 A、ε、b 求浓度 c(mol/L)。

    ε 的单位是 L·mol⁻¹·cm⁻¹，b 的单位是 cm。
    """
    if epsilon <= 0:
        raise ChemToolboxError("摩尔吸光系数 ε 必须大于 0")
    if path_length_cm <= 0:
        raise ChemToolboxError("光程 b 必须大于 0")
    if absorbance < 0:
        raise ChemToolboxError("吸光度不能为负数")
    return absorbance / (epsilon * path_length_cm)


# ===========================================================================
# 三、实验数据处理：统计量与线性拟合
# ===========================================================================

@dataclass
class StatsSummary:
    """一组数据的统计摘要。"""

    n: int
    mean: float
    std_sample: float       # 样本标准差（除以 n-1，实验报告一般用这个）
    std_population: float   # 总体标准差（除以 n）
    minimum: float
    maximum: float
    median: float
    range_: float


def summarise(values: Sequence[float]) -> StatsSummary:
    """计算均值、样本/总体标准差、极差、中位数。"""
    data = [float(v) for v in values]
    if not data:
        raise ChemToolboxError("数据为空，无法统计")

    n = len(data)
    mean = sum(data) / n

    if n > 1:
        std_sample = math.sqrt(sum((x - mean) ** 2 for x in data) / (n - 1))
    else:
        std_sample = 0.0
    std_population = math.sqrt(sum((x - mean) ** 2 for x in data) / n)

    ordered = sorted(data)
    if n % 2 == 1:
        median = ordered[n // 2]
    else:
        median = (ordered[n // 2 - 1] + ordered[n // 2]) / 2.0

    return StatsSummary(
        n=n,
        mean=mean,
        std_sample=std_sample,
        std_population=std_population,
        minimum=ordered[0],
        maximum=ordered[-1],
        median=median,
        range_=ordered[-1] - ordered[0],
    )


@dataclass
class LinearFit:
    """最小二乘线性拟合结果 y = slope·x + intercept。"""

    slope: float
    intercept: float
    r_squared: float
    n: int
    stderr_slope: float
    stderr_intercept: float


def linear_fit(xs: Sequence[float], ys: Sequence[float]) -> LinearFit:
    """最小二乘线性拟合，返回斜率、截距、R² 以及斜率和截距的标准误差。

    公式（最小二乘法）
    ------------------
    slope     = Σ((xi - x̄)(yi - ȳ)) / Σ((xi - x̄)²)
    intercept = ȳ - slope · x̄
    R²        = 1 - SS_res / SS_tot

    实验室常见用途：标准曲线（浓度 → 吸光度），斜率就是 ε·b。
    """
    x = [float(v) for v in xs]
    y = [float(v) for v in ys]
    if len(x) != len(y):
        raise ChemToolboxError(f"x 与 y 的数据点个数不一致：{len(x)} vs {len(y)}")
    n = len(x)
    if n < 3:
        raise ChemToolboxError("线性拟合至少需要 3 个数据点")

    mean_x = sum(x) / n
    mean_y = sum(y) / n

    sxx = sum((xi - mean_x) ** 2 for xi in x)
    if sxx == 0:
        raise ChemToolboxError("所有 x 值相同，无法拟合直线")

    sxy = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    slope = sxy / sxx
    intercept = mean_y - slope * mean_x

    residuals = [yi - (slope * xi + intercept) for xi, yi in zip(x, y)]
    ss_res = sum(r ** 2 for r in residuals)
    ss_tot = sum((yi - mean_y) ** 2 for yi in y)
    r_squared = 1.0 if ss_tot == 0 else 1.0 - ss_res / ss_tot

    # 残差标准误与参数标准误差（拟合了 2 个参数，自由度是 n-2）
    if n > 2:
        residual_se = math.sqrt(ss_res / (n - 2))
        stderr_slope = residual_se / math.sqrt(sxx)
        stderr_intercept = residual_se * math.sqrt(1.0 / n + mean_x ** 2 / sxx)
    else:
        stderr_slope = stderr_intercept = float("nan")

    return LinearFit(
        slope=slope,
        intercept=intercept,
        r_squared=r_squared,
        n=n,
        stderr_slope=stderr_slope,
        stderr_intercept=stderr_intercept,
    )


def concentration_from_curve(
    absorbance: float,
    slope: float,
    intercept: float,
) -> float:
    """由标准曲线 y = slope·x + intercept 反推未知样品的浓度 x。

    实验里的用法：先测一系列标准溶液的吸光度做拟合，再测未知样，
    把它的吸光度代进方程反解出浓度。
    """
    if slope == 0:
        raise ChemToolboxError("标准曲线斜率为 0，无法反推浓度")
    return (absorbance - intercept) / slope


# ===========================================================================
# 四、光谱数据处理：找峰与半高宽
# ===========================================================================

@dataclass
class Peak:
    """一个谱峰。"""

    index: int
    position: float          # 峰位（x 轴单位，如 nm 或 cm⁻¹）
    height: float            # 峰高（y 轴原始值）
    prominence: float        # 相对基线的凸出高度
    fwhm: float | None       # 半高全宽（在峰两侧按线性插值估计）


def find_peaks(
    xs: Sequence[float],
    ys: Sequence[float],
    threshold_sigma: float = 3.0,
    min_prominence: float | None = None,
) -> Tuple[List[Peak], float]:
    """在一维谱数据中找峰，返回 (峰列表, 噪声水平估计)。

    判定规则
    --------
    1. 局部极大：y[i] 同时大于左右邻居。
    2. 显著：峰相对基线的凸出量 ≥ max(threshold_sigma × 噪声, min_prominence)。
       噪声用相邻点差分的稳健估计：噪声 ≈ 1.4826 × median(|Δy|) / √2
       （系数 1.4826 把中位数绝对偏差折算成标准差；除以 √2 是因为
        相邻两点相减会把方差放大 2 倍。）
    3. 半高宽：在 y = 基线 + 凸出量/2 处向左右线性插值求交点，得到峰宽。

    为什么用中位数估基线而不是均值？—— 谱图里的峰本身会把均值抬高，
    中位数对少数强峰不敏感，估出来的基线更接近真实底噪。
    """
    x = [float(v) for v in xs]
    y = [float(v) for v in ys]
    if len(x) != len(y):
        raise ChemToolboxError(f"x 与 y 的数据点个数不一致：{len(x)} vs {len(y)}")
    if len(x) < 5:
        raise ChemToolboxError("找峰至少需要 5 个数据点")
    if not all(x[i] <= x[i + 1] for i in range(len(x) - 1)):
        raise ChemToolboxError("x 轴数据必须按从小到大排列（请先排序）")

    sorted_y = sorted(y)
    n = len(sorted_y)
    baseline = (
        sorted_y[n // 2] if n % 2 == 1
        else (sorted_y[n // 2 - 1] + sorted_y[n // 2]) / 2.0
    )

    # --- 噪声估计 ---
    diffs = sorted(abs(y[i + 1] - y[i]) for i in range(len(y) - 1))
    m = len(diffs)
    median_abs_diff = (
        diffs[m // 2] if m % 2 == 1
        else (diffs[m // 2 - 1] + diffs[m // 2]) / 2.0
    )
    noise = 1.4826 * median_abs_diff / math.sqrt(2.0)

    effective_min_prominence = threshold_sigma * noise
    if min_prominence is not None:
        effective_min_prominence = max(effective_min_prominence, min_prominence)

    peaks: List[Peak] = []
    for i in range(1, len(y) - 1):
        if not (y[i] > y[i - 1] and y[i] >= y[i + 1]):
            continue
        prominence = y[i] - baseline
        if prominence < effective_min_prominence:
            continue

        # --- 半高宽：在 y = 基线 + 凸出量/2 处左右插值 ---
        half_level = baseline + prominence / 2.0

        left_x = None
        for j in range(i, 0, -1):
            if y[j - 1] <= half_level <= y[j] and y[j] != y[j - 1]:
                ratio = (half_level - y[j - 1]) / (y[j] - y[j - 1])
                left_x = x[j - 1] + ratio * (x[j] - x[j - 1])
                break

        right_x = None
        for j in range(i, len(y) - 1):
            if y[j] >= half_level >= y[j + 1] and y[j] != y[j + 1]:
                ratio = (y[j] - half_level) / (y[j] - y[j + 1])
                right_x = x[j] + ratio * (x[j + 1] - x[j])
                break

        fwhm = None
        if left_x is not None and right_x is not None and right_x > left_x:
            fwhm = right_x - left_x

        peaks.append(Peak(index=i, position=x[i], height=y[i], prominence=prominence, fwhm=fwhm))

    return peaks, noise
