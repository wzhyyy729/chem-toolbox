"""chem_toolbox —— 化学专业学生的日常计算工具箱。

包结构
------
core.py     纯计算逻辑（无副作用，方便单独测试）
__main__.py 命令行界面（解析参数、打印结果）

两个模块刻意分开，是为了让核心算法能被单独导入和复用：
    >>> from chem_toolbox import analyse_formula
    >>> report = analyse_formula("C6H12O6")
    >>> round(report.molar_mass, 3)
    180.156
"""

from .core import (
    ATOMIC_MASS,
    KW_25C,
    ChemToolboxError,
    FormulaReport,
    LinearFit,
    StatsSummary,
    analyse_formula,
    concentration_from_curve,
    dilute,
    ka_from_pka,
    kb_from_ka,
    linear_fit,
    mass_for_molarity,
    mass_from_absorbance,
    molarity,
    moles,
    parse_formula,
    ph_strong_acid,
    ph_strong_base,
    ph_weak_acid,
    ph_weak_base,
    summarise,
)

__version__ = "0.1.0"

__all__ = [
    "ATOMIC_MASS",
    "KW_25C",
    "ChemToolboxError",
    "FormulaReport",
    "LinearFit",
    "StatsSummary",
    "analyse_formula",
    "concentration_from_curve",
    "dilute",
    "ka_from_pka",
    "kb_from_ka",
    "linear_fit",
    "mass_for_molarity",
    "mass_from_absorbance",
    "molarity",
    "moles",
    "parse_formula",
    "ph_strong_acid",
    "ph_strong_base",
    "ph_weak_acid",
    "ph_weak_base",
    "summarise",
    "__version__",
]
