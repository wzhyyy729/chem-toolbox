"""chem_toolbox —— 化学专业学生的日常计算工具箱。

三个子命令
----------
formula   解析分子式：摩尔质量、元素质量分数、不饱和度
solution  溶液配制、稀释、pH（强/弱酸弱碱）、朗伯-比尔定律
stats     实验数据统计：均值、标准差、线性拟合、R²

所有计算只用 Python 标准库，无需安装任何第三方包。

用法示例
--------
python -m chem_toolbox formula "CuSO4.5H2O"
python -m chem_toolbox solution ph --acid-kind weak --pka 4.76 --conc 0.1
python -m chem_toolbox stats fit --file data.csv --x-col 0 --y-col 1
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import List, Sequence, Tuple

from . import core

# Windows 控制台默认编码是 GBK，直接打印中文可能抛 UnicodeEncodeError。
# 统一切成 UTF-8 输出，遇到无法编码的字符退化成替代字符而不是崩溃。
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

__version__ = "0.1.0"


# ---------------------------------------------------------------------------
# 工具函数
# ---------------------------------------------------------------------------

def _fmt(value: float | None, digits: int = 4) -> str:
    """把数字格式化成便于阅读的字符串；None / NaN 显示为 '—'。"""
    if value is None:
        return "—"
    if isinstance(value, float) and math.isnan(value):
        return "—"
    if value != 0 and (abs(value) < 1e-4 or abs(value) >= 1e6):
        return f"{value:.{digits}e}"
    return f"{value:.{digits}f}"


def _load_table(path: str, delimiter: str | None = None) -> List[List[str]]:
    """读取文本数据表（CSV / TSV / 空白分隔都支持），返回二维字符串列表。

    会自动跳过空行和以 # 开头的注释行——实验仪器导出的文件经常带表头注释。
    """
    file_path = Path(path)
    if not file_path.exists():
        raise core.ChemToolboxError(f"找不到文件：{path}")

    text = file_path.read_text(encoding="utf-8-sig", errors="replace")
    lines = [
        line for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not lines:
        raise core.ChemToolboxError(f"文件 {path} 里没有有效数据行")

    if delimiter is None:
        # 自动猜分隔符：逗号 / 制表符 / 分号 / 空白
        sample = lines[0]
        if "," in sample:
            delimiter = ","
        elif "\t" in sample:
            delimiter = "\t"
        elif ";" in sample:
            delimiter = ";"
        else:
            delimiter = None  # None 表示按任意空白切分

    if delimiter is None:
        rows = [line.split() for line in lines]
    else:
        rows = list(csv.reader(lines, delimiter=delimiter))

    return [row for row in rows if row]


def _column(rows: Sequence[Sequence[str]], index: int, label: str) -> List[float]:
    """取出某一列并转成浮点数，遇到非数字（比如表头）就跳过。"""
    values: List[float] = []
    for row_number, row in enumerate(rows, start=1):
        if index >= len(row):
            continue
        cell = row[index].strip()
        if not cell:
            continue
        try:
            values.append(float(cell))
        except ValueError:
            # 表头或单位行，跳过即可
            continue
    if not values:
        raise core.ChemToolboxError(
            f"第 {index} 列（{label}）里没有解析到任何数字，请检查 --x-col / --y-col 是否选对"
        )
    return values


def _emit(payload: dict, as_json: bool) -> None:
    """统一输出：--json 给机器读，否则打印人类可读的文本。"""
    if as_json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))


# ---------------------------------------------------------------------------
# 子命令一：formula
# ---------------------------------------------------------------------------

def cmd_formula(args: argparse.Namespace) -> int:
    report = core.analyse_formula(args.formula)

    if args.json:
        _emit(
            {
                "formula": report.formula,
                "counts": report.counts,
                "molar_mass": report.molar_mass,
                "mass_fractions": report.mass_fractions,
                "dbe": report.dbe,
                "dbe_note": report.dbe_note,
            },
            True,
        )
        return 0

    print(f"分子式        : {report.formula}")
    composition = "  ".join(
        f"{element}{count}" for element, count in
        sorted(report.counts.items(), key=lambda kv: (-kv[1], kv[0]))
    )
    print(f"原子组成      : {composition}")
    print(f"摩尔质量      : {_fmt(report.molar_mass)} g/mol")
    print()
    print("元素质量分数")
    print(f"  {'元素':<6}{'原子数':>6}{'质量分数':>14}{'质量占比':>12}")
    for element, count in sorted(
        report.counts.items(), key=lambda kv: -report.mass_fractions[kv[0]]
    ):
        fraction = report.mass_fractions[element]
        print(
            f"  {element:<6}{count:>6}{fraction:>14.6f}{fraction * 100:>11.2f}%"
        )
    print()
    if report.dbe is not None:
        print(f"不饱和度 DBE  : {_fmt(report.dbe)}")
        if report.dbe < 0:
            print("                ⚠ 负值说明该组成不合法（氢或卤素太多了）")
        elif abs(report.dbe - round(report.dbe)) > 1e-9:
            print("                ⚠ 非整数：常规中性分子应得到整数，请检查分子式是否写错")
    elif report.dbe_note:
        print(f"不饱和度 DBE  : 不适用 —— {report.dbe_note}")
    return 0


# ---------------------------------------------------------------------------
# 子命令二：solution
# ---------------------------------------------------------------------------

def cmd_solution(args: argparse.Namespace) -> int:
    action = args.action

    # --- 配制：由质量求浓度，或由目标浓度求所需质量 ---
    if action == "prepare":
        if args.formula:
            molar_mass = core.analyse_formula(args.formula).molar_mass
            mass_source = f"由分子式 {args.formula} 算出"
        elif args.molar_mass:
            molar_mass = args.molar_mass
            mass_source = "由 --molar-mass 指定"
        else:
            raise core.ChemToolboxError("请用 --formula 或 --molar-mass 提供摩尔质量")

        if args.volume is None:
            raise core.ChemToolboxError("请用 --volume 指定溶液体积（单位 L）")

        if args.mass is not None:
            concentration = core.molarity(args.mass, molar_mass, args.volume)
            payload = {
                "action": "prepare",
                "molar_mass": molar_mass,
                "mass_g": args.mass,
                "volume_l": args.volume,
                "concentration_mol_per_l": concentration,
            }
            if args.json:
                _emit(payload, True)
                return 0
            print(f"摩尔质量      : {_fmt(molar_mass)} g/mol（{mass_source}）")
            print(f"称取质量      : {_fmt(args.mass)} g")
            print(f"定容体积      : {_fmt(args.volume)} L")
            print(f"物质的量浓度  : {_fmt(concentration)} mol/L")
        else:
            mass_needed = core.mass_for_molarity(args.conc, args.volume, molar_mass)
            payload = {
                "action": "prepare",
                "molar_mass": molar_mass,
                "target_concentration": args.conc,
                "volume_l": args.volume,
                "mass_needed_g": mass_needed,
            }
            if args.json:
                _emit(payload, True)
                return 0
            print(f"摩尔质量      : {_fmt(molar_mass)} g/mol（{mass_source}）")
            print(f"目标浓度      : {_fmt(args.conc)} mol/L")
            print(f"定容体积      : {_fmt(args.volume)} L")
            print(f"需称取质量    : {_fmt(mass_needed)} g")
        return 0

    # --- 稀释：c1V1 = c2V2 ---
    if action == "dilute":
        # 允许用户用 --volume-from / --volume-to 这类更直观的写法
        c1 = args.c1 if args.c1 is not None else args.conc
        result = core.dilute(c1, args.v1, args.c2, args.v2)
        payload = {"action": "dilute", **{k: v for k, v in result.items()}}
        if args.json:
            _emit(payload, True)
            return 0
        print("稀释定律 c1·V1 = c2·V2")
        for key, value in result.items():
            label = {"c1": "起始浓度 c1", "v1": "起始体积 V1", "c2": "最终浓度 c2", "v2": "最终体积 V2"}[key]
            unit = "mol/L" if key.startswith("c") else "L"
            print(f"  求得 {label:<12} = {_fmt(value)} {unit}")
        if result.get("v2") is not None and args.v1 is not None:
            print(f"  需加溶剂        = {_fmt(result['v2'] - args.v1)} L")
        return 0

    # --- pH 计算 ---
    if action == "ph":
        if args.conc is None or args.conc <= 0:
            raise core.ChemToolboxError("请用 --conc 指定浓度（mol/L，需大于 0）")

        payload: dict = {"action": "ph", "kind": args.acid_kind, "concentration": args.conc}

        if args.acid_kind == "strong-acid":
            ph = core.ph_strong_acid(args.conc)
            payload["ph"] = ph
            explanation = "强酸完全电离，[H+] = c，pH = -lg c"
        elif args.acid_kind == "strong-base":
            ph = core.ph_strong_base(args.conc)
            payload["ph"] = ph
            explanation = "强碱完全电离，[OH-] = c，pOH = -lg c，pH = 14 - pOH"
        elif args.acid_kind == "weak-acid":
            ka = args.ka if args.ka is not None else (
                core.ka_from_pka(args.pka) if args.pka is not None else None
            )
            if ka is None:
                raise core.ChemToolboxError("弱酸请用 --ka 或 --pka 提供电离常数")
            ph, alpha = core.ph_weak_acid(ka, args.conc)
            payload.update({"ka": ka, "ph": ph, "degree_of_ionisation": alpha})
            explanation = "一元弱酸精确解：x² + Ka·x - Ka·c = 0，[H+] = x"
        else:  # weak-base
            kb = args.kb if args.kb is not None else (
                core.kb_from_ka(core.ka_from_pka(args.pka)) if args.pka is not None else None
            )
            if kb is None:
                raise core.ChemToolboxError("弱碱请用 --kb 提供 Kb，或用 --pka 提供其共轭酸的 pKa")
            ph, alpha = core.ph_weak_base(kb, args.conc)
            payload.update({"kb": kb, "ph": ph, "degree_of_ionisation": alpha})
            explanation = "一元弱碱精确解：x² + Kb·x - Kb·c = 0，[OH-] = x"

        # 强酸强碱忽略水电离是有前提的，这里主动提醒
        warnings: List[str] = []
        if args.acid_kind in ("strong-acid", "strong-base") and args.conc < 1e-6:
            warnings.append("浓度低于 1e-6 mol/L，水的电离已不可忽略，本结果偏差较大")
        if args.acid_kind == "weak-acid" and payload.get("degree_of_ionisation", 0) > 0.05:
            warnings.append(
                f"电离度已达 {payload['degree_of_ionisation'] * 100:.1f}%，"
                f"用最简式 √(Ka·c) 会有明显误差；本结果用的是精确解，可直接使用"
            )

        payload["warnings"] = warnings
        if args.json:
            _emit(payload, True)
            return 0

        print(f"溶液类型      : {args.acid_kind}")
        print(f"浓度          : {_fmt(args.conc)} mol/L")
        if "ka" in payload:
            print(f"Ka            : {_fmt(payload['ka'])}")
        if "kb" in payload:
            print(f"Kb            : {_fmt(payload['kb'])}")
        print(f"pH            : {_fmt(payload['ph'], 3)}")
        if "degree_of_ionisation" in payload:
            print(f"电离度 α      : {payload['degree_of_ionisation'] * 100:.3f}%")
        print(f"依据          : {explanation}")
        for warning in warnings:
            print(f"⚠ 注意        : {warning}")
        return 0

    # --- 由标准曲线反推未知样品浓度 ---
    if action == "unknown":
        if args.absorbance is None:
            raise core.ChemToolboxError("请用 --absorbance 提供未知样品的吸光度")
        if args.slope is None:
            raise core.ChemToolboxError("请用 --slope 提供标准曲线斜率（可由 stats fit 得到）")
        intercept = args.intercept if args.intercept is not None else 0.0
        concentration = core.concentration_from_curve(args.absorbance, args.slope, intercept)
        payload = {
            "action": "unknown",
            "absorbance": args.absorbance,
            "slope": args.slope,
            "intercept": intercept,
            "concentration": concentration,
        }
        if args.json:
            _emit(payload, True)
            return 0
        print("由标准曲线反推未知样品")
        print(f"  标准曲线      : y = {_fmt(args.slope)}·x + {_fmt(intercept)}")
        print(f"  实测吸光度 A  = {_fmt(args.absorbance)}")
        print(f"  求得浓度 c    = {_fmt(concentration)}")
        if concentration < 0:
            print("  ⚠ 注意        : 算出负浓度，说明吸光度低于曲线截距，样品可能低于检出限")
        return 0

    # --- 朗伯-比尔定律 ---
    if action == "beer":
        if args.epsilon and args.absorbance is not None and args.path:
            concentration = core.mass_from_absorbance(args.absorbance, args.epsilon, args.path)
            payload = {
                "action": "beer",
                "absorbance": args.absorbance,
                "epsilon": args.epsilon,
                "path_cm": args.path,
                "concentration_mol_per_l": concentration,
            }
            if args.json:
                _emit(payload, True)
                return 0
            print("朗伯-比尔定律 A = ε·b·c")
            print(f"  吸光度 A      = {_fmt(args.absorbance)}")
            print(f"  摩尔吸光系数 ε = {_fmt(args.epsilon)} L·mol⁻¹·cm⁻¹")
            print(f"  光程 b        = {_fmt(args.path)} cm")
            print(f"  求得浓度 c    = {_fmt(concentration)} mol/L")
        else:
            raise core.ChemToolboxError(
                "求浓度请同时提供 --absorbance、--epsilon、--path；"
                "若想由标准曲线求 ε，请用：stats fit --file 数据.csv"
            )
        return 0

    raise core.ChemToolboxError(f"未知的 solution 动作：{action}")


# ---------------------------------------------------------------------------
# 子命令三：stats
# ---------------------------------------------------------------------------

def cmd_stats(args: argparse.Namespace) -> int:
    rows = _load_table(args.file, args.delimiter)

    if args.action == "summary":
        values = _column(rows, args.x_col, "数据列")
        summary = core.summarise(values)
        payload = {
            "action": "summary",
            "n": summary.n,
            "mean": summary.mean,
            "std_sample": summary.std_sample,
            "std_population": summary.std_population,
            "median": summary.median,
            "min": summary.minimum,
            "max": summary.maximum,
            "range": summary.range_,
        }
        if args.json:
            _emit(payload, True)
            return 0
        print(f"数据点数 n       : {summary.n}")
        print(f"平均值  x̄        : {_fmt(summary.mean)}")
        print(f"样本标准差 s     : {_fmt(summary.std_sample)}   ← 实验报告通常用这个（除以 n-1）")
        print(f"总体标准差 σ     : {_fmt(summary.std_population)}   ← 除以 n")
        print(f"中位数           : {_fmt(summary.median)}")
        print(f"最小值 / 最大值  : {_fmt(summary.minimum)} / {_fmt(summary.maximum)}")
        print(f"极差 R           : {_fmt(summary.range_)}")
        if summary.n > 1:
            sem = summary.std_sample / math.sqrt(summary.n)
            print(f"平均值的标准误   : {_fmt(sem)}")
            print(f"结果可写为       : {_fmt(summary.mean)} ± {_fmt(sem)}")
        return 0

    if args.action == "fit":
        xs = _column(rows, args.x_col, "x 列")
        ys = _column(rows, args.y_col, "y 列")
        fit = core.linear_fit(xs, ys)
        payload = {
            "action": "fit",
            "n": fit.n,
            "slope": fit.slope,
            "intercept": fit.intercept,
            "r_squared": fit.r_squared,
            "stderr_slope": fit.stderr_slope,
            "stderr_intercept": fit.stderr_intercept,
        }
        if args.json:
            _emit(payload, True)
            return 0
        print("最小二乘线性拟合  y = k·x + b")
        print(f"  数据点数 n    = {fit.n}")
        print(f"  斜率 k        = {_fmt(fit.slope)} ± {_fmt(fit.stderr_slope)}")
        print(f"  截距 b        = {_fmt(fit.intercept)} ± {_fmt(fit.stderr_intercept)}")
        print(f"  相关系数 R²   = {fit.r_squared:.6f}")
        if fit.r_squared >= 0.999:
            print("  线性评价      : 优秀（R² ≥ 0.999）")
        elif fit.r_squared >= 0.99:
            print("  线性评价      : 良好（R² ≥ 0.99）")
        else:
            print("  ⚠ 线性评价    : 偏差较大，检查是否有离群点或超出线性范围")
        print()
        print(f"  拟合方程      : y = {_fmt(fit.slope)}·x + {_fmt(fit.intercept)}")
        if args.label_x and args.label_y:
            print(f"  物理含义      : {args.label_y} 对 {args.label_x} 作图")
        return 0

    raise core.ChemToolboxError(f"未知的 stats 动作：{args.action}")


# ---------------------------------------------------------------------------
# 命令行入口
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="chem_toolbox",
        description="化学专业学生的日常计算工具箱（分子式 / 溶液 / 数据）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例：\n"
            "  python -m chem_toolbox formula \"CuSO4.5H2O\"\n"
            "  python -m chem_toolbox formula \"Ca(OH)2\"     注意 Windows 下给括号加引号\n"
            "  python -m chem_toolbox solution ph --acid-kind weak-acid --pka 4.76 --conc 0.1\n"
            "  python -m chem_toolbox stats fit --file data.csv --x-col 0 --y-col 1\n"
        ),
    )
    parser.add_argument("--version", action="version", version=f"chem_toolbox {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # formula
    p_formula = subparsers.add_parser("formula", help="解析分子式")
    p_formula.add_argument("formula", help="化学式，如 C6H12O6、Ca(OH)2、CuSO4.5H2O")
    p_formula.add_argument("--json", action="store_true", help="以 JSON 输出，便于程序调用")
    p_formula.set_defaults(func=cmd_formula)

    # solution
    p_solution = subparsers.add_parser("solution", help="溶液配制 / 稀释 / pH / 比尔定律")
    solution_sub = p_solution.add_subparsers(dest="action", required=True)

    p_prepare = solution_sub.add_parser("prepare", help="配制溶液：质量↔浓度换算")
    p_prepare.add_argument("--formula", help="溶质分子式，如 NaCl（自动算摩尔质量）")
    p_prepare.add_argument("--molar-mass", type=float, help="或直接给摩尔质量 g/mol")
    p_prepare.add_argument("--mass", type=float, help="已称取的质量 g（给了就求浓度）")
    p_prepare.add_argument("--conc", type=float, default=0.0, help="目标浓度 mol/L（没给 --mass 时用它求需称质量）")
    p_prepare.add_argument("--volume", type=float, help="溶液体积 L")
    p_prepare.add_argument("--json", action="store_true")
    p_prepare.set_defaults(func=cmd_solution)

    p_dilute = solution_sub.add_parser("dilute", help="稀释定律 c1V1 = c2V2，留一个未知不填")
    p_dilute.add_argument("--c1", type=float, help="起始浓度 mol/L")
    p_dilute.add_argument("--conc", type=float, help="同 --c1（二选一）")
    p_dilute.add_argument("--v1", type=float, help="起始体积 L")
    p_dilute.add_argument("--c2", type=float, help="最终浓度 mol/L")
    p_dilute.add_argument("--v2", type=float, help="最终体积 L")
    p_dilute.add_argument("--json", action="store_true")
    p_dilute.set_defaults(func=cmd_solution)

    p_ph = solution_sub.add_parser("ph", help="pH 计算（强酸/强碱/弱酸/弱碱）")
    p_ph.add_argument("--acid-kind", required=True,
                      choices=["strong-acid", "strong-base", "weak-acid", "weak-base"])
    p_ph.add_argument("--conc", type=float, required=True, help="浓度 mol/L")
    p_ph.add_argument("--ka", type=float, help="弱酸电离常数 Ka")
    p_ph.add_argument("--kb", type=float, help="弱碱电离常数 Kb")
    p_ph.add_argument("--pka", type=float, help="pKa（弱酸用；弱碱时视为共轭酸的 pKa）")
    p_ph.add_argument("--json", action="store_true")
    p_ph.set_defaults(func=cmd_solution)

    p_beer = solution_sub.add_parser("beer", help="朗伯-比尔定律 A = εbc")
    p_beer.add_argument("--absorbance", type=float, help="吸光度 A")
    p_beer.add_argument("--epsilon", type=float, help="摩尔吸光系数 ε（L·mol⁻¹·cm⁻¹）")
    p_beer.add_argument("--path", type=float, default=1.0, help="光程 b（cm，默认 1）")
    p_beer.add_argument("--json", action="store_true")
    p_beer.set_defaults(func=cmd_solution)

    p_unknown = solution_sub.add_parser("unknown", help="由标准曲线反推未知样品浓度")
    p_unknown.add_argument("--absorbance", type=float, required=True, help="未知样品的吸光度 A")
    p_unknown.add_argument("--slope", type=float, required=True, help="标准曲线斜率（用 stats fit 得到）")
    p_unknown.add_argument("--intercept", type=float, default=0.0, help="标准曲线截距（默认 0）")
    p_unknown.add_argument("--json", action="store_true")
    p_unknown.set_defaults(func=cmd_solution)

    # stats
    p_stats = subparsers.add_parser("stats", help="实验数据统计与线性拟合")
    stats_sub = p_stats.add_subparsers(dest="action", required=True)

    p_summary = stats_sub.add_parser("summary", help="均值 / 标准差 / 极差")
    p_summary.add_argument("--file", required=True, help="数据文件（CSV/TXT）")
    p_summary.add_argument("--x-col", type=int, default=0, help="数据所在列号，从 0 开始")
    p_summary.add_argument("--delimiter", help="自定义分隔符，如 ',' 或 '\\t'")
    p_summary.add_argument("--json", action="store_true")
    p_summary.set_defaults(func=cmd_stats)

    p_fit = stats_sub.add_parser("fit", help="最小二乘线性拟合与 R²")
    p_fit.add_argument("--file", required=True, help="数据文件（CSV/TXT）")
    p_fit.add_argument("--x-col", type=int, default=0, help="x 列号，从 0 开始")
    p_fit.add_argument("--y-col", type=int, default=1, help="y 列号，从 0 开始")
    p_fit.add_argument("--delimiter", help="自定义分隔符")
    p_fit.add_argument("--label-x", help="x 轴物理量名称，仅用于打印说明")
    p_fit.add_argument("--label-y", help="y 轴物理量名称")
    p_fit.add_argument("--json", action="store_true")
    p_fit.set_defaults(func=cmd_stats)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except core.ChemToolboxError as error:
        print(f"错误：{error}", file=sys.stderr)
        return 2
    except FileNotFoundError as error:
        print(f"错误：找不到文件 {error.filename}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("\n已取消", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
