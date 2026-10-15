# Chem Toolbox · 化学计算工具箱

[English](README.md) | **简体中文**

[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Dependencies](https://img.shields.io/badge/dependencies-none-3FB950)](#安装)
[![Tests](https://img.shields.io/badge/tests-87%20passed-3FB950)](#自检测试)
[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

化学专业学生日常要做的三类计算，收进一个命令行工具里：**分子式分析、溶液配制与 pH、实验数据统计与曲线拟合**。

只用 Python 标准库，**零第三方依赖** —— 拷到任何装了 Python 的电脑上都能直接跑，不用 pip install。

> 这个项目由化学专业大二学生发起，代码在 AI 协作下产出，但**每一个公式都经过手工核对，并有 87 项自动化断言覆盖**。写得不算花哨，但算出来的数字你可以拿去交实验报告。
>
> 方向聚焦**有机化学、物理化学、药学与材料**；分析化学方向的谱图处理工具已移除（详见[开发记录](#开发记录)）。

---

## 目录

- [能做什么](#能做什么)
- [为什么做这个](#为什么做这个)
- [安装](#安装)
- [使用示例](#使用示例)
- [背后的化学与数学](#背后的化学与数学)
- [自检测试](#自检测试)
- [已知限制](#已知限制)
- [目录结构](#目录结构)
- [开发记录](#开发记录)
- [License](#license)

---

## 能做什么

| 命令 | 功能 |
| :--- | :--- |
| `formula` | 解析化学式 → 摩尔质量、元素质量分数、不饱和度（DBE） |
| `solution prepare` | 配制溶液：质量 ↔ 浓度互算（可直接给分子式，自动查摩尔质量） |
| `solution dilute` | 稀释定律 c₁V₁ = c₂V₂，任意留一个未知 |
| `solution ph` | pH 计算：强酸、强碱、一元弱酸、一元弱碱（精确解，非最简式） |
| `solution beer` | 朗伯-比尔定律 A = εbc |
| `solution unknown` | 由标准曲线反推未知样品浓度 |
| `stats summary` | 均值、样本/总体标准差、标准误、极差、中位数 |
| `stats fit` | 最小二乘线性拟合：斜率、截距、R²、参数标准误差 |

**特别支持**：

- 嵌套括号化学式 `K3[Fe(CN)6]`、水合物 `CuSO4.5H2O`、前置系数 `2H2O`
- **79 种元素的相对原子质量**，覆盖有机、物化、药学（Pt/Sn/Sb/Te/As/Se）与材料（Zr/Hf/Nb/Ta/W/Mo/稀土/In/Ga/Ge）
- 带注释行的仪器导出数据（自动跳过 `#` 开头与空行），自动识别逗号/制表符/空格分隔

---

## 为什么做这个

做实验报告时反复遇到同一批计算：配溶液要算称多少克、平行滴定要算标准差、标准曲线要算 R² 和未知样浓度。用计算器一个个敲，既慢又容易抄错。

更关键的是：**这些计算里的每一个公式都有适用条件**。比如强酸 pH 在浓度低于 10⁻⁶ mol/L 时不能忽略水的电离；弱酸 pH 用最简式 √(Kₐc) 在电离度超过 5% 时会明显偏差。很多脚本和在线计算器不告诉你这些，只是闷头给个数。

所以这个工具的做法是：**给结果的同时给依据，并在前提不成立时主动警告**。

---

## 安装

需要 Python 3.10 或更高版本（用到了 `X | None` 这种类型标注语法）。

```bash
git clone https://github.com/wzhyyy729/chem-toolbox.git
cd chem-toolbox
python -m chem_toolbox --help
```

没有依赖要装。不需要 `pip install`，没有 `requirements.txt`。

---

## 使用示例

### 分子式分析

```console
$ python -m chem_toolbox formula "CuSO4.5H2O"
分子式        : CuSO4.5H2O
原子组成      : H10  O9  Cu1  S1
摩尔质量      : 249.6770 g/mol

元素质量分数
  元素       原子数          质量分数        质量占比
  O          9      0.576709      57.67%
  Cu         1      0.254513      25.45%
  S          1      0.128406      12.84%
  H         10      0.040372       4.04%

不饱和度 DBE  : 不适用 —— 分子含 Cu、S，不饱和度公式仅适用于 C/H/N/O/P/S/卤素组成的有机物
```

嵌套括号也认：

```console
$ python -m chem_toolbox formula "K3[Fe(CN)6]"
摩尔质量      : 329.2470 g/mol
```

有机物的不饱和度会正常算，硫和磷也计入：

```console
$ python -m chem_toolbox formula "C6H6"
不饱和度 DBE  : 4.0000        # 苯：3 个双键 + 1 个环

$ python -m chem_toolbox formula "C4H4S"
不饱和度 DBE  : 3.0000        # 噻吩：2 个双键 + 1 个环

$ python -m chem_toolbox formula "P(C6H5)3"
不饱和度 DBE  : 12.0000       # 三苯基膦：3 个苯环 × 4
```

### pH 计算

```console
$ python -m chem_toolbox solution ph --acid-kind weak-acid --pka 4.76 --conc 0.1
溶液类型      : weak-acid
浓度          : 0.1000 mol/L
Ka            : 1.7378e-05
pH            : 2.883
电离度 α      : 1.310%
依据          : 一元弱酸精确解：x² + Ka·x - Ka·c = 0，[H+] = x
```

前提不成立时会主动警告：

```console
$ python -m chem_toolbox solution ph --acid-kind strong-acid --conc 1e-8
pH            : 8.000
依据          : 强酸完全电离，[H+] = c，pH = -lg c
⚠ 注意        : 浓度低于 1e-6 mol/L，水的电离已不可忽略，本结果偏差较大
```

### 配制溶液

```console
$ python -m chem_toolbox solution prepare --formula NaCl --conc 0.1 --volume 0.5
摩尔质量      : 58.4400 g/mol（由分子式 NaCl 算出）
目标浓度      : 0.1000 mol/L
定容体积      : 0.5000 L
需称取质量    : 2.9220 g
```

### 标准曲线拟合

```console
$ python -m chem_toolbox stats fit --file examples/standard_curve.csv --x-col 0 --y-col 1
最小二乘线性拟合  y = k·x + b
  数据点数 n    = 8
  斜率 k        = 1.9957 ± 0.0096
  截距 b        = 0.0020 ± 0.0024
  相关系数 R²   = 0.999861
  线性评价      : 优秀（R² ≥ 0.999）
```

### 脚本调用

除了命令行，也可以当普通 Python 库用：

```python
from chem_toolbox import analyse_formula, ph_weak_acid, linear_fit

report = analyse_formula("C6H12O6")
print(report.molar_mass)            # 180.156

ph, alpha = ph_weak_acid(ka=1.75e-5, concentration=0.1)
print(round(ph, 4))                 # 2.8813

fit = linear_fit([0.1, 0.2, 0.3], [0.2, 0.4, 0.6])
print(fit.slope, fit.r_squared)     # 2.0 1.0
```

加 `--json` 可以让任何子命令输出 JSON，方便接到别的程序里：

```console
$ python -m chem_toolbox formula "H2O" --json
{
  "formula": "H2O",
  "counts": {"H": 2, "O": 1},
  "molar_mass": 18.015,
  "mass_fractions": {"H": 0.111899, "O": 0.888101},
  "dbe": 0.0,
  "dbe_note": ""
}
```

---

## 背后的化学与数学

这一节解释每个功能用的公式和它的适用条件 —— 也是我准备专业评估测试时整理的内容。

### 不饱和度 DBE

$$\mathrm{DBE} = \frac{2C + 2 + N + P - H - X}{2}$$

其中 C、N、P、H 分别是碳/氮/磷/氢原子数，X 是卤素原子总数；**氧和硫不计入**（二价，对不饱和度没有贡献）。

- 适用条件：只适用于 **C、H、N、O、P、S、卤素**组成的有机物。含金属或其他元素的式子（如 `K3[Fe(CN)6]`、`CuSO4`）公式不成立，程序会拒绝计算而不是给个错数。
- **磷要算进去**：三苯基膦 P(C₆H₅)₃ = C₁₈H₁₅P，代入得 (36+2+1−15)/2 = 12，正好等于 3 个苯环 × 4；漏掉磷会算出 11.5 这种非整数。
- **硫不用算**：噻吩 C₄H₄S 代入得 (8+2−4)/2 = 3，正确；若把硫也当三价加上去就会变成 3.5。
- 非整数结果说明分子式写错了（或该结构是盐/离子型）；负值说明氢或卤素过多，组成不合法。

### 一元弱酸 pH（精确解）

设电离出的 $[\mathrm{H^+}] = x$，由物料守恒与电荷守恒：

$$K_a = \frac{x^2}{c - x} \quad\Longrightarrow\quad x^2 + K_a x - K_a c = 0 \quad\Longrightarrow\quad x = \frac{-K_a + \sqrt{K_a^2 + 4K_a c}}{2}$$

- 对比最简式 $x \approx \sqrt{K_a c}$：当电离度 α = x/c 超过 5% 时最简式偏差明显，本工具用精确解并会提示。
- 强酸强碱用 $[\mathrm{H^+}] = c$，前提是 **c ≥ 10⁻⁶ mol/L**；更稀时水的电离不可忽略，程序会警告。
- 弱碱完全对称：$x^2 + K_b x - K_b c = 0$，再由 $\mathrm{pOH} = -\lg x$、$\mathrm{pH} = \mathrm{p}K_w - \mathrm{pOH}$ 换算。
- 共轭酸碱对：$K_a \cdot K_b = K_w$（25 ℃ 时 $K_w = 1.0\times10^{-14}$）。

### 最小二乘拟合

$$k = \frac{\sum (x_i - \bar{x})(y_i - \bar{y})}{\sum (x_i - \bar{x})^2}, \qquad b = \bar{y} - k\bar{x}, \qquad R^2 = 1 - \frac{SS_{\text{res}}}{SS_{\text{tot}}}$$

参数标准误差由残差标准误 $s = \sqrt{SS_{\text{res}}/(n-2)}$ 推得：

$$s_k = \frac{s}{\sqrt{\sum (x_i-\bar{x})^2}}, \qquad s_b = s\sqrt{\frac{1}{n} + \frac{\bar{x}^2}{\sum (x_i-\bar{x})^2}}$$

---

## 自检测试

```console
$ python -m chem_toolbox.tests
==============================================================
chem_toolbox 自检测试
==============================================================
  [通过] 分子式解析
  [通过] 药学与材料分子
  [通过] 溶液与 pH
  [通过] 统计与拟合
--------------------------------------------------------------
断言通过 87 项，失败 0 项
全部通过
```

测试特点：**每个期望值都能用教科书、药典或手算核对**，不是"跑起来不报错就算过"。例如：

- 葡萄糖 180.156 g/mol、胆矾 249.68 g/mol、`K3[Fe(CN)6]` 329.24 g/mol
- 阿司匹林 180.158、对乙酰氨基酚 151.163、布洛芬 206.281、顺铂 300.05
- 二氧化钛 79.866、钛酸钡 233.192、氧化锆 123.222、碳化钨 195.85
- 苯 DBE = 4、环己烷 = 1、乙烷 = 0、噻吩 = 3、三苯基膦 = 12、三甲基膦 = 0；含金属时拒绝计算
- 0.1 mol/L 醋酸（pKₐ 4.76）pH = 2.8813 ± 0.0005；0.1 mol/L 氨水（K_b 1.8×10⁻⁵）pH = 11.1248
- 非法输入（空分子式、小写元素符号、括号不匹配、负浓度、x 全相同）必须抛错，而不是悄悄算错

不需要 pytest，直接 `python -m chem_toolbox.tests` 就能跑。

**自动验证**：每次推送都会由 GitHub Actions 在 **Ubuntu / Windows / macOS × Python 3.10–3.13** 上重复跑一遍这些断言，运行记录公开可查：[Actions](https://github.com/wzhyyy729/chem-toolbox/actions)。这意味着"测试通过"不是我自己说的，是 GitHub 说的。

> 为什么要专门测 Windows：项目里处理过 Windows 控制台默认 GBK 编码导致的 `UnicodeEncodeError`（打印中文和 `✓` 时会崩）。这类问题只在 Windows 上暴露，所以在 CI 里固定保留 Windows 环境。

---

## 已知限制

写清楚边界，比假装什么都能算更负责：

- **pH 只支持一元酸碱**。多元酸（如 H₃PO₄、H₂CO₃）需要逐级电离常数联立求解，暂未实现。
- **不计算活度系数**，用浓度代替活度。高离子强度（如 ≥ 0.1 mol/L 的强电解质）下结果会有偏差。
- **不处理氧化还原与电化学计算**（能斯特方程、电极电势），暂未实现。
- **不做谱图处理**（找峰、峰面积、基线校正）。原因见下节开发记录。
- **不画图**。这是刻意取舍：绘图会引入 matplotlib 依赖，破坏"零依赖"这个最大的便利性。需要作图请把 CSV 丢给 Excel 或 Origin。

---

## 目录结构

```
chem-toolbox/
├── .github/
│   └── workflows/
│       └── tests.yml    CI：三系统 × 四版本 Python 自动跑 87 项断言
├── chem_toolbox/
│   ├── __init__.py     对外接口（方便 from chem_toolbox import ...）
│   ├── __main__.py     命令行界面：参数解析与结果打印
│   ├── core.py         纯计算逻辑，无打印、无文件读写，便于单独测试
│   └── tests.py        87 项自检断言，无需 pytest
├── examples/
│   ├── standard_curve.csv    分光光度法标准曲线示例
│   └── titration.csv         平行滴定数据示例
├── LICENSE
└── README.md
```

`core.py` 与 `__main__.py` 刻意分开：算法不含任何输入输出，所以既能被命令行调用，也能被别的脚本 import，还能被测试直接断言。

---

## 开发记录

**2026-10-16 · v0.1.0**

- 首个版本：分子式分析、溶液与 pH、统计与拟合。
- 元素表从最初约 40 种扩到 **79 种**，补上药学与材料方向常用的 Pt、Sn、Sb、Te、As、Se、Zr、Hf、Nb、Ta、W、Mo、稀土、In、Ga、Ge。
- **修正了一个化学错误**：不饱和度公式原先漏掉了磷。三苯基膦 P(C₆H₅)₃ 用旧公式会得到 11.5（非整数），加上磷之后正确得到 12。硫仍然不计入（噻吩 C₄H₄S 用含硫公式会错算成 3.5）。已补测试锁定这个行为。
- **移除了谱图找峰工具**：最初实现过 `peak` 子命令（基于中位数基线与 1.4826×MAD 稳健噪声估计的找峰、半高宽插值），但它属于分析化学方向，与本项目聚焦的有机/物化/药学/材料不一致，因此删除。想找回它：

  ```bash
  git show c5d0484:chem_toolbox/core.py        # 查看当时的实现
  git log --oneline                            # 查看完整提交历史
  ```

  代码在 git 历史里完整保留，需要时随时可以恢复。

---

## License

[MIT](LICENSE)
