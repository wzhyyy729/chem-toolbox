# Chem Toolbox · 化学计算工具箱

[English](README.md) | **简体中文**

[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Dependencies](https://img.shields.io/badge/dependencies-none-3FB950)](#安装)
[![Tests](https://img.shields.io/badge/tests-78%20passed-3FB950)](#自检测试)
[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

化学专业学生日常要做的四类计算，收进一个命令行工具里：**分子式分析、溶液配制与 pH、实验数据统计与拟合、谱图找峰**。

只用 Python 标准库，**零第三方依赖** —— 拷到任何装了 Python 的电脑上都能直接跑，不用 pip install。

> 这个项目由化学专业大二学生发起，代码在 AI 协作下产出，但**每一个公式都经过手工核对，并有 78 项自动化断言覆盖**。写得不算花哨，但算出来的数字你可以拿去交实验报告。

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
| `peak` | 谱图找峰：峰位、峰高、信噪比、半高宽（FWHM） |

**特别支持**：嵌套括号化学式（`K3[Fe(CN)6]`）、水合物（`CuSO4.5H2O`）、带注释行的仪器导出数据（自动跳过 `#` 开头和空行）、自动识别逗号/制表符/空格分隔。

---

## 为什么做这个

做实验报告时反复遇到同一批计算：配溶液要算称多少克、平行滴定要算标准差、标准曲线要算 R² 和未知样浓度、紫外谱图要读峰位。用计算器一个个敲，既慢又容易抄错。

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

不饱和度 DBE  : 不适用 —— 分子含 Cu、S，不饱和度公式仅适用于 C/H/N/O/卤素组成的有机物
```

嵌套括号也认：

```console
$ python -m chem_toolbox formula "K3[Fe(CN)6]"
摩尔质量      : 329.2470 g/mol
不饱和度 DBE  : 不适用 —— 分子含 Fe、K，不饱和度公式仅适用于 C/H/N/O/卤素组成的有机物
```

有机物的不饱和度会正常算：

```console
$ python -m chem_toolbox formula "C6H6"
不饱和度 DBE  : 4.0000        # 苯：3 个双键 + 1 个环
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

### 谱图找峰

```console
$ python -m chem_toolbox peak --file examples/uvvis_spectrum.csv --x-col 0 --y-col 1
噪声水平估计    : 0.0199
判定阈值        : 3.0000 × 噪声 = 0.0598

共找到 2 个峰
  序号             峰位            峰高          凸出高度       信噪比         半高宽
  1         455.000        0.4420        0.3080      15.5      37.635
  2         605.000        0.2850        0.1510       7.6      44.359
```

### 脚本调用

除了命令行，也可以当普通 Python 库用：

```python
from chem_toolbox import analyse_formula, ph_weak_acid, linear_fit

report = analyse_formula("C6H12O6")
print(report.molar_mass)            # 180.156
print(report.dbe)                   # None（无碳以外... 实际为 0.0，葡萄糖无环无双键）

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

$$\mathrm{DBE} = \frac{2C + 2 + N - H - X}{2}$$

其中 C、N、H 是碳/氮/氢原子数，X 是卤素原子总数；氧和硫不计入。

- 适用条件：只适用于 **C、H、N、O、卤素**组成的有机物。含金属的式子（如 `K3[Fe(CN)6]`）公式不成立，程序会拒绝计算而不是给个错数。
- 非整数结果说明分子式写错了；负值说明氢或卤素过多，组成不合法。

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

### 谱图找峰

1. **基线**用中位数估计，而不是均值 —— 峰本身会把均值抬高，中位数对少数强峰不敏感。
2. **噪声**用相邻点差分的稳健估计：$\hat{\sigma} \approx 1.4826 \times \mathrm{median}(|\Delta y|) / \sqrt{2}$。系数 1.4826 把中位数绝对偏差折算为标准差，除以 $\sqrt{2}$ 是因为相邻两点相减会把方差放大 2 倍。
3. **判定**：局部极大且凸出高度 ≥ 3σ 才算峰（可用 `--sigma` 调整）。
4. **半高宽 FWHM**：在 $y = \text{基线} + \text{凸出高度}/2$ 处向两侧线性插值求交点。对高斯峰，理论上 $\mathrm{FWHM} = 2\sqrt{2\ln 2}\,\sigma \approx 2.3548\sigma$。

---

## 自检测试

```console
$ python -m chem_toolbox.tests
==============================================================
chem_toolbox 自检测试
==============================================================
  [通过] 分子式解析
  [通过] 溶液与 pH
  [通过] 统计与拟合
  [通过] 谱图找峰
--------------------------------------------------------------
断言通过 78 项，失败 0 项
全部通过
```

测试特点：**每个期望值都能用教科书或手算核对**，不是"跑起来不报错就算过"。例如：

- 葡萄糖摩尔质量 180.156 g/mol、胆矾 249.68 g/mol、`K3[Fe(CN)6]` 329.24 g/mol
- 苯的 DBE = 4、环己烷 = 1、乙烷 = 0；含金属时拒绝计算
- 0.1 mol/L 醋酸（pKₐ 4.76）pH = 2.8813 ± 0.0005；0.1 mol/L 氨水（K_b 1.8×10⁻⁵）pH = 11.1248
- 人造高斯峰的 FWHM 应落在 2.3548σ ± 1.5 内
- 非法输入（空分子式、小写元素符号、括号不匹配、负浓度、斜率全相同）必须抛错，而不是悄悄算错

不需要 pytest，直接 `python -m chem_toolbox.tests` 就能跑。

---

## 已知限制

写清楚边界，比假装什么都能算更负责：

- **原子质量表只收了约 40 种常见元素**（H 到 Pb 的常用项）。需要别的元素，请在 `chem_toolbox/core.py` 的 `ATOMIC_MASS` 里补一行即可。
- **pH 只支持一元酸碱**。多元酸（如 H₃PO₄、H₂CO₃）需要逐级电离常数联立求解，暂未实现。
- **不计算活度系数**，用浓度代替活度。高离子强度（如 ≥ 0.1 mol/L 的强电解质）下结果会有偏差。
- **找峰不做峰形拟合**，只做局部极大 + 阈值判定，且要求 x 轴严格递增（会自动反转递减数据，但不会插值）。
- **不画图**。这是刻意取舍：绘图会引入 matplotlib 依赖，破坏"零依赖"这个最大的便利性。需要作图请把 CSV 丢给 Excel 或 Origin。

---

## 目录结构

```
chem-toolbox/
├── chem_toolbox/
│   ├── __init__.py     对外接口（方便 from chem_toolbox import ...）
│   ├── __main__.py     命令行界面：参数解析与结果打印
│   ├── core.py         纯计算逻辑，无打印、无文件读写，便于单独测试
│   └── tests.py        78 项自检断言，无需 pytest
├── examples/
│   ├── standard_curve.csv    分光光度法标准曲线示例
│   ├── titration.csv         平行滴定数据示例
│   └── uvvis_spectrum.csv    紫外-可见吸收光谱示例（双峰）
├── LICENSE
└── README.md
```

`core.py` 与 `__main__.py` 刻意分开：算法不含任何输入输出，所以既能被命令行调用，也能被别的脚本 import，还能被测试直接断言。

---

## License

[MIT](LICENSE)
