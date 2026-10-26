# t检验配对设计与Welch检验

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26核对原始与官方资料；本页Python微型计算已执行；人工数据不代表专业软件或真实研究性能。

t检验将均值差除以其标准误，再与适当自由度的t分布比较。先判断独立两组还是同一对象前后配对，比选择哪个函数更重要。读完应能手算均值差与SE，解释Welch为何不要求两组方差相等，并保留原始测量单位。

前置为[假设检验与双侧单侧选择](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%81%87%E8%AE%BE%E6%A3%80%E9%AA%8C%E4%B8%8E%E5%8F%8C%E4%BE%A7%E5%8D%95%E4%BE%A7%E9%80%89%E6%8B%A9)与[置信区间与重采样](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%BD%AE%E4%BF%A1%E5%8C%BA%E9%97%B4%E4%B8%8E%E9%87%8D%E9%87%87%E6%A0%B7)。下文将效应统一为B−A；函数调用顺序需与此对应。

## 独立两组与配对在估计不同的不确定性

独立两组Welch统计量为$t=(\bar B-\bar A)/\sqrt{s_A^2/n_A+s_B^2/n_B}$。s²是样本方差，n是独立观测数。Welch–Satterthwaite自由度为$\nu=(v_A+v_B)^2/[v_A^2/(n_A-1)+v_B^2/(n_B-1)]$，其中vA=sA²/nA、vB=sB²/nB。区间为均值差±t分位数×SE。[SciPy官方说明](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html)中`equal_var=False`选择Welch；默认值不能凭记忆假定。

配对设计先按ID计算di=Bi−Ai，再对d做单样本t检验：$t=\bar d/(s_d/\sqrt n)$，自由度n−1。配对的价值来自保留同一对象的相关结构；把配对拆成独立两组会丢失该信息。反过来，不能为了增加显著性随意将无对应关系的样本配对。

## 一份小数据分别计算

A=[1,2,3,4]，B=[2,4,5,7]。两组均值2.5与4.5，差2；样本方差约1.6667与4.3333，所以Welch SE=√1.5≈1.224745，t≈1.632993，自由度约5.0103。若四对确实对应同一对象，差值[1,2,2,3]，SE≈0.408248，配对t≈4.898979。两者P不同，是独立性模型不同，不是任选更小者。

> 本页例子在项目根目录的Python中运行；标准库例子无需安装包，含NumPy/SciPy的例子先按[R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83)建立项目环境。Windows命令使用`python`，macOS/Linux可使用`python3`；代码块可复制进该环境的Python交互会话，或保存为项目根目录下的临时脚本再用该环境解释器运行；不要保存到examples目录或任意外部目录后假定项目导入路径仍相同。除特别说明外，只向屏幕输出，不写文件。

```python
import numpy as np
from scipy import stats
A = np.array([1.,2.,3.,4.])
B = np.array([2.,4.,5.,7.])
effect = B.mean()-A.mean()
se = np.sqrt(A.var(ddof=1)/len(A)+B.var(ddof=1)/len(B))
welch = stats.ttest_ind(B,A,equal_var=False,alternative='two-sided')
paired = stats.ttest_rel(B,A,alternative='two-sided')
print('effect_SE_t', round(effect,6), round(se,6), round(float(welch.statistic),6))
print('Welch_p_df', round(float(welch.pvalue),6), round(float(welch.df),6))
print('paired_t_p', round(float(paired.statistic),6), round(float(paired.pvalue),6))
assert np.isclose(welch.statistic,effect/se)
assert np.isclose(paired.statistic,np.sqrt(24))
```

数值检查：Welch双侧P约0.163279、df约5.010309；配对双侧P约0.016277。这里同时列出两种计算只为教学比较，真实研究由采样设计决定其中哪种适用。项目[合成数据贯通统计推断](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%90%88%E6%88%90%E6%95%B0%E6%8D%AE%E8%B4%AF%E9%80%9A%E7%BB%9F%E8%AE%A1%E6%8E%A8%E6%96%AD)提供更多独立与配对演示。

## 正态性究竟要求谁

经典小样本t推断依赖各组误差或配对差值的正态假设；不是把所有分组混在一起测“总体正态”。Welch放宽等方差条件，但不自动处理相关样本、严重离群或很小样本的极端偏态。大样本近似也不能靠“n>30”一条规则保证；看分布、独立单位与有效信息。

原始RNA计数、二元结果和比例边界数据通常有更合适的分布模型，见[广义线性模型与逻辑回归](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%B9%BF%E4%B9%89%E7%BA%BF%E6%80%A7%E6%A8%A1%E5%9E%8B%E4%B8%8E%E9%80%BB%E8%BE%91%E5%9B%9E%E5%BD%92)。t检验也不能同时调整年龄、批次等多个协变量；需要回归或混合模型。

## 输入、输出与排错

输入应包含样本ID、数值、组别和必要配对ID。先检查重复ID、单位、缺失、每组有效n、配对缺失是否整对移除。输出保留B−A、SE、df、CI、原始P及多重校正范围。不要用柱状图SE替代原始个体分布。

| 症状 | 可能原因 | 怎样检查 |
|---|---|---|
| 差为正而t为负 | 调用顺序与报告方向不同 | 手算B均值−A均值 |
| 配对后n比预期少 | ID连接缺失或重复 | 检查一对一键与整对有效样本 |
| NaN或无穷t | 所有差值相同、缺失、无变异 | 输出差值与样本方差 |
| 极小P但只有一个供体 | 技术重复被当独立 | 重建生物学单位 |

## 自测与方法比较

**问：Welch是否意味着任何样本量都不怕离群？** 答：不是，它主要处理两组方差不等。**问：同一受试者前后差值为常数，t为无穷就代表结论无限可靠吗？** 答：先检查测量精度、数据构造和方差模型，不能跳过设计。**问：非参数秩检验能直接替代均值差问题吗？** 答：其目标通常不同，参见[秩检验符号检验与非参数方法](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%A7%A9%E6%A3%80%E9%AA%8C%E7%AC%A6%E5%8F%B7%E6%A3%80%E9%AA%8C%E4%B8%8E%E9%9D%9E%E5%8F%82%E6%95%B0%E6%96%B9%E6%B3%95)；应先确定要解释的效应。


## 接着读什么

[假设检验与双侧单侧选择](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%81%87%E8%AE%BE%E6%A3%80%E9%AA%8C%E4%B8%8E%E5%8F%8C%E4%BE%A7%E5%8D%95%E4%BE%A7%E9%80%89%E6%8B%A9) → [置信区间与重采样](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%BD%AE%E4%BF%A1%E5%8C%BA%E9%97%B4%E4%B8%8E%E9%87%8D%E9%87%87%E6%A0%B7) → [合成数据贯通统计推断](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%90%88%E6%88%90%E6%95%B0%E6%8D%AE%E8%B4%AF%E9%80%9A%E7%BB%9F%E8%AE%A1%E6%8E%A8%E6%96%AD) → [广义线性模型与逻辑回归](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%B9%BF%E4%B9%89%E7%BA%BF%E6%80%A7%E6%A8%A1%E5%9E%8B%E4%B8%8E%E9%80%BB%E8%BE%91%E5%9B%9E%E5%BD%92) → [秩检验符号检验与非参数方法](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%A7%A9%E6%A3%80%E9%AA%8C%E7%AC%A6%E5%8F%B7%E6%A3%80%E9%AA%8C%E4%B8%8E%E9%9D%9E%E5%8F%82%E6%95%B0%E6%96%B9%E6%B3%95) → [R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83) → [线性模型与混合模型](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%BA%BF%E6%80%A7%E6%A8%A1%E5%9E%8B%E4%B8%8E%E6%B7%B7%E5%90%88%E6%A8%A1%E5%9E%8B) → [方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。学习顺序与练习入口见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。


## 运行本页的完整小例

先下载项目并进入含README与examples的根目录，环境准备见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。用项目环境的Python运行以下命令（macOS/Linux路径示例；Windows用`.venv\Scripts\python.exe`）：

```bash
.venv/bin/python examples/methods_foundations_demo.py --case t-tests --output examples/output/methods-t-tests
```

输出summary.json保存本页计算与检查结果。目录非空时拒绝覆盖，再练习请更换`--output`目录；页面内短代码只打印结果。标准库case也可用系统Python直接运行。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
