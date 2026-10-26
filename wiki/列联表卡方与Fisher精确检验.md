# 列联表卡方与Fisher精确检验

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26核对原始与官方资料；本页Python微型计算已执行；人工数据不代表专业软件或真实研究性能。

列联表用独立对象的分类计数检查两类变量是否有关联，例如两组样本中阳性/阴性比例是否相同。Pearson卡方依靠近似分布，Fisher精确检验在固定边际条件下计算表的概率。读完应能检查输入是不是计数、手算期望频数，并说明OR与风险比的不同。

前置为[概率条件概率与贝叶斯定理](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%A6%82%E7%8E%87%E6%9D%A1%E4%BB%B6%E6%A6%82%E7%8E%87%E4%B8%8E%E8%B4%9D%E5%8F%B6%E6%96%AF%E5%AE%9A%E7%90%86)和[假设检验与双侧单侧选择](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%81%87%E8%AE%BE%E6%A3%80%E9%AA%8C%E4%B8%8E%E5%8F%8C%E4%BE%A7%E5%8D%95%E4%BE%A7%E9%80%89%E6%8B%A9)。若同一对象前后两次分类配对，应考虑McNemar等方法，而不是把四格表当独立两组。

## 卡方从哪里来

对r×c表，H0为行列独立。每格期望计数$E_{ij}=R_iC_j/N$，Ri/Cj为行列总数，N为总数。Pearson统计量$X^2=\sum_{ij}(O_{ij}-E_{ij})^2/E_{ij}$；在适当大样本条件下与自由度(r−1)(c−1)卡方分布比较。近似是否可靠取决于期望计数结构等，不能只看总N。[SciPy卡方列联表文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.chi2_contingency.html)说明期望计数与校正选项。

观察到的比例或百分数不能直接冒充原始计数。相同50%可以来自1/2或500/1000，不确定性不同。单个对象多次出现也会破坏独立计数假设。

## 用四格表完整核对

组A阳性8、阴性2；组B阳性1、阴性9。行总数都是10，阳性总9、阴性总11，故每行期望为4.5与5.5。未经连续性校正的X²≈9.89899，自由度1。OR=(8×9)/(2×1)=36；风险比(8/10)/(1/10)=8；风险差0.7。它们的单位和含义不相同，OR=36不等于“概率36倍”。

> 本页例子在项目根目录的Python中运行；标准库例子无需安装包，含NumPy/SciPy的例子先按[R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83)建立项目环境。Windows命令使用`python`，macOS/Linux可使用`python3`；代码块可复制进该环境的Python交互会话，或保存为项目根目录下的临时脚本再用该环境解释器运行；不要保存到examples目录或任意外部目录后假定项目导入路径仍相同。除特别说明外，只向屏幕输出，不写文件。

```python
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
T = np.array([[8,2],[1,9]])
chi = chi2_contingency(T,correction=False)
fisher = fisher_exact(T,alternative='two-sided')
print('expected', chi.expected_freq.tolist())
print('chi2_df_p', round(float(chi.statistic),6), chi.dof, round(float(chi.pvalue),6))
print('OR_Fisher_p', float(fisher.statistic), round(float(fisher.pvalue),6))
print('risk_ratio', (8/10)/(1/10))
assert np.allclose(chi.expected_freq,[[4.5,5.5],[4.5,5.5]])
assert fisher.statistic == 36
```

本例输出卡方P约0.001654、Fisher双侧P约0.005477，可连同期望频数和OR=36核对。此例特意给出`correction=False`便于与手算对应，不表示所有数据都应关闭校正。两种P不必相同，因为使用的参考计算不同。真正报告时选定方法并说明，不是比较后挑小者。

## Fisher检验“精确”在哪里

2×2表在给定行列总数、H0优势比为1时，左上角计数服从超几何分布。枚举所有边际相同的表，按预设单双侧规则累加概率，即可得到P。常见双侧定义累加概率不高于观测表的表，通常不简单等于较小单侧P的两倍；见[Fisher官方说明](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html)。

精确指条件参考分布的计算，不保证样本独立、测量无偏或研究因果成立。零格使普通样本OR为0或无穷，仍可能计算检验P；区间、校正估计及所用OR定义要另作说明，不应随意加0.5后当成原始表。

## 富集与其他场景的边界

基因集过度代表分析常形成“入选/未入选×属于/不属于基因集”表，但背景必须是可被检出的全集，不能默认整基因组。基因依赖、长度和注释偏差也影响解释，见[富集分析与GSEA](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%AF%8C%E9%9B%86%E5%88%86%E6%9E%90%E4%B8%8EGSEA)。多位点、多通路检验还需多重校正。

如需调整年龄、批次或连续协变量，四格表不能完成，需逻辑回归等模型。若一个对象可以属于多个类别，行列可能不是互斥计数，应重新定义统计单位和问题。

## 排错与自测

| 症状 | 优先查什么 |
|---|---|
| 期望计数为0 | 是否有全零行/列，类别是否无观测 |
| 与其他软件P不同 | 连续性校正、单双侧、精确P定义 |
| OR方向相反 | 行列顺序，阳性/阴性编码 |
| P很小但比例差微小 | 样本量与效应尺度，不能只报显著 |

**问：把所有格计数乘100会怎样？** 答：OR和比例不变，但独立样本量如果真实增加，统计精度会改变；复制同一批对象则不是真的增加信息。**问：四格表支持干预因果吗？** 答：单独不支持，还需随机化、混杂控制与设计证据。


## 接着读什么

[概率条件概率与贝叶斯定理](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%A6%82%E7%8E%87%E6%9D%A1%E4%BB%B6%E6%A6%82%E7%8E%87%E4%B8%8E%E8%B4%9D%E5%8F%B6%E6%96%AF%E5%AE%9A%E7%90%86) → [假设检验与双侧单侧选择](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%81%87%E8%AE%BE%E6%A3%80%E9%AA%8C%E4%B8%8E%E5%8F%8C%E4%BE%A7%E5%8D%95%E4%BE%A7%E9%80%89%E6%8B%A9) → [富集分析与GSEA](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%AF%8C%E9%9B%86%E5%88%86%E6%9E%90%E4%B8%8EGSEA) → [R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83) → [合成数据贯通统计推断](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%90%88%E6%88%90%E6%95%B0%E6%8D%AE%E8%B4%AF%E9%80%9A%E7%BB%9F%E8%AE%A1%E6%8E%A8%E6%96%AD) → [方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。学习顺序与练习入口见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。


## 运行本页的完整小例

先下载项目并进入含README与examples的根目录，环境准备见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。用项目环境的Python运行以下命令（macOS/Linux路径示例；Windows用`.venv\Scripts\python.exe`）：

```bash
.venv/bin/python examples/methods_foundations_demo.py --case contingency --output examples/output/methods-contingency
```

输出summary.json保存本页计算与检查结果。目录非空时拒绝覆盖，再练习请更换`--output`目录；页面内短代码只打印结果。标准库case也可用系统Python直接运行。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
