# FarmCPU与BLINK原理

> 更新：2026-10-26；状态：可用；验证范围：官方/原始来源核对；标准库教学算例已实跑；专业流程范围见正文。

一条标记的单独检验可能同时捡到自身效应、附近因果位点及群体背景。把所有标记一起放进固定效应模型又会在m远大于n时欠定。FarmCPU与BLINK试图选择少量背景代理位点，反复在“背景选择”和“逐标记检验”之间更新。目标是理解这种自适应条件检验，而不是把多位点方法等同于因果发现。

前置为[线性模型与混合模型](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%BA%BF%E6%80%A7%E6%A8%A1%E5%9E%8B%E4%B8%8E%E6%B7%B7%E5%90%88%E6%A8%A1%E5%9E%8B)、[连锁不平衡与LD剪枝](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%BF%9E%E9%94%81%E4%B8%8D%E5%B9%B3%E8%A1%A1%E4%B8%8ELD%E5%89%AA%E6%9E%9D)及GWAS样本对齐。输入通常包括n个样本的表型y、协变量C、n×m剂量矩阵M和标记染色体/位置。y、C、M必须逐ID对应；同一模型内的样本剔除与等位计数方向也要一致。

## 先理解条件检验

把本轮选中的背景标记放入矩阵Q，检测第j点时可用简化固定模型 $y=C\gamma+Qb+x_j\beta_j+e$。γ是协变量系数，b是背景标记系数，βj是待检验点的条件效应。若xj与某个Q列几乎相同，两者效应难以分开；若漏掉重要背景，xj又可能承接其信号。

这说明pseudo-QTN是**被选作背景的代理**，不是已经确认的因果核苷酸。选择来自当前数据和迭代规则，不能把Q加入模型以后得到的小P值都视为预先指定检验。被选位点本身的处理、相关位点替换及P值更新须依具体实现，不能用上述一行OLS冒充完整算法。

## FarmCPU的循环

原方法交替进行固定效应检验与随机效应背景评估。固定部分利用当前pseudo-QTN进行逐点扫描；随机部分根据候选背景构建关系并选择配置，原方案涉及按基因组分箱挑代表点以及背景规模/箱宽评估。得到新Q后再扫描，直到候选稳定或达到停止条件。[FarmCPU原始研究](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005767)

读日志时应顺着“本轮候选从哪里来→哪些候选被保留→背景如何更新→为什么停止”追踪。分箱尺度改变了哪些位点竞争同一个代表位置；它不是一个无关紧要的画图参数。将FarmCPU仅称为“更快的MLM”会遗漏候选选择与固定/随机交替这两个关键环节。

## BLINK改变了什么

BLINK使用LD关系减少背景候选的冗余，并以BIC等固定模型比较选择背景规模，替换FarmCPU中的相应随机模型与分箱环节。相近的两个位点若并不处于高LD，不必仅因距离近而挤在一个箱里；反之，远处高度相关位点仍可能冗余。BIC、LD阈值、候选排序与迭代次数共同决定最后背景。[BLINK原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6365300/)

## 手算一次“拟合改善是否值得”

在相同n与响应、独立高斯残差的简化回归下，忽略模型间相同常数，可比较 $BIC=n\log(RSS/n)+k\log n$。RSS为残差平方和，k在此计回归参数；共同残差方差参数的惩罚也作为相同常数略去。截距、协变量等回归参数必须用一致规则计数。

设n=100，候选背景1的RSS=100、k=1，背景2的RSS=80、k=2，BIC分别4.605和−13.104，较低者为背景2。若复杂模型的RSS只降到99，改善可能抵不过惩罚。这个算式只解释复杂度取舍，实际BLINK的候选搜索、检验更新和停止机制不能由一次BIC比较替代。

## 怎样组织一次可审查的分析

先冻结同一组样本、标记、性状变换和协变量，保存GLM/MLM基线；再单独运行FarmCPU或BLINK，保留每轮Q、参数和最终结果。比较共同标记的方向、效应尺度、QQ图、失败和敏感性，而不是按显著数排名。对显著更多的情况，先问是否背景改变、过滤不同或实际假阳性增加。

GAPIT提供两种方法的入口；rMVP包含FarmCPU并不意味着也有同样的BLINK接口。专业运行按[GAPIT多模型GWAS实践](https://github.com/SpiralKernel/bioinfo-wiki/wiki/GAPIT%E5%A4%9A%E6%A8%A1%E5%9E%8BGWAS%E5%AE%9E%E8%B7%B5)与[rMVP关联分析实践](https://github.com/SpiralKernel/bioinfo-wiki/wiki/rMVP%E5%85%B3%E8%81%94%E5%88%86%E6%9E%90%E5%AE%9E%E8%B7%B5)核对版本。真实功效和错误控制要用合适模拟、保留结构的验证或独立队列评估；本页没有新跑两种专业GWAS。

## 动手：运行后先对照什么

在下载项目的根目录执行；只需Python 3标准库，不安装专业统计包：

```bash
python3 examples/methods_genetics_demo.py --case background
```

输出两候选的BIC约为`[4.605170,-13.104015]`。将第二个RSS改成99可先手算复杂度惩罚；不能把较低BIC读成某标记的P值。

该分支只有BIC算术，不生成pseudo-QTN、不实现FarmCPU/BLINK迭代，也不声称重新执行专业软件。

命令只向终端打印JSON。若要保存所有22个小例，可用`python3 examples/methods_genetics_demo.py --output examples/output/methods-genetics-practice1`；目录必须新建或为空，已有结果不会被覆盖。源代码见[标准库方法算例](https://github.com/SpiralKernel/bioinfo-wiki/blob/main/examples/methods_genetics_demo.py)。

## 排错：先定位原因再改参数

| 症状 | 可能原因 | 下一步检查 |
|---|---|---|
| 每轮候选反复变化 | LD冗余、近似等价模型或停止条件 | 保存候选轨迹，调整有依据的LD/背景设置 |
| 加背景后某点消失 | 共线或条件效应改变 | 检查与背景标记的LD、设计秩和效应方向 |
| 显著数突然大增 | 输入、结构控制或自适应选择改变 | 固定共同样本/位点，核对基线与独立证据 |

## 自测与答案

1. **pseudo-QTN与已证实QTN有什么差别？** 前者是当前统计过程选出的背景代理；需要额外定位和实验才能讨论因果。

2. **BIC更小能证明假阳性更少吗？** 不能。它只在给定模型族和数据下折衷拟合与复杂度，并非错误发现率保证。

## 原有实践与延伸入口

[BLINK原始研究](https://doi.org/10.1093/gigascience/giy154)；[GAPIT官方实现](https://github.com/jiabowang/GAPIT)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
