# RNA交互作用与LRT比较

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26官方/原始来源核对；标准库教学演示实际执行；专业验证范围见正文。

## 先把科学问题写成可比较的效应

“处理改变表达吗”和“不同基因型的处理效应是否不同”不是同一检验。学完应能用四个组均值算交互，并解释整体LRT为什么不对应一列LFC。先读[RNA设计矩阵配对与对比方向](https://github.com/SpiralKernel/bioinfo-wiki/wiki/RNA%E8%AE%BE%E8%AE%A1%E7%9F%A9%E9%98%B5%E9%85%8D%E5%AF%B9%E4%B8%8E%E5%AF%B9%E6%AF%94%E6%96%B9%E5%90%91)、[线性模型与混合模型](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%BA%BF%E6%80%A7%E6%A8%A1%E5%9E%8B%E4%B8%8E%E6%B7%B7%E5%90%88%E6%A8%A1%E5%9E%8B)和[似然最大似然与MAP估计](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E4%BC%BC%E7%84%B6%E6%9C%80%E5%A4%A7%E4%BC%BC%E7%84%B6%E4%B8%8EMAP%E4%BC%B0%E8%AE%A1)。

## 四种组合的模型

设$G$为基因型指示量、$T$为处理指示量，二者参考均为0。用对数均值写

$$
\log\mu=\log s+\beta_0+\beta_GG+\beta_TT+\beta_{GT}GT.
$$

$\mu$为期望计数，$s$为文库尺度，$\beta_0$为参考组合截距，$\beta_G,\beta_T$是另一因子在参考水平时的效应；$\beta_{GT}$才是处理效应的差。RNA计数通常还需负二项离散度与样本设计。自然对数系数除以$\log2$才变成log2尺度，不能把模型内log与报告LFC直接混写。[DESeq2复杂设计及LRT说明](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)给出了对应接口与检验含义。

四组标准化均值为G0未处理10、G0处理20、G1未处理20、G1处理80。G0处理倍数2，G1倍数4，处理效应比值为4/2=2，交互log2效应为1。基因型主效应只描述未处理时20/10=2，不能把它解读为所有处理条件下固定不变的效应。

## LRT是“少掉一组参数会怎样”

完整模型含G、T、G:T；简化模型只含G、T。若两者在同一数据与相容设定下拟合，似然比统计量为 $D=2(\ell_{full}-\ell_{reduced})$，$\ell$为最大化后的对数似然。D较大说明删去这组项使拟合损失更多。适当正则条件下可用参数自由度差构造渐近参考分布，但不能只因为D为正就宣称显著。

本页toy将四格数当作独立Poisson计数，简化的无交互期望是行总数×列总数/总数，得到约6.923、23.077、23.077、76.923。饱和模型完全拟合四格，D约2.18174。它没有RNA生物学重复或离散度，故只示范似然损失，不输出RNA差异P值。专业负二项LRT不能用此四格计算替代。

## 从设计到可报告结果

先列出每个组合的样本数、配对及批次，确认模型矩阵可估计。完整模型包含的时间项若有多列，删掉所有时间相关项是在问整体时间效应，不是某一个时间点。报告应写清完整/简化公式、自由度差、检验家族及随后做了哪些具体对比。

“A组P<0.05、B组P>0.05”不能推出交互显著；两组效应可能一样而SE不同。应直接检验效应差及其不确定性。先按LRT筛基因，再挑最显著时间点展示，也涉及额外选择；保存完整对比而非只展示最有利一列。

| 症状 | 原因候选 | 检查 |
|---|---|---|
| 设计矩阵不满秩 | 组合缺失、批次混杂 | 组别交叉表和实际设计列 |
| 主效应方向出乎意料 | 参考水平或尺度改变 | 因子水平、resultsNames及对比 |
| LRT显著但展示LFC小 | 整体检验与单系数不同 | 被删项及全部时间/组间效应 |
| 交互不稳 | 重复少、离群或稀疏组合 | 供体图、覆盖与拟合诊断 |

## 自测与答案

1. 两个处理倍数都等于2时交互是多少？答：log2(2/2)=0。
2. 完整模型和简化模型能用不同样本吗？答：直接LRT要求可比似然，不能随意换数据。
3. 三时间点整体LRT可由一列LFC解释完吗？答：通常不能，需查看时间模式及具体对比。

项目已有[公开airway数据贯通DESeq2](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%85%AC%E5%BC%80airway%E6%95%B0%E6%8D%AE%E8%B4%AF%E9%80%9ADESeq2)只验证加性配对Wald分析，本轮没有运行专业交互/LRT；不能把下方Poisson演示改称真实RNA完整验证。

## 在电脑上核对这一小步

先[下载项目ZIP](https://github.com/SpiralKernel/bioinfo-wiki/archive/refs/heads/main.zip)并解压，进入同时含README.md和examples的项目根目录。不熟悉终端和当前目录时读[R与Python路径文件与项目组织](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E4%B8%8EPython%E8%B7%AF%E5%BE%84%E6%96%87%E4%BB%B6%E4%B8%8E%E9%A1%B9%E7%9B%AE%E7%BB%84%E7%BB%87)。以下在**系统终端**执行，需Python 3.10+，只用标准库，无需pip安装；Windows可将python3替换为已安装的python或py。

```bash
python3 examples/methods_omics_demo.py --case interaction
```

输入是脚本内置人工数字，不需下载专业数据。成功会打印JSON，并生成`examples/output/methods-omics/interaction/summary.json`；`cases.interaction`下的字段对应本页演算。非空输出目录会拒绝覆盖，第二次练习可在命令后加`--output examples/output/methods-omics/interaction-practice2`。不要为了重跑删除不明来源的数据。

**检查点：**interaction_log2=1，ratio_of_ratios=2，poisson_deviance约2.181742；reduced_expected为无交互四格期望。

阅读源码时先找到同名函数，逐一对照输入、运算与输出；要修改数字，先复制整个解压项目并运行副本，再手算新的预期。脚本含若干固定结果断言，修改教学条件时也需重新核对相应断言；断言失败不自动说明专业方法失效。本轮已实际运行全部22个标准库模式，Python 3.14.8；本页范围为四格Poisson与效应差算术；专业RNA交互/LRT未执行。演示没有为真实样本产生专业P值或性能结论。


## 来源与下一步

本轮查阅日期为2026-10-26；正文原始来源支持对应技术事实，手算与代码是本项目人工教学设计。方法页解释模型，安装及大数据运行交给相应工具/工作流页。教学组织参考[GTN](https://training.galaxyproject.org/training-material/)的步骤与检查点、[Bioinformatics Workbook](https://bioinformaticsworkbook.org/)的问题驱动练习；未复制其教程正文。

保留的相关入口：[DESeq2复杂设计指南](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)、[RNA设计矩阵配对与对比方向](https://github.com/SpiralKernel/bioinfo-wiki/wiki/RNA%E8%AE%BE%E8%AE%A1%E7%9F%A9%E9%98%B5%E9%85%8D%E5%AF%B9%E4%B8%8E%E5%AF%B9%E6%AF%94%E6%96%B9%E5%90%91)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
