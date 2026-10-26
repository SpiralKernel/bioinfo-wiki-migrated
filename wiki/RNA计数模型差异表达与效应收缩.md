# RNA计数模型差异表达与效应收缩

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26官方/原始来源核对；标准库教学演示实际执行；专业验证范围见正文。

## 问题与前置

差异表达不是简单比较两列计数，而是问条件间变化是否超过样本内/样本间的随机变异。学完应能区分归一化、离散度、检验和效应收缩。先读[RNA文库大小因子与组成偏差](https://github.com/SpiralKernel/bioinfo-wiki/wiki/RNA%E6%96%87%E5%BA%93%E5%A4%A7%E5%B0%8F%E5%9B%A0%E5%AD%90%E4%B8%8E%E7%BB%84%E6%88%90%E5%81%8F%E5%B7%AE)、[负二项离散度与均值方差](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%B4%9F%E4%BA%8C%E9%A1%B9%E7%A6%BB%E6%95%A3%E5%BA%A6%E4%B8%8E%E5%9D%87%E5%80%BC%E6%96%B9%E5%B7%AE)和[RNA设计矩阵配对与对比方向](https://github.com/SpiralKernel/bioinfo-wiki/wiki/RNA%E8%AE%BE%E8%AE%A1%E7%9F%A9%E9%98%B5%E9%85%8D%E5%AF%B9%E4%B8%8E%E5%AF%B9%E6%AF%94%E6%96%B9%E5%90%91)。

## 把每个符号对应到数据

对基因$g$、样本$j$，常见负二项模型为 $K_{gj}\sim NB(\mu_{gj},\alpha_g)$，$\mu_{gj}=s_jq_{gj}$，$\log q_{gj}=x_j^T\beta_g$。$K$是观察计数，$s_j$是相对文库尺度，$x_j$是设计矩阵一行，$\beta_g$描述条件及协变量，$\alpha_g$是离散度。其方差为 $\mu+\alpha\mu^2$；均值100、离散度0.1时方差1100，远大于Poisson的100。

DESeq2与edgeR使用相应计数模型和跨基因信息；limma-voom对变换后的计数估计均值—方差权重再拟合线性模型。它们对低计数、离散度和权重的处理不同，不是换个名称的普通TPM t检验。[DESeq2方法与教程](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)、[limma用户指南](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)支持这些方法区别。

## 先手算尺度，再理解模型

两组各两个样本，计数20、40、80、160，已知尺度1、2、1、2；除以尺度得到20、20、80、80，组均值倍数为4，log2FC=2。若只比较某两列40与80，会得到倍数2，混进文库尺度差。本例故意给定尺度，不是在四个数上估计DESeq2 size factor，也没有足够变异结构去检验真实显著性。

归一化不负责移除所有批次和混杂；设计矩阵定义比较与条件化关系。若批次1全是A、批次2全是B，负二项模型也不能从观测中唯一拆出两个效应。转录本估计导入、长度修正和整数计数接口各有约定，不能随意把TPM四舍五入冒充原始计数。

## 为什么要收缩LFC

一个低信息基因的原始LFC可能极端而SE很大。用独立的正态示意：观测效应$\hat\beta\sim N(\beta,s^2)$，先验$\beta\sim N(0,\tau^2)$，后验均值为 $\hat\beta\tau^2/(\tau^2+s^2)$。$\hat\beta=3,s=2,\tau=1$时得到0.6。信息弱时更靠近0，信息强时更靠近原估计。

这不是apeglm算法：真实实现可用不同先验、估计与区间。收缩便于稳定排名和展示，不等于把任意收缩LFC除以原SE产生新P值。检验“非零”和“超过有意义阈值”也不同，需调用支持相应原假设的接口。

## 操作链与检查点

输入为基因×样本计数、唯一且正确对齐的样本表、设计和明确对比。先核对计数来源、样本对应与重复，再定预过滤、拟合、诊断、对比及多重校正，最后按需要收缩效应。输出保留计数基线、LFC、SE/区间、原始P、padj、NA原因和完整检验家族；top50只适合预览。

| 症状 | 候选原因 | 检查 |
|---|---|---|
| 大量NA | 全零、异常、独立过滤等 | 对照baseMean与具体模型状态 |
| 低计数LFC极大 | 信息弱或单样本驱动 | 原计数、SE与收缩前后 |
| 组方向颠倒 | 参考和contrast改变 | 显式写B−A及因子水平 |
| 所有基因同向 | 组成偏差或设计问题 | size factor、样本图及控制依据 |

## 自测与答案

1. 标准化值相同意味着原始计数相同吗？答：不一定，尺度可以不同。
2. 方差1100能用Poisson方差100替代吗？答：会忽略相应过度离散。
3. 收缩后不显眼是否证明基因无效？答：只是估计与信息的综合，需结合区间和目标判断。

[公开airway数据贯通DESeq2](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%85%AC%E5%BC%80airway%E6%95%B0%E6%8D%AE%E8%B4%AF%E9%80%9ADESeq2)已有原生R的DESeq2/apeglm公开计数实测；本轮新增是下方模型算术，edgeR、voom和完整上游RNA流程不因此变成已测。

## 在电脑上核对这一小步

先[下载项目ZIP](https://github.com/SpiralKernel/bioinfo-wiki/archive/refs/heads/main.zip)并解压，进入同时含README.md和examples的项目根目录。不熟悉终端和当前目录时读[R与Python路径文件与项目组织](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E4%B8%8EPython%E8%B7%AF%E5%BE%84%E6%96%87%E4%BB%B6%E4%B8%8E%E9%A1%B9%E7%9B%AE%E7%BB%84%E7%BB%87)。以下在**系统终端**执行，需Python 3.10+，只用标准库，无需pip安装；Windows可将python3替换为已安装的python或py。

```bash
python3 examples/methods_omics_demo.py --case rna
```

输入是脚本内置人工数字，不需下载专业数据。成功会打印JSON，并生成`examples/output/methods-omics/rna/summary.json`；`cases.rna`下的字段对应本页演算。非空输出目录会拒绝覆盖，第二次练习可在命令后加`--output examples/output/methods-omics/rna-practice2`。不要为了重跑删除不明来源的数据。

**检查点：**normalized=[20,20,80,80]、normalized_log2_fold_change=2、NB_variance=1100、normal_prior_toy.posterior_mean=0.6。

阅读源码时先找到同名函数，逐一对照输入、运算与输出；要修改数字，先复制整个解压项目并运行副本，再手算新的预期。脚本含若干固定结果断言，修改教学条件时也需重新核对相应断言；断言失败不自动说明专业方法失效。本轮已实际运行全部22个标准库模式，Python 3.14.8；本页范围为已知尺度、负二项方差和正态先验演算；不重新运行DESeq2/apeglm，未运行edgeR/voom。演示没有为真实样本产生专业P值或性能结论。


## 来源与下一步

本轮查阅日期为2026-10-26；正文原始来源支持对应技术事实，手算与代码是本项目人工教学设计。方法页解释模型，安装及大数据运行交给相应工具/工作流页。教学组织参考[GTN](https://training.galaxyproject.org/training-material/)的步骤与检查点、[Bioinformatics Workbook](https://bioinformaticsworkbook.org/)的问题驱动练习；未复制其教程正文。

保留的相关入口：[RNA-seq从计数矩阵到差异表达](https://github.com/SpiralKernel/bioinfo-wiki/wiki/RNA-seq%E4%BB%8E%E8%AE%A1%E6%95%B0%E7%9F%A9%E9%98%B5%E5%88%B0%E5%B7%AE%E5%BC%82%E8%A1%A8%E8%BE%BE)、[转录组差异与qPCR验证设计](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%BD%AC%E5%BD%95%E7%BB%84%E5%B7%AE%E5%BC%82%E4%B8%8EqPCR%E9%AA%8C%E8%AF%81%E8%AE%BE%E8%AE%A1)、[公开airway数据贯通DESeq2](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%85%AC%E5%BC%80airway%E6%95%B0%E6%8D%AE%E8%B4%AF%E9%80%9ADESeq2)、[DESeq2官方教程](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)、[voom原始论文](https://doi.org/10.1186/gb-2014-15-2-r29)、[edgeR原始论文](https://doi.org/10.1093/bioinformatics/btp616)、[广义线性模型与逻辑回归](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%B9%BF%E4%B9%89%E7%BA%BF%E6%80%A7%E6%A8%A1%E5%9E%8B%E4%B8%8E%E9%80%BB%E8%BE%91%E5%9B%9E%E5%BD%92)、[线性模型与混合模型](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%BA%BF%E6%80%A7%E6%A8%A1%E5%9E%8B%E4%B8%8E%E6%B7%B7%E5%90%88%E6%A8%A1%E5%9E%8B)、[概率条件概率与贝叶斯定理](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%A6%82%E7%8E%87%E6%9D%A1%E4%BB%B6%E6%A6%82%E7%8E%87%E4%B8%8E%E8%B4%9D%E5%8F%B6%E6%96%AF%E5%AE%9A%E7%90%86)、[负二项离散度与均值方差](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%B4%9F%E4%BA%8C%E9%A1%B9%E7%A6%BB%E6%95%A3%E5%BA%A6%E4%B8%8E%E5%9D%87%E5%80%BC%E6%96%B9%E5%B7%AE)、[RNA文库大小因子与组成偏差](https://github.com/SpiralKernel/bioinfo-wiki/wiki/RNA%E6%96%87%E5%BA%93%E5%A4%A7%E5%B0%8F%E5%9B%A0%E5%AD%90%E4%B8%8E%E7%BB%84%E6%88%90%E5%81%8F%E5%B7%AE)、[RNA效应收缩与阈值检验](https://github.com/SpiralKernel/bioinfo-wiki/wiki/RNA%E6%95%88%E5%BA%94%E6%94%B6%E7%BC%A9%E4%B8%8E%E9%98%88%E5%80%BC%E6%A3%80%E9%AA%8C)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
