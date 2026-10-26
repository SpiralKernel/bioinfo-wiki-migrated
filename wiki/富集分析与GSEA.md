# 富集分析与GSEA

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26官方/原始来源核对；标准库教学演示实际执行；专业验证范围见正文。

## 列表问题和排序问题先分开

过度表示分析ORA问“选中列表中某注释是否过多”；GSEA类排序分析问“一个基因集是否集中在完整排序的一端”。学完应能写出背景分母、手算一个ORA概率，并解释富集分数与校正后证据的区别。前置读[基因标识符与同源基因](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%9F%BA%E5%9B%A0%E6%A0%87%E8%AF%86%E7%AC%A6%E4%B8%8E%E5%90%8C%E6%BA%90%E5%9F%BA%E5%9B%A0)、[P值效应量与多重检验](https://github.com/SpiralKernel/bioinfo-wiki/wiki/P%E5%80%BC%E6%95%88%E5%BA%94%E9%87%8F%E4%B8%8E%E5%A4%9A%E9%87%8D%E6%A3%80%E9%AA%8C)和[GO术语注释与证据](https://github.com/SpiralKernel/bioinfo-wiki/wiki/GO%E6%9C%AF%E8%AF%AD%E6%B3%A8%E9%87%8A%E4%B8%8E%E8%AF%81%E6%8D%AE)。

## ORA的四个数

设背景有$N$个实际有机会测量、检验并映射的基因，某基因集在背景中有$K$个；选中列表有$n$个，其中$k$个属于该集。零假设下的超几何概率为

$$
P(X\ge k)=\sum_{i=k}^{\min(K,n)}\frac{{K\choose i}{N-K\choose n-i}}{{N\choose n}}.
$$

这是在固定列表大小下从背景抽取的模型；并不意味着基因表达或选择过程真的随机且完全可交换。长基因/高表达更易入选等选择偏差仍需处理。[GO官方富集说明](https://www.geneontology.org/docs/go-enrichment-analysis/)强调合适的参考列表与映射，不能机械使用全物种已注释基因。

手算$N=20,K=5,n=4,k=3$：恰3个命中有$\binom53\binom{15}1=150$种，恰4个有5种，所有抽法$\binom{20}4=4845$，故上尾P=155/4845≈0.03199。这是单个集合原始P，检验许多集合后必须按声明的家族校正；不能把它直接称FDR或“通路已激活”。

## GSEA的运行和分数

按有方向的统计量给全部合格基因排序；对某集合从头走到尾，遇集合内基因加分，集合外减分。无权toy中，8个基因、命中排名1/2/6，每次命中加1/3，未命中减1/5，运行和在第2位达到2/3，最后回到0。真实加权GSEA可按排名统计量绝对值赋权，ES取相应极端偏离，NES再按零分布等规则归一化；toy的2/3不能直接当NES。

[GSEA原始论文](https://doi.org/10.1073/pnas.0506580102)及[用户指南](https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/)解释排序、置换和分数。表型置换、基因集/预排序零分布对应不同假设，配对或重复测量不能随意交换所有样本标签。只保留显著基因再排序，已经删掉零分布所需的信息，不再是同一个完整排序问题。

## 实际输入与输出

先冻结基因ID、版本、映射规则和重复处理；ORA保存背景与选中阈值，GSEA保存完整可用排序、正方向、统计量及集合大小限制。不能为每个条目单独挑最有利转录本或把一个基因一对多映射重复计数。输出保留集合ID/版本、大小、命中/leading edge、原始与调整证据、ES/NES及失败原因。

RNA差异表中的padj NA不一定都同一原因，应按可检验和可映射目标定义背景；[RNA差异结果报告与富集背景](https://github.com/SpiralKernel/bioinfo-wiki/wiki/RNA%E5%B7%AE%E5%BC%82%E7%BB%93%E6%9E%9C%E6%8A%A5%E5%91%8A%E4%B8%8E%E5%AF%8C%E9%9B%86%E8%83%8C%E6%99%AF)与[公开airway数据贯通DESeq2](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%85%AC%E5%BC%80airway%E6%95%B0%E6%8D%AE%E8%B4%AF%E9%80%9ADESeq2)提供完整结果入口。top50是预览，不是默认背景，也不是完整排名。

| 症状 | 可能原因 | 检查 |
|---|---|---|
| 换ID后结论骤变 | 映射损失/一对多 | 各阶段唯一基因数和丢失表 |
| 上百GO项都显著 | 层级重叠、同一批驱动基因 | 命中交集，不按条目数算机制数 |
| 同一通路正负方向混乱 | 排名方向或混合集合 | 指标符号、leading edge和基因层效应 |
| P很小却只有少数基因 | 集合小或背景不匹配 | 四个数、大小过滤及校正家族 |

## 自测与答案

1. 背景从20改200而K不动会影响ORA吗？答：会改变抽样比例问题，不能任意换分母。
2. ES=2/3是否就是FDR=1/3？答：不是，分数与零分布/多重校正不同。
3. 富集显著能确认因果调控吗？答：只支持集合关联，仍需独立机制证据。

## 在电脑上核对这一小步

先[下载项目ZIP](https://github.com/SpiralKernel/bioinfo-wiki/archive/refs/heads/main.zip)并解压，进入同时含README.md和examples的项目根目录。不熟悉终端和当前目录时读[R与Python路径文件与项目组织](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E4%B8%8EPython%E8%B7%AF%E5%BE%84%E6%96%87%E4%BB%B6%E4%B8%8E%E9%A1%B9%E7%9B%AE%E7%BB%84%E7%BB%87)。以下在**系统终端**执行，需Python 3.10+，只用标准库，无需pip安装；Windows可将python3替换为已安装的python或py。

```bash
python3 examples/methods_omics_demo.py --case enrichment
```

输入是脚本内置人工数字，不需下载专业数据。成功会打印JSON，并生成`examples/output/methods-omics/enrichment/summary.json`；`cases.enrichment`下的字段对应本页演算。非空输出目录会拒绝覆盖，第二次练习可在命令后加`--output examples/output/methods-omics/enrichment-practice2`。不要为了重跑删除不明来源的数据。

**检查点：**ORA.upper_tail_p约0.031991744；positive_ES=2/3；running score最终约0（浮点误差可到1e-16）。

阅读源码时先找到同名函数，逐一对照输入、运算与输出；要修改数字，先复制整个解压项目并运行副本，再手算新的预期。脚本含若干固定结果断言，修改教学条件时也需重新核对相应断言；断言失败不自动说明专业方法失效。本轮已实际运行全部22个标准库模式，Python 3.14.8；本页范围为超几何上尾和无权运行和；未运行专业GSEA的置换、NES或FDR。演示没有为真实样本产生专业P值或性能结论。


## 来源与下一步

本轮查阅日期为2026-10-26；正文原始来源支持对应技术事实，手算与代码是本项目人工教学设计。方法页解释模型，安装及大数据运行交给相应工具/工作流页。教学组织参考[GTN](https://training.galaxyproject.org/training-material/)的步骤与检查点、[Bioinformatics Workbook](https://bioinformaticsworkbook.org/)的问题驱动练习；未复制其教程正文。

保留的相关入口：[基因标识符与同源基因](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%9F%BA%E5%9B%A0%E6%A0%87%E8%AF%86%E7%AC%A6%E4%B8%8E%E5%90%8C%E6%BA%90%E5%9F%BA%E5%9B%A0)、[GO 官方富集说明](https://www.geneontology.org/docs/go-enrichment-analysis/)、[GSEA 原始论文](https://doi.org/10.1073/pnas.0506580102)、[GSEA 用户指南](https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/)、[GO术语注释与证据](https://github.com/SpiralKernel/bioinfo-wiki/wiki/GO%E6%9C%AF%E8%AF%AD%E6%B3%A8%E9%87%8A%E4%B8%8E%E8%AF%81%E6%8D%AE)、[RNA差异结果报告与富集背景](https://github.com/SpiralKernel/bioinfo-wiki/wiki/RNA%E5%B7%AE%E5%BC%82%E7%BB%93%E6%9E%9C%E6%8A%A5%E5%91%8A%E4%B8%8E%E5%AF%8C%E9%9B%86%E8%83%8C%E6%99%AF)、[公开airway数据贯通DESeq2](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%85%AC%E5%BC%80airway%E6%95%B0%E6%8D%AE%E8%B4%AF%E9%80%9ADESeq2)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
