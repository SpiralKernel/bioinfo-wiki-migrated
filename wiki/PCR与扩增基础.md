# PCR与扩增基础

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26官方/原始来源核对；标准库教学演示实际执行；专业验证范围见正文。

## 学完应能回答什么

PCR把已有模板复制成更多可检测拷贝，不能增加原始样本数。本页目标是算清增长与偏差、解释为什么扩增产量不能直接当起始分子量。先读[从DNA到RNA和蛋白质](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E4%BB%8EDNA%E5%88%B0RNA%E5%92%8C%E8%9B%8B%E7%99%BD%E8%B4%A8)、[引物扩增子与特异性](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%BC%95%E7%89%A9%E6%89%A9%E5%A2%9E%E5%AD%90%E4%B8%8E%E7%89%B9%E5%BC%82%E6%80%A7)；RNA先逆转录与DNA直接扩增是不同输入。本页不提供针对具体物种的实验温度、引物序列或优化配方。

## 从一次复制理解指数增长

变性分开双链，引物界定可延伸的起点，聚合酶沿模板合成新链，反复循环使目标积累。这是[NHGRI原理说明](https://www.genome.gov/about-genomics/fact-sheets/Polymerase-Chain-Reaction-Fact-Sheet)中的基本复制过程。特异性来自模板、引物与反应条件共同作用，看到一条条带不等于已经验证序列身份。

用简化模型 $N_c=N_0(1+E)^c$。$N_0$为起始可扩增目标拷贝，$c$为循环数，$E$为每轮新增拷贝相对现有拷贝的比例；$E=1$才是每轮翻倍。模型假设效率恒定、没有竞争或平台期。真实体系后期资源限制与抑制使这些假设失效，不能把全程强套成一条指数曲线。

## 手算增长，也手算偏差

10份模板、10轮、理想效率1，结果是 $10×2^{10}=10240$份。若效率0.8，则是 $10×1.8^{10}≈3570.47$，不是8192。两种目标初始一样多，效率分别0.9和0.8，20轮后的比值为 $(1.9/1.8)^{20}≈2.95$；末端接近三倍差异可以纯由扩增效率产生。这说明跨目标或等位定量需要测量效率和验证假设。

拷贝的谱系也影响证据。假设一个原始分子早期出现复制错误，它的许多后代会共同带错；它们不是许多个独立的错误事件。UMI若用于识别原始分子，必须在相应扩增前引入并考虑碰撞、测序错误及分组规则，不能只在事后给reads任意编号就消除偏差。

## 怎样把原理用于读实验记录

先明确PCR出现在靶向选择、文库扩增还是平台信号生成阶段，再记录原始输入、对照、扩增记录和测序布局。无模板对照帮助发现部分污染；阳性对照检查相应链路；阴性结果仍可能来自抑制、引物不匹配或输入不足。对照的作用不同，不能互相替代。

| 现象 | 可能原因 | 先核对什么 |
|---|---|---|
| 终点信号相似，起始量却不同 | 平台期抹平差异 | 是否使用有效定量区间和效率信息 |
| 两个等位比例明显偏移 | 引物结合差异、低分子抽样 | 原始分子量、等位覆盖及独立验证 |
| 阴性对照也有产物 | 污染或非特异扩增 | 对照身份、产物序列及处理批次 |
| 测序深度高但独立分子少 | 扩增拷贝堆积 | 去重定义、UMI与文库复杂度 |

增加循环通常增加信号，也可能放大偏差；不能只按产量选“最佳”条件。终点PCR、实时qPCR与数字PCR的观测量及定量模型不同，真实报告可按[MIQE 2.0](https://pubmed.ncbi.nlm.nih.gov/40272429/)核对实验与分析信息，接[RT-PCR实时PCR与数字PCR的区别](https://github.com/SpiralKernel/bioinfo-wiki/wiki/RT-PCR%E5%AE%9E%E6%97%B6PCR%E4%B8%8E%E6%95%B0%E5%AD%97PCR%E7%9A%84%E5%8C%BA%E5%88%AB)。

## 自测与答案

1. 一份原始分子扩成一万份，生物学重复数变成一万吗？答：没有，仍来自原来分子与样本的谱系。
2. 增长慢一半是否一定说明起始模板少一半？答：不一定，效率、抑制和目标身份也会影响结果。
3. PCR-free是否表示平台完全没有信号放大？答：需说明阶段，文库处理与簇/DNB等模板呈现不是同一概念。

## 在电脑上核对这一小步

先[下载项目ZIP](https://github.com/SpiralKernel/bioinfo-wiki/archive/refs/heads/main.zip)并解压，进入同时含README.md和examples的项目根目录。不熟悉终端和当前目录时读[R与Python路径文件与项目组织](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E4%B8%8EPython%E8%B7%AF%E5%BE%84%E6%96%87%E4%BB%B6%E4%B8%8E%E9%A1%B9%E7%9B%AE%E7%BB%84%E7%BB%87)。以下在**系统终端**执行，需Python 3.10+，只用标准库，无需pip安装；Windows可将python3替换为已安装的python或py。

```bash
python3 examples/methods_omics_demo.py --case pcr
```

输入是脚本内置人工数字，不需下载专业数据。成功会打印JSON，并生成`examples/output/methods-omics/pcr/summary.json`；`cases.pcr`下的字段对应本页演算。非空输出目录会拒绝覆盖，第二次练习可在命令后加`--output examples/output/methods-omics/pcr-practice2`。不要为了重跑删除不明来源的数据。

**检查点：**ideal=10240，efficiency_0_8约3570.467，relative_bias_after_20_cycles约2.94864。

阅读源码时先找到同名函数，逐一对照输入、运算与输出；要修改数字，先复制整个解压项目并运行副本，再手算新的预期。脚本含若干固定结果断言，修改教学条件时也需重新核对相应断言；断言失败不自动说明专业方法失效。本轮已实际运行全部22个标准库模式，Python 3.14.8；本页范围为指数增长和效率差异算术；未执行PCR湿实验或验证定量试剂。演示没有为真实样本产生专业P值或性能结论。


## 来源与下一步

本轮查阅日期为2026-10-26；正文原始来源支持对应技术事实，手算与代码是本项目人工教学设计。方法页解释模型，安装及大数据运行交给相应工具/工作流页。教学组织参考[GTN](https://training.galaxyproject.org/training-material/)的步骤与检查点、[Bioinformatics Workbook](https://bioinformaticsworkbook.org/)的问题驱动练习；未复制其教程正文。

保留的相关入口：[RT-PCR实时PCR与数字PCR的区别](https://github.com/SpiralKernel/bioinfo-wiki/wiki/RT-PCR%E5%AE%9E%E6%97%B6PCR%E4%B8%8E%E6%95%B0%E5%AD%97PCR%E7%9A%84%E5%8C%BA%E5%88%AB)、[引物扩增子与特异性](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%BC%95%E7%89%A9%E6%89%A9%E5%A2%9E%E5%AD%90%E4%B8%8E%E7%89%B9%E5%BC%82%E6%80%A7)、[扩增子测序与靶向测序设计](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%89%A9%E5%A2%9E%E5%AD%90%E6%B5%8B%E5%BA%8F%E4%B8%8E%E9%9D%B6%E5%90%91%E6%B5%8B%E5%BA%8F%E8%AE%BE%E8%AE%A1)、[NHGRI PCR](https://www.genome.gov/about-genomics/fact-sheets/Polymerase-Chain-Reaction-Fact-Sheet)、[MIQE 2.0](https://pubmed.ncbi.nlm.nih.gov/40272429/)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
