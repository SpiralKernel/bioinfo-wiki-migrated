# PacBioSMRT与HiFi测序

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26官方/原始来源核对；标准库教学演示实际执行；专业验证范围见正文。

## 本页解决的混淆

长read、单次通过、共识read与独立分子不是一回事。学完要能解释HiFi为何同时具有长片段信息与多次观测，也能说明多次pass没有增加生物学样本数。前置读[文库接头索引与条形码](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%87%E5%BA%93%E6%8E%A5%E5%A4%B4%E7%B4%A2%E5%BC%95%E4%B8%8E%E6%9D%A1%E5%BD%A2%E7%A0%81)、[测序深度覆盖度与质量分数](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%B5%8B%E5%BA%8F%E6%B7%B1%E5%BA%A6%E8%A6%86%E7%9B%96%E5%BA%A6%E4%B8%8E%E8%B4%A8%E9%87%8F%E5%88%86%E6%95%B0)。

## 从SMRTbell到共识

SMRT记录聚合酶合成相关信号；SMRTbell用接头将目标片段组织成可循环读取的模板。polymerase read涉及一次聚合酶观测过程，subread对应相应片段通过，共识CCS/HiFi则综合同一模板的重复观测。名词不能混作一个计数分母。[PacBio原理说明](https://www.pacb.com/technology/hifi-sequencing/how-it-works/)和[HiFi原始研究](https://doi.org/10.1038/s41587-019-0217-9)给出这一逻辑；产品化处理仍需按软件版本核对。

直觉是“同一页纸读五遍”，可利用多次观测识别偶发误读，但这不等于“抽到五页不同原稿”。覆盖度用于群体/变异推断时，应按最终独立分子与有效比对范围计算，不能把同一分子的所有pass当新模板叠加。

## 一个刻意简化的共识算术

设每次观测只有对/错两种状态，错误独立且概率$p=0.1$，五次多数决出错要求至少三次错。因此

$$
P_{\mathrm{maj}}=\sum_{k=3}^{5}{5\choose k}p^k(1-p)^{5-k}=0.00856.
$$

$k$为错误次数，组合数表示错误位置的选法。数值比0.1小，说明重复观测有用；但真实测序有四类替换、Indel、序列上下文、相关误差和质量模型，CCS也不是这个二元多数投票。尤其系统误差可在多次观测中重复，不能用这个0.856%当任何产品的准确率。

若总可用聚合酶读取长度粗略为$L_p$、插入长度为$L_i$，忽略接头及损耗时，pass数量大致随$L_p/L_i$变化。这个比例解释长度与重复观测之间的制约，不是为真实文库精确预测pass。过度追求长片段，可能改变可用重复次数及最终质量/产量；实际取舍应看运行和过滤结果。

## 输入输出和下游能做什么

保存模板来源、片段长度分布、试剂/软件、最终HiFi序列与质量及任务所需的BAM标签。动力学信息可支撑相应修饰模型，但是否保留原生DNA、是否经过扩增、模型与覆盖均影响可用性。Iso-Seq通常涉及cDNA来源，不等于直接RNA的原生电信号路线。

HiFi有助于跨越某些重复、支持单倍型和SV，但不能跨越所有重复，也不消除组装错误。历史CLR、subreads与最终HiFi的错误模型不同，工具预设应按输入类型选择，而非都套“PacBio”。

| 症状 | 可能原因 | 检查 |
|---|---|---|
| 共识数少于原始观测数 | 多pass合并或质量筛选 | 输入/输出计数单位和过滤日志 |
| 长片段质量分布改变 | 可用pass与运行制约 | 长度、pass、质量联合分布 |
| 覆盖异常高 | 把同一模板重复观测累计 | 最终read与分子定义 |
| 修饰信号缺失 | 模板处理或标签未保留 | 原生/扩增来源、数据和模型 |

## 自测与答案

1. 五次pass等于五个独立样本吗？答：此处只有一个模板分子。
2. toy算出的错误率能否直接写为HiFi实测？答：不能，假设和算法均不同。
3. 最长read是不是最佳输入？答：取决于任务，需同时看质量、覆盖、重复跨越和有效分子数。

## 在电脑上核对这一小步

先[下载项目ZIP](https://github.com/SpiralKernel/bioinfo-wiki/archive/refs/heads/main.zip)并解压，进入同时含README.md和examples的项目根目录。不熟悉终端和当前目录时读[R与Python路径文件与项目组织](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E4%B8%8EPython%E8%B7%AF%E5%BE%84%E6%96%87%E4%BB%B6%E4%B8%8E%E9%A1%B9%E7%9B%AE%E7%BB%84%E7%BB%87)。以下在**系统终端**执行，需Python 3.10+，只用标准库，无需pip安装；Windows可将python3替换为已安装的python或py。

```bash
python3 examples/methods_omics_demo.py --case hifi
```

输入是脚本内置人工数字，不需下载专业数据。成功会打印JSON，并生成`examples/output/methods-omics/hifi/summary.json`；`cases.hifi`下的字段对应本页演算。非空输出目录会拒绝覆盖，第二次练习可在命令后加`--output examples/output/methods-omics/hifi-practice2`。不要为了重跑删除不明来源的数据。

**检查点：**passes=5，majority_error=0.00856，original_molecules=1。

阅读源码时先找到同名函数，逐一对照输入、运算与输出；要修改数字，先复制整个解压项目并运行副本，再手算新的预期。脚本含若干固定结果断言，修改教学条件时也需重新核对相应断言；断言失败不自动说明专业方法失效。本轮已实际运行全部22个标准库模式，Python 3.14.8；本页范围为独立二元错误多数决；未生成真实CCS/HiFi或执行修饰分析。演示没有为真实样本产生专业P值或性能结论。


## 来源与下一步

本轮查阅日期为2026-10-26；正文原始来源支持对应技术事实，手算与代码是本项目人工教学设计。方法页解释模型，安装及大数据运行交给相应工具/工作流页。教学组织参考[GTN](https://training.galaxyproject.org/training-material/)的步骤与检查点、[Bioinformatics Workbook](https://bioinformaticsworkbook.org/)的问题驱动练习；未复制其教程正文。

保留的相关入口：[HiFi官方原理](https://www.pacb.com/technology/hifi-sequencing/how-it-works/)、[Wenger等原始研究](https://www.nature.com/articles/s41587-019-0217-9)、[长读长测序分析入门](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E9%95%BF%E8%AF%BB%E9%95%BF%E6%B5%8B%E5%BA%8F%E5%88%86%E6%9E%90%E5%85%A5%E9%97%A8)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
