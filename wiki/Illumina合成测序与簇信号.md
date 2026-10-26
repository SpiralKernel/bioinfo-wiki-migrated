# Illumina合成测序与簇信号

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26官方/原始来源核对；标准库教学演示实际执行；专业验证范围见正文。

## 学习目标与前置

本页解释簇为何增强信号、每轮信号怎样变成碱基，以及为什么读长、质量和独立分子数不能互换。先读[文库接头索引与条形码](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%87%E5%BA%93%E6%8E%A5%E5%A4%B4%E7%B4%A2%E5%BC%95%E4%B8%8E%E6%9D%A1%E5%BD%A2%E7%A0%81)和[单端双端读长与插入片段](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%8D%95%E7%AB%AF%E5%8F%8C%E7%AB%AF%E8%AF%BB%E9%95%BF%E4%B8%8E%E6%8F%92%E5%85%A5%E7%89%87%E6%AE%B5)。这是一篇原理课，拆分样本和测序质控的具体运行放在工具与交接页。

## 从一个模板到一串碱基

适配文库在流动槽相应位置形成可读簇，同一簇中相关模板增强可检测信号。经典桥式扩增有助理解，但不同flow cell与产品的簇形成方式不完全相同，不能把所有型号写成同一桥式流程。文库PCR-free通常描述建库环节，并不表示仪器中所有模板呈现过程都没有放大。

SBS即边合成边识别：在循环中加入相应核苷酸，采集信号，处理后进入下一轮。碱基识别器根据多个观测与上下文估计序列及质量。通道编码随化学变化，四种碱基不要求始终各有一个独立颜色；具体实现以[Illumina SBS说明](https://www.illumina.com/science/technology/next-generation-sequencing/sequencing-technology.html)及对应运行文档为准。

## 同步为何影响长read

假设一个簇中模板每轮有比例$f$脱离理想同步，且只统计“从未脱离”的那部分，则$c$轮后比例为 $(1-f)^c$。$f=0.01,c=100$时约0.366。这里不是说真实仪器100轮后只有36.6%可读：实际phasing/prephasing、校正与信号模型更复杂。这个模型只说明小的逐轮偏差会累积，不能把早期和末端循环视为同样条件。

质量分数是另一层：$Q=-10\log_{10}p_e$，$p_e$为模型给出的碱基错误概率。Q20对应0.01，Q30对应0.001。150个均为Q30的碱基，期望错误数为150×0.001=0.15；这是期望值，不是每条read必有0.15个错误，也不证明错误独立或质量已经校准。Q30比例不涵盖污染、样本串扰、比对歧义和覆盖缺口。

## R1、索引、R2的账不能混

R1与R2观察同一文库片段的两端，不是两个独立取样的分子。若片段250bp、两端各150bp，中间会有50bp重叠；产生300个测序碱基不代表300bp独立基因组范围。index read主要用于样本拆分，某些文库还在read中携带UMI或细胞条码；应按结构解析，不能一律全当作目标序列比对。[官方索引入口](https://support-docs.illumina.com/SHARE/IndexedSeq/Content/SHARE/IndexedSequencing/IndexedSeqIntro.htm)承接不同布局。

## 参数、产物和排错

输入记录文库/试剂、单双端布局、索引和样本表；输出保留拆分FASTQ、每循环质量及有效read数。循环数增加可延伸read，也可能遇到片段末端接头和质量下降；装载与文库多样性问题不能仅靠延长read解决。

| 症状 | 候选原因 | 先查 |
|---|---|---|
| 末端质量下降 | 同步、化学或接头读穿 | 每循环质量和片段长度分布 |
| R1好而R2差 | 两端读取及文库差异 | 分端报告，不只合并Q值 |
| 低多样性文库表现异常 | 早期信号分布或样本结构 | 运行设置和实际文库设计 |
| 高比对率却覆盖不全 | 深测少数片段、重复 | 覆盖广度、重复定义及可调用范围 |

## 自测与答案

1. Q30是“某条read整体99.9%正确”吗？答：单碱基概率与整条read正确概率不同。
2. 一个簇含很多模板是否增加生物学重复？答：它们通常是相关模板，不能代替独立样本。
3. 将末端低质碱基剪掉必然改善一切吗？答：会丢长度与定位信息，应比较有效比对和任务结果，不能只追求质量数字。

## 在电脑上核对这一小步

先[下载项目ZIP](https://github.com/SpiralKernel/bioinfo-wiki/archive/refs/heads/main.zip)并解压，进入同时含README.md和examples的项目根目录。不熟悉终端和当前目录时读[R与Python路径文件与项目组织](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E4%B8%8EPython%E8%B7%AF%E5%BE%84%E6%96%87%E4%BB%B6%E4%B8%8E%E9%A1%B9%E7%9B%AE%E7%BB%84%E7%BB%87)。以下在**系统终端**执行，需Python 3.10+，只用标准库，无需pip安装；Windows可将python3替换为已安装的python或py。

```bash
python3 examples/methods_omics_demo.py --case illumina
```

输入是脚本内置人工数字，不需下载专业数据。成功会打印JSON，并生成`examples/output/methods-omics/illumina/summary.json`；`cases.illumina`下的字段对应本页演算。非空输出目录会拒绝覆盖，第二次练习可在命令后加`--output examples/output/methods-omics/illumina-practice2`。不要为了重跑删除不明来源的数据。

**检查点：**synchronous_fraction约0.366032；Q20/Q30错误概率为0.01/0.001；expected_errors_150bp_Q30=0.15。

阅读源码时先找到同名函数，逐一对照输入、运算与输出；要修改数字，先复制整个解压项目并运行副本，再手算新的预期。脚本含若干固定结果断言，修改教学条件时也需重新核对相应断言；断言失败不自动说明专业方法失效。本轮已实际运行全部22个标准库模式，Python 3.14.8；本页范围为同步衰减与Phred换算；未处理真实成像信号或运行拆分软件。演示没有为真实样本产生专业P值或性能结论。


## 来源与下一步

本轮查阅日期为2026-10-26；正文原始来源支持对应技术事实，手算与代码是本项目人工教学设计。方法页解释模型，安装及大数据运行交给相应工具/工作流页。教学组织参考[GTN](https://training.galaxyproject.org/training-material/)的步骤与检查点、[Bioinformatics Workbook](https://bioinformaticsworkbook.org/)的问题驱动练习；未复制其教程正文。

保留的相关入口：[文库接头索引与条形码](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%87%E5%BA%93%E6%8E%A5%E5%A4%B4%E7%B4%A2%E5%BC%95%E4%B8%8E%E6%9D%A1%E5%BD%A2%E7%A0%81)、[单端双端读长与插入片段](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%8D%95%E7%AB%AF%E5%8F%8C%E7%AB%AF%E8%AF%BB%E9%95%BF%E4%B8%8E%E6%8F%92%E5%85%A5%E7%89%87%E6%AE%B5)、[SBS官方说明](https://www.illumina.com/science/technology/next-generation-sequencing/sequencing-technology.html)、[索引说明](https://support-docs.illumina.com/SHARE/IndexedSeq/Content/SHARE/IndexedSequencing/IndexedSeqIntro.htm)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
