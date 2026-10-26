# Sanger测序原理与峰图判读

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26官方/原始来源核对；标准库教学演示实际执行；专业验证范围见正文。

## 从测量对象开始

Sanger适合有限目标片段的序列确认，但峰图同时混合来自许多相关模板的信号。学完应能区分单点混峰、Indel后的错位混峰与低质量，不把输出字母当成全部证据。先读[正负链坐标与转录方向](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%AD%A3%E8%B4%9F%E9%93%BE%E5%9D%90%E6%A0%87%E4%B8%8E%E8%BD%AC%E5%BD%95%E6%96%B9%E5%90%91)、[SNVIndelMNV与复杂变异](https://github.com/SpiralKernel/bioinfo-wiki/wiki/SNVIndelMNV%E4%B8%8E%E5%A4%8D%E6%9D%82%E5%8F%98%E5%BC%82)。

## 链终止怎样变成有顺序的峰

链终止核苷酸ddNTP缺少继续正常延伸所需的3′羟基。产生在不同位置终止的片段后，按长度分离并识别末端碱基，便可恢复延伸链序列。[Sanger等1977原始论文](https://pubmed.ncbi.nlm.nih.gov/271968/)提出链终止方法；现代常见毛细管荧光读出把不同长度产物转为色谱轨迹，历史方法与当前仪器不能完全等同。

峰图横轴是迁移/检测位置，经过模型得到碱基索引，不能直接视为参考基因组坐标。四通道的峰形、峰间距、背景及质量共同支持判读。一个主峰较高不表示模板一定纯合，另一等位可能扩增不足；小次峰也不必然是低频真实变异。

## 手算一个错位混合

参考延伸序列为ACGTTACG，另一等位少一个T成为ACGTACG。按两条read的合成序号相叠，前4位一样，第5位开始分别出现T/A、A/C、C/G。持续混峰可由Indel造成的错位解释，不必假设连续三个独立杂合SNV。真实重复序列中可能暂时重合，所以不是“Indel后每位必定双峰”。

用信号示意 $I_b(t)=\sum_h a_h F_b(t;h)+\epsilon_b(t)$。$b$是碱基通道，$t$是检测位置，$h$为模板类型，$a_h$为相应信号权重，$F$为片段形成与分离后的响应，$\epsilon$为背景。权重受扩增、染料和序列影响，并非直接等于细胞或等位比例，因此不能仅凭峰高1:1保证基因型，也不能用toy信号建立统一低频检出限。

## 实际峰图的阅读顺序

先看整个AB1轨迹及可判读范围，避开引物附近或末端的低质区域；再定位候选附近的单峰/混峰与质量；将正反向结果统一到参考正链，检查反向互补；最后核对模板来源、PCR产物特异性和独立证据。[官方排查指南](https://www.thermofisher.com/TFS-Assets/LSG/manuals/MAN0014435_Trbleshoot_Sanger_seq_data_UB.pdf)可按症状继续定位问题，本页未读取真实AB1。

| 症状 | 解释候选 | 核对动作 |
|---|---|---|
| 某点两峰，前后清楚 | 杂合或混模板 | 双向结果、独立扩增与比对 |
| 某点后持续重叠 | Indel错位或混产物 | 两等位对齐、产物身份和方向 |
| 全程背景高 | 非特异产物、污染或样本问题 | 原始峰图、对照与模板记录 |
| 单向支持而反向不清楚 | 方向/局部质量差异 | 可读区间与反向互补，保留未知 |

read修剪不能依据“为了让候选看起来干净”任意选择；先定义质量与可判读规则。保存AB1、basecall、引物方向、参考版本、候选坐标、不确定区和验证是否独立。重复读取同一PCR产物与另一样本/技术验证的证据层级不同。

## 自测与答案

1. 输出ACG是否证明底层只有一个模板？答：否，basecaller可能用主信号选择字母。
2. 反向read的A能直接与参考正链A比较吗？答：需先统一方向，反向互补后再核对坐标。
3. 看到低峰即可确认低频突变吗？答：需排除背景和误差，并建立该实验的检出范围。

## 在电脑上核对这一小步

先[下载项目ZIP](https://github.com/SpiralKernel/bioinfo-wiki/archive/refs/heads/main.zip)并解压，进入同时含README.md和examples的项目根目录。不熟悉终端和当前目录时读[R与Python路径文件与项目组织](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E4%B8%8EPython%E8%B7%AF%E5%BE%84%E6%96%87%E4%BB%B6%E4%B8%8E%E9%A1%B9%E7%9B%AE%E7%BB%84%E7%BB%87)。以下在**系统终端**执行，需Python 3.10+，只用标准库，无需pip安装；Windows可将python3替换为已安装的python或py。

```bash
python3 examples/methods_omics_demo.py --case sanger
```

输入是脚本内置人工数字，不需下载专业数据。成功会打印JSON，并生成`examples/output/methods-omics/sanger/summary.json`；`cases.sanger`下的字段对应本页演算。非空输出目录会拒绝覆盖，第二次练习可在命令后加`--output examples/output/methods-omics/sanger-practice2`。不要为了重跑删除不明来源的数据。

**检查点：**aligned_by_read_index中第1–4项mixed为false，第5–7项为true；反向互补为CGTAACGT。

阅读源码时先找到同名函数，逐一对照输入、运算与输出；要修改数字，先复制整个解压项目并运行副本，再手算新的预期。脚本含若干固定结果断言，修改教学条件时也需重新核对相应断言；断言失败不自动说明专业方法失效。本轮已实际运行全部22个标准库模式，Python 3.14.8；本页范围为字符串错位与方向核对；没有AB1峰图、真实荧光或湿实验验证。演示没有为真实样本产生专业P值或性能结论。


## 来源与下一步

本轮查阅日期为2026-10-26；正文原始来源支持对应技术事实，手算与代码是本项目人工教学设计。方法页解释模型，安装及大数据运行交给相应工具/工作流页。教学组织参考[GTN](https://training.galaxyproject.org/training-material/)的步骤与检查点、[Bioinformatics Workbook](https://bioinformaticsworkbook.org/)的问题驱动练习；未复制其教程正文。

保留的相关入口：[Sanger等1977论文](https://pubmed.ncbi.nlm.nih.gov/271968/)、[官方峰图排查指南](https://www.thermofisher.com/TFS-Assets/LSG/manuals/MAN0014435_Trbleshoot_Sanger_seq_data_UB.pdf)、[变异与群体分析的验证设计](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%8F%98%E5%BC%82%E4%B8%8E%E7%BE%A4%E4%BD%93%E5%88%86%E6%9E%90%E7%9A%84%E9%AA%8C%E8%AF%81%E8%AE%BE%E8%AE%A1)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
