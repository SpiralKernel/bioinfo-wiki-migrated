# SV合并匹配与群体再分型

> 更新：2026-10-26；状态：可用；验证范围：官方/原始来源核对；标准库教学算例已实跑；专业流程范围见正文。

跨样本SV矩阵最危险的错误，是把“该样本没有报告记录”写成“已经测得纯合参考”。可靠流程要先定义候选事件，再按各样本证据推断状态，最后组织文件。本页目标是区分事件匹配、群体调用、强制分型和样本列合并，并保留多等位及不可判定状态。

## 三个步骤和三个输入对象

候选合并从各样本发现结果得到事件集合；再分型用调用器支持的read或中间证据重新评估这些候选；文件合并只把兼容结果组织到统一样本列。普通VCF、Sniffles SNF及GATK gVCF各含不同信息，不能改后缀互换。

输入应冻结参考、坐标和样本ID，并记录每个工具的类型、长度与等位序列。输出不只是矩阵，还需原始记录到统一事件的映射，以及每个样本的GT、质量、支持和缺失原因。缺失保留未知；候选发现本身有偏差，再分型并不自动消除未发现事件的遗漏。

## 相近不等于同一个等位

同一事件在重复区可有不同锚点，独立插入也可能相距很近。匹配要联合类型、两端/方向、长度与可用序列，BND不能仅看一端POS。位置阈值设大可降低重复记录，却也会丢失真实多等位；阈值太小则使同一等位在群体里被拆成多个低频候选。[Truvari事件合并入口](https://github.com/ACEnglish/truvari)

一个纯几何反例：三个插入候选位置1000、1060、1120，距离阈值70。A近B、B近C，单链接会把三个全连成一个簇，但A到C距离120，超过阈值。软件使用的聚类或代表点策略决定是否发生这种链式合并；不能仅从“阈值70”推断簇内任意两点都小于70。

## DELLY与SNF的两种组织方式

[DELLY结构变异与群体分型](https://github.com/SpiralKernel/bioinfo-wiki/wiki/DELLY%E7%BB%93%E6%9E%84%E5%8F%98%E5%BC%82%E4%B8%8E%E7%BE%A4%E4%BD%93%E5%88%86%E5%9E%8B)给出先合并候选、逐样本用相同集合分型、再用兼容输出合并样本列的模板。最后的bcftools merge没有读取BAM推断GT，所以不负责把无记录变成参考。[DELLY官方流程](https://github.com/dellytools/delly)

Sniffles在单样本时生成SNF，群体模式整合SNF候选及中间信息；force模式则从BAM和已知VCF推断给定事件。它们不是同一个步骤的两种文件后缀。[Sniffles官方群体与force模式](https://github.com/fritzsedlazeck/Sniffles)

## 实测反例：缺失与低质量参考不能混在一起

[人工长读长贯通SV调用与评估](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E4%BA%BA%E5%B7%A5%E9%95%BF%E8%AF%BB%E9%95%BF%E8%B4%AF%E9%80%9ASV%E8%B0%83%E7%94%A8%E4%B8%8E%E8%AF%84%E4%BC%B0)的V01发现DEL/INS，R01单样本为0候选。SNF群体结果保留R01的`./.`；另在R01 BAM对V01候选force得到`0/0`，但GQ=0。前者未被群体模式判为参考，后者也不能重新描述为高置信阴性；应分别保留其证据和模型输出。

本机Sniffles2.8.1 force还发现FORMAT键少于样本值、样本表头继承候选的问题。项目保存raw，仅在版本及已知4键/5值末尾点模式符合时，派生.compat补PS；经BAM RG核对后记录样本改名。GT/GQ/DR/DV与INFO不改，原生与派生SHA可追溯。此修正不是通用VCF修复器。

force的INFO/SUPPORT、VAF、COVERAGE可来自候选原记录；当前样本证据需读FORMAT/DR、DV与GQ。仅看到PASS或INFO支持高就把R01判成高支持变异，会颠倒数据含义。

## 矩阵构建与敏感性分析

在事件映射稳定后，检查样本缺失率是否集中于平台、批次或某类长度；核对多等位拆分和剂量，不能把所有INS统一编码成一个0/1。对难以分型的位点可保留支持状态或unknown，但不要为了矩阵完整而自动填参考。

改变匹配距离、序列相似性和支持门槛时，比较事件数、多等位保留、缺失模式及个别read证据。重新计算频率要基于实际可用基因型/倍性，不能沿用合并前不同样本集的INFO/AF。用于PAV时，还需定义特征与剂量如何对应存在/缺失。

## 动手：运行后先对照什么

在下载项目的根目录执行；只需Python 3标准库，不安装专业统计包：

```bash
python3 examples/methods_genetics_demo.py --case merge
```

输出位置1000/1060/1120、相邻距离60/60、两端距离120，展示单链接可以形成一个簇。另保留`absent_record_gt=null`与`observed_force_gt=0/0、force_GQ=0`的不同语义。

脚本只演示链式邻近与状态表示，不合并真实VCF或运行任何SV分型器。实测SNF/force结果以链接案例及原始文件为准。

命令只向终端打印JSON。若要保存所有22个小例，可用`python3 examples/methods_genetics_demo.py --output examples/output/methods-genetics-practice1`；目录必须新建或为空，已有结果不会被覆盖。源代码见[标准库方法算例](https://github.com/SpiralKernel/bioinfo-wiki/blob/main/examples/methods_genetics_demo.py)。

## 排错：先定位原因再改参数

| 症状 | 可能原因 | 下一步检查 |
|---|---|---|
| 事件数减少但多等位消失 | 容差过大或链式合并 | 查看原ID到代表事件映射与序列差异 |
| 某平台几乎全是0/0 | 无记录被误填参考 | 回查分型证据与原生缺失 |
| force INFO与GT矛盾 | INFO可能继承候选 | 查当前样本FORMAT及genotypeCommand来源 |

## 自测与答案

1. **bcftools merge能从没有记录的样本推断参考GT吗？** 不能，它组织文件，不替代read层再分型。

2. **低GQ的0/0与./.应该统一填0吗？** 不应。二者是不同输出状态，仍须结合质量规则决定可用性。

## 原有实践与延伸入口

[Truvari工具与论文入口](https://github.com/ACEnglish/truvari)；[Sniffles长读长结构变异实践](https://github.com/SpiralKernel/bioinfo-wiki/wiki/Sniffles%E9%95%BF%E8%AF%BB%E9%95%BF%E7%BB%93%E6%9E%84%E5%8F%98%E5%BC%82%E5%AE%9E%E8%B7%B5)；[SV事件匹配与一对一评估](https://github.com/SpiralKernel/bioinfo-wiki/wiki/SV%E4%BA%8B%E4%BB%B6%E5%8C%B9%E9%85%8D%E4%B8%8E%E4%B8%80%E5%AF%B9%E4%B8%80%E8%AF%84%E4%BC%B0)；[Truvari事件比较与评估参数](https://github.com/SpiralKernel/bioinfo-wiki/wiki/Truvari%E4%BA%8B%E4%BB%B6%E6%AF%94%E8%BE%83%E4%B8%8E%E8%AF%84%E4%BC%B0%E5%8F%82%E6%95%B0)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
