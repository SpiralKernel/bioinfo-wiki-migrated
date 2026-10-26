# DNBSEQ测序与滚环扩增

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26官方/原始来源核对；标准库教学演示实际执行；专业验证范围见正文。

## 从问题进入

一条DNA太少而不易被光学读出时，可以把同一模板的信息重复呈现。DNBSEQ使用DNA纳米球DNB作为可读单元。本页要分清原始分子、圈内重复与reads，解释为何DNB形成阶段不等于整条实验链路没有PCR。前置读[文库接头索引与条形码](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%87%E5%BA%93%E6%8E%A5%E5%A4%B4%E7%B4%A2%E5%BC%95%E4%B8%8E%E6%9D%A1%E5%BD%A2%E7%A0%81)、[PCR与扩增基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/PCR%E4%B8%8E%E6%89%A9%E5%A2%9E%E5%9F%BA%E7%A1%80)。

## 圆环、串联复制和阵列分别做什么

带适配结构的DNA经相应处理成为环状模板；滚环扩增RCA沿同一个环重复复制，形成包含多份模板信息的产物，再形成DNB并在阵列上读出。环化、DNB生成、装载、循环信号采集与basecalling是不同环节。[MGI技术说明](https://us.global-mgitech.com/technologies/dnbseq-technology/)给出了这一结构及cPAS读出入口。cPAS是相关化学体系名称，实际试剂版本的信号编码与处理应查该产品说明，不能把所有历史方案合并为一种固定流程。

把环想成一张绕回起点的纸带。复印很多遍可以增强同一序列的信号，却没有抽到新的原始分子。常规PCR以新产物继续参与扩增的指数谱系，与沿环状模板反复复制的方式不同；但是建库在到达环化以前完全可能经历PCR或靶向扩增。

## 一个能手算的分子账本

原样本只有4个独立圆环，3个带A、1个带G。假设每个DNB都含100份模板重复，则有400份模板信息，但独立原始分子仍是4，G比例仍为1/4。若建库选择使A来源变成6份、G只有1份，之后再理想地放大每个环，G比例也只是1/7≈0.143；DNB本身不会知道原样本应为1/4并把比例纠正回来。

可把这一抽象写成 $S_{ib}=a_i n_i u_{ib}+\epsilon_{ib}$：$S_{ib}$是第$i$个位置在某轮碱基通道$b$的观测信号，$n_i$是有效模板重复量，$u_{ib}$是与当前碱基及化学有关的理想响应，$a_i$综合装载/成像响应，$\epsilon$是噪声。它只是解释信号增强的线性示意，不是MGI专有basecaller。增加$n_i$不保证所有误差都下降，背景、混合模板和上游选择可能保持。

## 输入输出及实际检查顺序

输入需要平台兼容的文库结构、索引方案和样本清单；输出常见为拆分后的R1/R2 FASTQ及运行质量报告。FASTQ统一了序列与质量格式，并未抹去接头、索引方向和误差上下文差异。第三方文库转换应逐项确认结构与拆分规则，不凭“都叫双端150”判断兼容。

实际分析先核对样本/索引分配，再看每循环质量、接头及有效长度，最后检查比对、覆盖、重复和下游可调用区域。重复率还受文库复杂度、测序深度与重复定义影响，不能直接把一个数字归因于平台原理。

| 现象 | 可能解释 | 检查动作 |
|---|---|---|
| 大量reads未分配样本 | 索引表或方向不匹配 | 核对试剂布局和拆分软件日志 |
| Q值不错但目标覆盖差 | 文库选择或比对偏差 | 看目标范围、GC及可比对区域 |
| G等位少于预期 | 上游扩增/捕获或抽样 | 对照原始分子量和等位支持 |
| 重复率高 | 复杂度低、深测或定义差异 | 分开报告原始reads、去重后及分子计数 |

## 自测与答案

1. 400份圈内重复是否等于400条独立生物学证据？答：此例只有4个原始圆环，谱系必须保留。
2. DNB阶段不用指数PCR，能否写整实验PCR-free？答：只有核对前序建库后才能说明具体阶段。
3. 比较DNBSEQ和其他平台应先匹配什么？答：输入样本、文库、有效数据量、目标区域和评估规则；不是先排列宣传准确率。

## 在电脑上核对这一小步

先[下载项目ZIP](https://github.com/SpiralKernel/bioinfo-wiki/archive/refs/heads/main.zip)并解压，进入同时含README.md和examples的项目根目录。不熟悉终端和当前目录时读[R与Python路径文件与项目组织](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E4%B8%8EPython%E8%B7%AF%E5%BE%84%E6%96%87%E4%BB%B6%E4%B8%8E%E9%A1%B9%E7%9B%AE%E7%BB%84%E7%BB%87)。以下在**系统终端**执行，需Python 3.10+，只用标准库，无需pip安装；Windows可将python3替换为已安装的python或py。

```bash
python3 examples/methods_omics_demo.py --case rca
```

输入是脚本内置人工数字，不需下载专业数据。成功会打印JSON，并生成`examples/output/methods-omics/rca/summary.json`；`cases.rca`下的字段对应本页演算。非空输出目录会拒绝覆盖，第二次练习可在命令后加`--output examples/output/methods-omics/rca-practice2`。不要为了重跑删除不明来源的数据。

**检查点：**original_circles=4、total_template_repeats=400，allele_fraction_G=0.25；上游偏差后biased_fraction_G=1/7。

阅读源码时先找到同名函数，逐一对照输入、运算与输出；要修改数字，先复制整个解压项目并运行副本，再手算新的预期。脚本含若干固定结果断言，修改教学条件时也需重新核对相应断言；断言失败不自动说明专业方法失效。本轮已实际运行全部22个标准库模式，Python 3.14.8；本页范围为分子账本模型；未运行DNB制备、cPAS识别或平台实测。演示没有为真实样本产生专业P值或性能结论。


## 来源与下一步

本轮查阅日期为2026-10-26；正文原始来源支持对应技术事实，手算与代码是本项目人工教学设计。方法页解释模型，安装及大数据运行交给相应工具/工作流页。教学组织参考[GTN](https://training.galaxyproject.org/training-material/)的步骤与检查点、[Bioinformatics Workbook](https://bioinformaticsworkbook.org/)的问题驱动练习；未复制其教程正文。

保留的相关入口：[MGI DNBSEQ](https://us.global-mgitech.com/technologies/dnbseq-technology/)、[第二代测序原理与平台比较](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%AC%AC%E4%BA%8C%E4%BB%A3%E6%B5%8B%E5%BA%8F%E5%8E%9F%E7%90%86%E4%B8%8E%E5%B9%B3%E5%8F%B0%E6%AF%94%E8%BE%83)、[测序错误重复偏差与嵌合体](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%B5%8B%E5%BA%8F%E9%94%99%E8%AF%AF%E9%87%8D%E5%A4%8D%E5%81%8F%E5%B7%AE%E4%B8%8E%E5%B5%8C%E5%90%88%E4%BD%93)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
