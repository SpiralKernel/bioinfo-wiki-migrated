# Nanopore测序与直接RNA

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26官方/原始来源核对；标准库教学演示实际执行；专业验证范围见正文。

## 先问测到的是哪种分子

Nanopore将分子过孔时的电流变化转换为序列。原生DNA、逆转录cDNA与直接RNA是不同观测链路，不能仅凭FASTQ里出现T/U就判断实验。学完应能说明原始信号、basecall和修饰预测的区别。先读[从DNA到RNA和蛋白质](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E4%BB%8EDNA%E5%88%B0RNA%E5%92%8C%E8%9B%8B%E7%99%BD%E8%B4%A8)、[FASTA与FASTQ](https://github.com/SpiralKernel/bioinfo-wiki/wiki/FASTA%E4%B8%8EFASTQ)。

## 电流为何不等于一个碱基一档电压

孔中局部序列共同影响电信号，马达控制运动，停留时间和噪声也进入观测。真实basecaller学习连续信号到序列的映射，而不是逐点查四种固定电流。官方[测序机制说明](https://nanoporetech.com/platform/technology)和[Dorado模型文档](https://software-docs.nanoporetech.com/dorado/latest/)分别承接物理读出与软件处理；模型必须与孔、化学和数据类型配套。

用人工两碱基窗口解释“上下文”：设AA、AC、CG的理想信号分别为10、15、20，观测为10.2、14.8、20.1。候选AACG的连续窗口是AA/AC/CG，平方误差和为 $0.2^2+(-0.2)^2+0.1^2=0.09$。候选AAAC给出AA/AA/AC，误差和49.09，故前者更接近。窗口共享碱基，所以不是把每个电流单独翻译后随意拼接。

形式上可写 $\hat{s}=\arg\max_s P(s\mid I,M)$，$s$为候选序列，$I$为整段信号，$M$为选定模型。上面的平方误差只在特定高斯噪声假设下对应似然比较；没有真实采样率、停留时间、神经网络或校准，不能拿它当Dorado替代品。

## DNA、cDNA、直接RNA分别保留什么

原生DNA路线可保留某些修饰相关信号；PCR复制会改变可观测的原生修饰信息。cDNA路线先逆转录，读取DNA产物；PCR-cDNA另有扩增谱系。直接RNA读取RNA相关信号，流程可包含辅助cDNA合成，但被测分子的身份仍要按该试剂流程说明。不同路线在端点覆盖、截短、转录本及修饰解释上不能混用。

直接RNA不是自动的全转录组无偏普查，也不保证每条RNA读到两端。降解、捕获、接头连接、长度和运行过程都影响哪些分子被看见。比较组间转录本时还应核对输入RNA质量、建库路线与可观测范围，而不只比较read N50。

## 文件与参数的交接

POD5保存原始信号；FASTQ保存basecall序列与质量；BAM还可能携带模型输出的修饰/运动等标签。仅保存FASTQ通常不足以重做依赖信号的分析。[POD5官方格式](https://github.com/nanoporetech/pod5-file-format)说明其用途。simplex和duplex利用的证据不同，两条链的组合不能算成两份独立原始DNA覆盖。

| 症状 | 可能原因 | 优先检查 |
|---|---|---|
| 换模型后质量/变异明显改变 | 模型或化学不匹配 | 孔、试剂、模型名和版本 |
| 长read很多但全长RNA少 | 截短、RNA质量或端点限制 | 转录本覆盖和两端支持 |
| FASTQ无法重做修饰分析 | 缺少原始信号/标签 | POD5、模型与BAM保留策略 |
| duplex数量变化被当产量翻倍 | 分子与read计数混用 | 配对来源及独立分子分母 |

## 自测与答案

1. 信号分类只需四个均值吗？答：局部上下文、时间与噪声使模型远复杂于四类查表。
2. 辅助生成cDNA就说明直接RNA读的是cDNA吗？答：需看实际传感分子，不能仅由辅助步骤推断。
3. Q值提高能证明修饰预测正确吗？答：碱基识别与修饰模型是不同目标，均需匹配验证。

## 在电脑上核对这一小步

先[下载项目ZIP](https://github.com/SpiralKernel/bioinfo-wiki/archive/refs/heads/main.zip)并解压，进入同时含README.md和examples的项目根目录。不熟悉终端和当前目录时读[R与Python路径文件与项目组织](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E4%B8%8EPython%E8%B7%AF%E5%BE%84%E6%96%87%E4%BB%B6%E4%B8%8E%E9%A1%B9%E7%9B%AE%E7%BB%84%E7%BB%87)。以下在**系统终端**执行，需Python 3.10+，只用标准库，无需pip安装；Windows可将python3替换为已安装的python或py。

```bash
python3 examples/methods_omics_demo.py --case nanopore
```

输入是脚本内置人工数字，不需下载专业数据。成功会打印JSON，并生成`examples/output/methods-omics/nanopore/summary.json`；`cases.nanopore`下的字段对应本页演算。非空输出目录会拒绝覆盖，第二次练习可在命令后加`--output examples/output/methods-omics/nanopore-practice2`。不要为了重跑删除不明来源的数据。

**检查点：**chosen=AACG；candidate_squared_error约0.09对49.09。不要把人工信号单位解释成真实pA。

阅读源码时先找到同名函数，逐一对照输入、运算与输出；要修改数字，先复制整个解压项目并运行副本，再手算新的预期。脚本含若干固定结果断言，修改教学条件时也需重新核对相应断言；断言失败不自动说明专业方法失效。本轮已实际运行全部22个标准库模式，Python 3.14.8；本页范围为人工上下文信号匹配；未运行Dorado、POD5读取、真实RNA或修饰识别。演示没有为真实样本产生专业P值或性能结论。


## 来源与下一步

本轮查阅日期为2026-10-26；正文原始来源支持对应技术事实，手算与代码是本项目人工教学设计。方法页解释模型，安装及大数据运行交给相应工具/工作流页。教学组织参考[GTN](https://training.galaxyproject.org/training-material/)的步骤与检查点、[Bioinformatics Workbook](https://bioinformaticsworkbook.org/)的问题驱动练习；未复制其教程正文。

保留的相关入口：[ONT技术](https://nanoporetech.com/platform/technology/)、[Dorado](https://software-docs.nanoporetech.com/dorado/latest/)、[POD5](https://github.com/nanoporetech/pod5-file-format)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
