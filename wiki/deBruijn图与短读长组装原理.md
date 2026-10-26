# deBruijn图与短读长组装原理

> 更新：2026-10-26；状态：可用；验证范围：官方/原始资料核对；本页标准库Python人工例已执行；未运行生产软件或真实生物数据。

de Bruijn图将序列片段的重叠关系压缩为图结构，是理解短读长组装的重要入口。它不依赖先有完整参考，但重复、错误和覆盖空缺会让图产生歧义。本页用6bp人工序列构图和拼接，再说明为什么实际组装远比“沿箭头走”复杂。

前置为[图与最短路径的生信入门](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%9B%BE%E4%B8%8E%E6%9C%80%E7%9F%AD%E8%B7%AF%E5%BE%84%E7%9A%84%E7%94%9F%E4%BF%A1%E5%85%A5%E9%97%A8)、[排序查找与哈希在生信中的应用](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%8E%92%E5%BA%8F%E6%9F%A5%E6%89%BE%E4%B8%8E%E5%93%88%E5%B8%8C%E5%9C%A8%E7%94%9F%E4%BF%A1%E4%B8%AD%E7%9A%84%E5%BA%94%E7%94%A8)。读完应能明确节点/边约定，判断一个分叉可能来自哪些来源，而不把所有分叉都称为真实变异。

## 先固定一种图定义

本页以每个(k−1)-mer为节点，每个观测k-mer为从其前k−1字符到后k−1字符的有向边。k=3时，ATG给AT→TG。另一常见表示以k-mer为节点、重叠为边；两者相关但下标不同，读文献时需先确认约定。重复出现的相同k-mer可保存为一条带计数边，或多重边，不可静默丢掉覆盖信息。

理想情况下若所有必要边与次数准确、序列无歧义，可以寻找使用每条边规定次数的欧拉路径恢复序列；真实read覆盖会重复观察同一真实k-mer，因此观察次数不等于应在组装中走的次数。专业组装要估计覆盖、清理错误与解析重复，不能直接把所有read边重复走一遍。原始应用见[Velvet研究](https://doi.org/10.1101/gr.074492.107)和[作者仓库与手册](https://github.com/dzerbino/velvet)。

## ATGACG怎样变成路径

k=3，窗口依次ATG、TGA、GAC、ACG，对应AT→TG→GA→AC→CG。从起始节点AT开始，每条边追加目标节点最后字符，得到ATGACG。节点本来已有k−1字符，所以长度=边数+k−1=4+2=6。

```python
from collections import Counter
sequence, k = 'ATGACG', 3
kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
edges = [(word[:-1],word[1:]) for word in kmers]
counts = Counter(edges)
assembled = edges[0][0]+''.join(v[-1] for u,v in edges)
print('kmers',kmers)
print('edges',edges)
print('assembled',assembled,'length',len(assembled))
assert assembled==sequence
branch_words = ['ATG','ATC']
branches = {w[1:] for w in branch_words if w[:-1]=='AT'}
print('AT_outneighbors',sorted(branches))
assert len(branches)==2
```

这段代码的边顺序来自已知原序列，所以只验证“由已知路径拼接字符串”，没有从打乱reads真正求解组装。后一例ATG和ATC给同一前缀AT的两个出边，真实read图若有该分叉，需要额外证据决定路径。

## k为什么同时影响连通与重复

较小k使不同重复更容易共用节点，图更连通但歧义可能多；较大k有助区分重复，却更容易受覆盖不足和错误影响而断开。单个碱基错误可影响最多k个跨过该位置的k-mer，边界处更少。低频k-mer可能是错误，也可能是真实低覆盖或低丰度序列，不能一律删除。

长于read或缺乏跨越证据的重复，单靠短片段可能无法唯一解析。双端距离、长读长、连锁或其他信息可帮助连接，但不能凭图路径“看起来顺”就认定是真实染色体。异质二倍体会有单倍型分支，宏基因组还有物种覆盖差，不能按单一均匀覆盖假设清理。

## 常见图形意味着什么

| 图现象 | 可能来源 | 要补的证据 |
|---|---|---|
| tip短死端 | 错误、低覆盖、真实末端 | 碱基质量、覆盖、端点支持 |
| bubble分叉后汇合 | SNP/Indel、错误、重复 | 支持量、相位与长度 |
| 长重复节点多进多出 | 多拷贝重复 | 跨越reads、配对与覆盖 |
| 图断开 | 覆盖缺口、k过大、真实不同片段 | read分布、k敏感性与污染 |

“bubble一定是杂合”“tip一定是错”都过度简化。清图阈值可能去掉真实低频事件；保留所有边又可能导致错误拼接。参数与输入质量需要一起诊断。

## 从图到contig再到质量判断

无分叉路径可压缩为unitig，图结构简化后更容易处理。contig较长或N50更大不保证更准确，错误拼接也可能抬高长度。应结合read回贴、覆盖、保守基因完整性、参考比较及相位/污染检查。参考比较也受参考本身结构差异影响，见[组装比较共线性与变异方向](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%BB%84%E8%A3%85%E6%AF%94%E8%BE%83%E5%85%B1%E7%BA%BF%E6%80%A7%E4%B8%8E%E5%8F%98%E5%BC%82%E6%96%B9%E5%90%91)。

本页只执行人工构图算术，没有安装Velvet或运行真实组装。专业工作流需要单独记录read类型、错误校正、k策略、资源与结果QC；图的教学实现不能替代这些验证。

## 自测与排错

**问：把相同k-mer观测十次就应在最终序列中出现十次吗？** 答：不一定，覆盖重复不是基因组拷贝数真值。**问：ATG和ATC的分叉足以确认一个SNP吗？** 答：不足，要排除错误、重复和多来源。**问：增大k一定提高组装质量吗？** 答：不一定，可能打断低覆盖区域。若拼接长度不等于边数+k−1，先查起始节点长度、每次追加字符数与是否误重复了整个k-mer。

本页Python代码在项目根目录执行，仅用标准库且只打印结果，不生成文件；Windows用`python`，macOS/Linux可用`python3`。学习顺序见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。


## 运行本页的完整小例

先下载项目并进入含README与examples的根目录，环境准备见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。用项目环境的Python运行以下命令（macOS/Linux路径示例；Windows用`.venv\Scripts\python.exe`）：

```bash
.venv/bin/python examples/methods_foundations_demo.py --case debruijn --output examples/output/methods-debruijn
```

输出summary.json保存本页计算与检查结果。目录非空时拒绝覆盖，再练习请更换`--output`目录；页面内短代码只打印结果。标准库case也可用系统Python直接运行。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
