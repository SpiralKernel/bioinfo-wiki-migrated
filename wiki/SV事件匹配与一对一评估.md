# SV事件匹配与一对一评估

> 更新：2026-10-26；状态：可用；验证范围：官方/原始来源核对；标准库教学算例已实跑；专业流程范围见正文。

同一变异可以有不同锚点和序列表示，两个不同变异也可能非常靠近。因此“两个VCF共有多少行”通常不是SV准确度。评估先定义什么叫同一事件，再处理一个候选对应多个真值的歧义。本页要求能画候选匹配关系、手算TP/FP/FN，并解释匹配、基因型一致和序列准确各自的边界。

输入包括同样本/参考下的base真值、comp候选、范围BED及类型/长度和质量规则。先确认真值是否覆盖相应SV类型；范围外没有真值的记录不能自动判为假阳性。原始VCF、过滤后集合和被排除原因应分别保存。

## 从规则建立候选匹配边

对简单同类型事件，可检查两端距离、长度相似性与序列。示意长度相似度为 $s_L=\min(L_b,L_c)/\max(L_b,L_c)$；100bp与80bp为0.8。两个100bp插入即使位置完全相同，序列全不相同也不必视为同一等位。DEL区间重叠和INS/BND的连接关系不是同一规则。

Truvari使用明确的距离、大小、序列等参数和表示处理；序列比较并不是把原始ALT逐字符相等当全部定义。重复旋转或复合表示可影响匹配，相关功能需按版本核对。[Truvari官方bench说明](https://github.com/ACEnglish/truvari/wiki/bench)

## 一个候选不能被无限复用

设真值位置1000与1100，唯一候选1050，其他类型/长度条件均满足，距离容差60。候选与两条真值都够近。若对每条真值独立找“有没有近点”，会错误数出2个真阳性；在此一对一教学规则中最多匹配1对，所以TP=1、FP=0、FN=1，precision=1、recall=0.5。

可以把真值与候选画成二分图，边表示允许配对；选择不能重复使用节点的一组边。我们的toy用最大匹配演示计数约束，不实现Truvari的评分与挑选策略。真实复杂匹配可能有多种同样大小的配对，还需考察配对质量；不能只报告最高分而丢掉歧义记录。

## TP-base与TP-comp为什么分别报告

precision通常从候选侧衡量被确认的比例，recall从真值侧衡量被恢复的比例。一对一简单例两边TP相同，但分解、合并、多重匹配或refine之后的统计对象可能不同，不能把它们永远强行设为同一个整数。应读取工具的TP-base、TP-comp及对应VCF，而不是仅凭一个总分推导全部数量。

事件已匹配而GT从0/1变1/1，存在性与剂量结论不同；两个GT相同也不意味着边界或插入序列完全正确。相位、断点偏移、长度偏差、支持和过滤状态应按目标另行评价。未定义分母保留未定义，不写100%。

## 一次可复查的实际配置

[人工长读长贯通SV调用与评估](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E4%BA%BA%E5%B7%A5%E9%95%BF%E8%AF%BB%E9%95%BF%E8%B4%AF%E9%80%9ASV%E8%B0%83%E7%94%A8%E4%B8%8E%E8%AF%84%E4%BC%B0)使用Truvari5.5.0：`-b truth.vcf.gz -c V01.vcf.gz -f reference.fa --includebed callable.bed --passonly -s 50 -S 50 -r 100 -p 0.7 -P 0.7 --pick single`。其中距离100bp、序列与大小相似度0.7，两边相关长度门槛明确设为50；范围为人工chrToySV的BED[10000,90000)。

实际TP-base=TP-comp=2、FP/FN=0、匹配GT一致2。它验证两个已知易例的链路，不是全基因组准确率。改变阈值可能新增匹配，也可能错误合并相邻等位；应同时比较匹配对、被排除集合与错误类型，不挑得分最高的一次作为无条件性能。

## 参数敏感性与边界

模糊深度CNV和序列精确SV需要不同的容差解释；BND还需远端、朝向及mate关系。1对1记录比较难以涵盖所有复合事件，一条复杂真值与多条候选并不能靠增大距离无限解决。可按需要使用表示协调或更适合的事件级评估，并明确改变了什么对象。

结果表至少保存真值ID、候选ID、规则和匹配质量、GT比较、类型/长度/区域、FP/FN及排除原因。只有precision/recall两数字不足以发现表示或样本标签错误；更多检查见[变异评估真值集合与可检测区域](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%8F%98%E5%BC%82%E8%AF%84%E4%BC%B0%E7%9C%9F%E5%80%BC%E9%9B%86%E5%90%88%E4%B8%8E%E5%8F%AF%E6%A3%80%E6%B5%8B%E5%8C%BA%E5%9F%9F)。

## 动手：运行后先对照什么

在下载项目的根目录执行；只需Python 3标准库，不安装专业统计包：

```bash
python3 examples/methods_genetics_demo.py --case matching
```

输出一对一`TP=1、FP=0、FN=1、recall=0.5`；`naive_independent_base_hits=2`展示复用同一候选产生的错误计数。

这是仅按位置建立边的二分图算例，不解析VCF、不比较插入序列，也不是Truvari评分或--pick算法的替代实现。

命令只向终端打印JSON。若要保存所有22个小例，可用`python3 examples/methods_genetics_demo.py --output examples/output/methods-genetics-practice1`；目录必须新建或为空，已有结果不会被覆盖。源代码见[标准库方法算例](https://github.com/SpiralKernel/bioinfo-wiki/blob/main/examples/methods_genetics_demo.py)。

## 排错：先定位原因再改参数

| 症状 | 可能原因 | 下一步检查 |
|---|---|---|
| 候选比真值少却recall全满 | 同一候选被重复使用 | 检查逐事件匹配表和允许的一对多规则 |
| 相同位置被判不匹配 | 类型、序列或长度条件不同 | 逐项核对规则与原生表示 |
| 换阈值分数剧变 | 事件密集、模糊边界或筛选改变 | 比较新增/丢失匹配，不仅看汇总分数 |

## 自测与答案

1. **precision=1、recall=0.5矛盾吗？** 不矛盾；报出的候选都正确，但可能漏掉一半真值。

2. **GT一致是否等于插入序列完全正确？** 不是，基因型比较和等位序列准确性是不同层级。

## 原有实践与延伸入口

[Truvari事件比较与评估参数](https://github.com/SpiralKernel/bioinfo-wiki/wiki/Truvari%E4%BA%8B%E4%BB%B6%E6%AF%94%E8%BE%83%E4%B8%8E%E8%AF%84%E4%BC%B0%E5%8F%82%E6%95%B0)；[官方bench说明](https://github.com/ACEnglish/truvari/wiki/bench)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
