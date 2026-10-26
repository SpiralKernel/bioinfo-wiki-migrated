# BLAST原理与种子延伸

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26核对原始与官方资料；本页Python微型计算已执行；人工数据不代表专业软件或真实研究性能。

BLAST用于在序列数据库中寻找局部相似片段。它先找可能相似的短词，再延伸和评分，减少对所有位置做完整动态规划的成本。本页解释速度来自哪里、何时会漏命中，以及为什么“最高分”与“最显著”必须结合检索背景理解。

前置为[Smith-Waterman局部比对](https://github.com/SpiralKernel/bioinfo-wiki/wiki/Smith-Waterman%E5%B1%80%E9%83%A8%E6%AF%94%E5%AF%B9)和[排序查找与哈希在生信中的应用](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%8E%92%E5%BA%8F%E6%9F%A5%E6%89%BE%E4%B8%8E%E5%93%88%E5%B8%8C%E5%9C%A8%E7%94%9F%E4%BF%A1%E4%B8%AD%E7%9A%84%E5%BA%94%E7%94%A8)。读完应能用人工反例说明种子长度的灵敏度取舍，并理解E值的基本含义。

## 从短词命中到候选比对

可将流程理解为：定义查询和数据库→生成短词及索引→寻找词命中→延伸高潜力区域→合并/细化候选→计算分数与统计量。现代BLAST有不同task、蛋白邻域词、双命中等策略；不能把“所有任务都找完全相同DNA词”当通用事实。[NCBI BLAST官方帮助](https://blast.ncbi.nlm.nih.gov/doc/blast-help/)说明不同程序和参数；原始方法见[Altschul等1990](https://doi.org/10.1016/S0022-2836(05)80360-2)。

下面的教学程序只建立完全相同k字符词的字典，然后作不含缺口的双向延伸；未实现完整BLAST、替换矩阵、统计显著性、反向互补或数据库格式。

## 为什么种子可能漏掉最优片段

query=`ABCD`，subject=`AXCY`，k=2时query的词AB、BC、CD，在subject的AX、XC、CY中一个也没有，故没有候选。但单个A或C仍能产生局部正分；精确局部DP能够找到，种子法可能完全跳过。减小k会增加命中与计算负担，同时更容易接近弱相似片段。较大k通常减少偶然词命中，但也可能漏更远缘或高错误率序列。

> 本页例子在项目根目录的Python中运行；标准库例子无需安装包，含NumPy/SciPy的例子先按[R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83)建立项目环境。Windows命令使用`python`，macOS/Linux可使用`python3`；代码块可复制进该环境的Python交互会话，或保存为项目根目录下的临时脚本再用该环境解释器运行；不要保存到examples目录或任意外部目录后假定项目导入路径仍相同。除特别说明外，只向屏幕输出，不写文件。

```python
from examples.alignment_demo import seed_extend, smith_waterman
hit = seed_extend('ABCD', 'XXABCDYY', k=2)[0]
print('simple_hit', hit)
print('seed_miss', seed_extend('ABCD', 'AXCY', k=2))
print('exact_local', smith_waterman('ABCD', 'AXCY')['score'])
assert hit['score'] == 8
assert seed_extend('ABCD', 'AXCY', k=2) == []
assert smith_waterman('ABCD', 'AXCY')['score'] > 0
```

检查点：第一例8分、query[0:4]对subject[2:6]；第二例种子列表为空，而精确局部分数大于0。两种算法没有计算同一个完整候选空间，输出不同不必是程序bug。

## 延伸为什么不能无限进行

设延伸的当前累计分数为S，见过的最佳分数为Sbest。当Sbest−S超过X-drop阈值，继续延伸的潜力被认为较小，可停止。阈值小则更快但更易在暂时低相似区断开；阈值大则允许跨过更长低分区。本项目无缺口小实现采用这一思想，并只保留最佳延伸端点；不要据此解释所有BLAST task的内部默认参数。

匹配+2、错配−1时，两个匹配先到4，再出现三个错配到1，分数下降3。X-drop=2会在下降超过2时停止；X-drop=3是否恰好继续取决于采用“大于”还是“大于等于”的边界，这也是实现需记录的细节。

## 读懂输出时必须区分的量

| 输出 | 回答什么 | 不回答什么 |
|---|---|---|
| identity | 已对齐列中相同字符的比例，具体分母按格式核对 | 是否覆盖完整查询 |
| query coverage | 查询有多大部分被比对覆盖 | 命中功能是否相同 |
| raw/bit score | 评分模型下的局部相似强度及标准化分值 | 后验功能概率 |
| E值 | 该检索背景下预期偶然同等或更高分命中的数量 | “这个命中为假的概率” |
| HSP与坐标 | 哪些局部片段相似、方向如何 | 必然一对一正交关系 |

常用直觉关系为E约随有效搜索空间增长、随分数指数下降；理想化写作$E\approx Kmn e^{-\lambda S}$。m、n是有效查询与库长度，S是原始分数，K、λ与评分/背景有关；实际工具还有校正。相同分数查更大的库，E值可能变大。不能自行给toy分数套未经标定的K、λ产生“BLAST E值”。

## 实操选项与排错

真实检索先选择核酸对核酸、蛋白对蛋白或翻译搜索；确认查询、库版本、task、低复杂度过滤和输出字段。命中太少时先查序列类型、方向、task是否过于严格；不要第一步就取消全部过滤。命中很多且短时查低复杂度、重复和覆盖度。输出前若被设置最大返回数截断，未出现的候选不能被解释成不存在。

**自测：100% identity但只覆盖20bp的命中，足以注释一条2kb基因吗？** 答：不足，至少还需覆盖、结构和功能证据。**自测：减小种子一定得到更好的生物学结论吗？** 答：只改变搜索灵敏度和成本；仍要评价特异性、统计背景与研究问题。生产操作入口保留在相关工具与工作流页。


## 接着读什么

[Smith-Waterman局部比对](https://github.com/SpiralKernel/bioinfo-wiki/wiki/Smith-Waterman%E5%B1%80%E9%83%A8%E6%AF%94%E5%AF%B9) → [排序查找与哈希在生信中的应用](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%8E%92%E5%BA%8F%E6%9F%A5%E6%89%BE%E4%B8%8E%E5%93%88%E5%B8%8C%E5%9C%A8%E7%94%9F%E4%BF%A1%E4%B8%AD%E7%9A%84%E5%BA%94%E7%94%A8) → [R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83) → [比对相似性同源性与覆盖率](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%AF%94%E5%AF%B9%E7%9B%B8%E4%BC%BC%E6%80%A7%E5%90%8C%E6%BA%90%E6%80%A7%E4%B8%8E%E8%A6%86%E7%9B%96%E7%8E%87) → [BLAST本地检索实践](https://github.com/SpiralKernel/bioinfo-wiki/wiki/BLAST%E6%9C%AC%E5%9C%B0%E6%A3%80%E7%B4%A2%E5%AE%9E%E8%B7%B5) → [方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。学习顺序与练习入口见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。


## 运行本页的完整小例

先下载项目并进入含README与examples的根目录，环境准备见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。用项目环境的Python运行以下命令（macOS/Linux路径示例；Windows用`.venv\Scripts\python.exe`）：

```bash
.venv/bin/python examples/methods_foundations_demo.py --case blast --output examples/output/methods-blast
```

输出summary.json保存本页计算与检查结果。目录非空时拒绝覆盖，再练习请更换`--output`目录；页面内短代码只打印结果。标准库case也可用系统Python直接运行。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
