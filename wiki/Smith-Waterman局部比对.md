# Smith-Waterman局部比对

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26核对原始与官方资料；本页Python微型计算已执行；人工数据不代表专业软件或真实研究性能。

局部比对寻找两条序列中得分最高的一对相似片段，允许忽略两端不相似部分。Smith–Waterman在给定评分模型下求精确局部最优，适合讲清“最佳片段”的定义；它本身不是快速全数据库检索，也不自动给出统计显著性。

前置为[动态规划与序列比对基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%8A%A8%E6%80%81%E8%A7%84%E5%88%92%E4%B8%8E%E5%BA%8F%E5%88%97%E6%AF%94%E5%AF%B9%E5%9F%BA%E7%A1%80)。读完应能解释公式里的0、知道从哪里回溯，并能检查局部区间是否对应实际对齐字符串。

## 相比全局比对，改了哪三个地方

定义H(i,j)为结束在两前缀处的最佳局部得分：

$$
H(i,j)=\max\{0,H(i-1,j-1)+s(a_i,b_j),H(i-1,j)+g,H(i,j-1)+g\}.
$$

其中a_i、b_j是两序列第i/j个字符，s是字符对的评分函数，g<0为线性缺口罚分；i/j从前缀长度计，0表示空前缀。边界行列都设0；允许在任一位置重新起步；回溯从整个矩阵的最大值开始，遇0停止。0的作用是丢弃会拖低得分的负前缀，而不是说所有错配都免费。已选局部片段内部的错配和缺口仍要付代价。[Smith–Waterman原始论文](https://doi.org/10.1016/0022-2836(81)90087-5)讨论了局部同源序列识别。

全局比对固定从左上走到右下，局部比对可在矩阵内部起止；两者目标不同，所以不能只看谁分数高。若匹配给正分、错配/缺口给负分，完全没有有利片段时最佳值是0，对齐可以为空。

## 先手算，再看程序

`ABCD`对`XXABCDYY`，匹配+2、错配−1、缺口−2。A与第3个字符A配对使分数从0起为2，随后B、C、D依次4、6、8。两侧XX和YY不进入最优局部片段；0基半开区间为query[0:4]、subject[2:6]。同分多个片段时，本教学实现按遍历顺序保留一个，不保证列出所有局部命中。

> 本页例子在项目根目录的Python中运行；标准库例子无需安装包，含NumPy/SciPy的例子先按[R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83)建立项目环境。Windows命令使用`python`，macOS/Linux可使用`python3`；代码块可复制进该环境的Python交互会话，或保存为项目根目录下的临时脚本再用该环境解释器运行；不要保存到examples目录或任意外部目录后假定项目导入路径仍相同。除特别说明外，只向屏幕输出，不写文件。

```python
from examples.alignment_demo import smith_waterman
r = smith_waterman('ABCD', 'XXABCDYY')
print(r)
assert r['score'] == 8
assert r['start'] == [0, 2] and r['end'] == [4, 6]
assert r['query'] == r['subject'] == 'ABCD'
gapped = smith_waterman('ABCDE', 'ABXCDE')
print('gap_example', gapped['score'], gapped['query'], gapped['subject'])
assert gapped['score'] == 8
assert smith_waterman('AAAA', 'BBBB')['score'] == 0
```

第一个输出为8分和两端区间；第二个为8分，最优含一个缺口，5次匹配得10再减2。去掉缺口应恢复各自选中的子串，而非整条原序列。字母是人工符号，不是已用蛋白替换矩阵分析的真实序列。

## 评分参数改变的实际含义

较重的错配罚分会让低相似区域更早断开；较轻的缺口罚分可能把两个匹配块接成更长片段。蛋白通常使用如BLOSUM的替换矩阵，不同氨基酸替换的分值不同；DNA也可能需要转移/颠换、碱基质量等模型。教学函数只有恒定匹配/错配和线性缺口，未实现仿射缺口。

输出至少保留原始score、起止区间、对齐字符串、方向和评分设定。端点在不同工具中可能用1基闭区间，复制本页0基端点到外部工具前必须换算。反向互补要明确比较并保存方向，算法不会从一条正链输入自动理解所有生物学关系。

## 算法精确，不等于命中有生物意义

同一评分下，DP对输入求得最优局部分数，但随机序列也可能产生正分；序列越长、检索库越大，偶然高分机会通常越多。显著性需要合适统计模型和背景。低复杂度片段可能带来高分而缺乏特异性；局部相似也不能单独证明正交、同功能或完整基因匹配。数据库检索的E值及启发式见[BLAST原理与种子延伸](https://github.com/SpiralKernel/bioinfo-wiki/wiki/BLAST%E5%8E%9F%E7%90%86%E4%B8%8E%E7%A7%8D%E5%AD%90%E5%BB%B6%E4%BC%B8)。

时间与矩阵存储都是O(mn)，所以它适合小片段精确核对，却不直接等价于生产软件的全基因组搜索效率。带状、SIMD、GPU等加速还涉及不同边界和实现，不能从本例耗时外推。

## 异常与自测

| 现象 | 可能原因 | 检查方法 |
|---|---|---|
| 输出为空或0分 | 没有正收益片段、参数太重 | 人工加入一段相同字符串 |
| 比对两端缺失很多 | 局部目标正常排除端点 | 核对是否其实需要全局目标 |
| 原始序列位置偏1 | 坐标基准不同 | 用已知ABCD区间核对 |
| 得分对但字符串错误 | 回溯/同分处理有问题 | 逐列重算、对照选中子串 |

**问：最大值一定在右下角吗？** 答：不一定，右下角表示固定端点，会漏掉真正局部最优。**问：8分相对于4分是否“显著两倍”？** 答：不是；原始分数需结合评分、长度和统计背景解释。


## 接着读什么

[动态规划与序列比对基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%8A%A8%E6%80%81%E8%A7%84%E5%88%92%E4%B8%8E%E5%BA%8F%E5%88%97%E6%AF%94%E5%AF%B9%E5%9F%BA%E7%A1%80) → [BLAST原理与种子延伸](https://github.com/SpiralKernel/bioinfo-wiki/wiki/BLAST%E5%8E%9F%E7%90%86%E4%B8%8E%E7%A7%8D%E5%AD%90%E5%BB%B6%E4%BC%B8) → [R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83) → [方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。学习顺序与练习入口见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。


## 运行本页的完整小例

先下载项目并进入含README与examples的根目录，环境准备见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。用项目环境的Python运行以下命令（macOS/Linux路径示例；Windows用`.venv\Scripts\python.exe`）：

```bash
.venv/bin/python examples/methods_foundations_demo.py --case smith-waterman --output examples/output/methods-smith-waterman
```

输出summary.json保存本页计算与检查结果。目录非空时拒绝覆盖，再练习请更换`--output`目录；页面内短代码只打印结果。标准库case也可用系统Python直接运行。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
