# BWT与FM索引的精确匹配

> 更新：2026-10-26；状态：可用；验证范围：官方/原始资料核对；本页标准库Python人工例已执行；未运行生产软件或真实生物数据。

BWT把字符串字符可逆地重新排列，FM索引利用其顺序与计数支持检索。它解释了许多短读长比对器为何可以重用参考索引，避免对每条read扫描整条基因组。本页只做精确子串匹配，不把它当成完整的带错配/缺口比对器。

前置为[排序查找与哈希在生信中的应用](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%8E%92%E5%BA%8F%E6%9F%A5%E6%89%BE%E4%B8%8E%E5%93%88%E5%B8%8C%E5%9C%A8%E7%94%9F%E4%BF%A1%E4%B8%AD%E7%9A%84%E5%BA%94%E7%94%A8)与[算法输入输出复杂度与数据结构](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%AE%97%E6%B3%95%E8%BE%93%E5%85%A5%E8%BE%93%E5%87%BA%E5%A4%8D%E6%9D%82%E5%BA%A6%E4%B8%8E%E6%95%B0%E6%8D%AE%E7%BB%93%E6%9E%84)。读完应能手算一个BWT、跟踪反向搜索区间，并区分“计数命中”和“定位命中”。

## 从后缀顺序得到BWT

给文本`banana`加唯一终止符`$`，约定`$`小于其他字符。后缀按字典序排列，其起点为SA=[6,5,3,1,0,4,2]。BWT的每项取该后缀起点前一个字符，起点0时取末尾`$`，得到`annb$aa`。SA按后缀而非原始坐标排序，所以命中在SA中连续，并不意味着原基因组位置连续。

BWT长度与原加终止符文本相同，本身只是排列，不自动压缩；相似上下文可使字符聚集，配合编码才可能压缩。作者讲义[Langmead：BWT与FM index](https://www.cs.jhu.edu/~langmea/resources/bwt_fm.pdf)详细定义SA与反向检索；[Bowtie原始研究](https://doi.org/10.1186/gb-2009-10-3-r25)说明它如何扩展到实际短读长匹配。

## 反向搜索的两个计数

令C(c)为全文中严格小于字符c的字符数；Occ(c,k)为BWT前k个字符中c的个数，采用半开前缀[0,k)。当前SA区间[l,r)表示匹配已处理查询后缀的所有行。向左再加字符c时：

$$
l'=C(c)+Occ(c,l),\qquad r'=C(c)+Occ(c,r).
$$

初始[0,N)，从查询最右字符向左处理。区间空则无精确命中；最终r−l为命中数。对ana：先a得[1,4)，再n得[5,7)，再a得[2,4)，对应SA位置3和1，即原文0基起点1与3。

```python
T = 'banana$'
SA = sorted(range(len(T)),key=lambda i:T[i:])
BWT = ''.join(T[(i-1)%len(T)] for i in SA)
C = {c:sum(x<c for x in T) for c in set(T)}
def occ(c,k):
    return BWT[:k].count(c)
l, r = 0, len(T)
for c in reversed('ana'):
    l, r = C[c]+occ(c,l), C[c]+occ(c,r)
    print('prepend',c,'interval',l,r)
hits = sorted(SA[l:r])
print('SA',SA,'BWT',BWT,'positions',hits)
assert BWT=='annb$aa' and hits==[1,3]
assert hits==[i for i in range(len(T)) if T.startswith('ana',i)]
```

检查点为上述三个区间、`annb$aa`和[1,3]。代码用切片计数与后缀字符串排序，易读但效率差；只演示正确性，没有实现压缩rank查询或高效建索引。

## 真正节省资源还依靠哪些结构

如果Occ每次都扫描前缀，查询仍很慢。实际索引利用检查点、rank结构或压缩表示快速计数；定位坐标可借助采样SA和LF映射，而不是保存所有后缀字符串。每步成本取决于rank实现、采样和字母表，不能只看公式就保证任意Python实现O(查询长度)。

索引通常针对固定参考构建一次，许多reads重复使用。索引参考版本、contig命名和序列必须与后续分析一致；相同文件名不意味着相同参考。磁盘索引也不等于运行峰值内存，两阶段应分别估计。

## 精确匹配不是完整read mapping

测序错误与真实变异会导致查询没有完整精确命中。实际比对器还需要允许错配/缺口、种子、回溯或局部延伸、双端约束、方向与MAPQ等规则。FM索引检索出多个候选，也不自动确定正确比对位置；重复区域尤需解释。长读长软件又可能采用minimizer等不同索引/种子策略。

因此本页不能输出VCF、SAM或“变异已确认”。连接[BLAST原理与种子延伸](https://github.com/SpiralKernel/bioinfo-wiki/wiki/BLAST%E5%8E%9F%E7%90%86%E4%B8%8E%E7%A7%8D%E5%AD%90%E5%BB%B6%E4%BC%B8)比较候选生成思想，连接[短读长与长读长变异证据比较](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%9F%AD%E8%AF%BB%E9%95%BF%E4%B8%8E%E9%95%BF%E8%AF%BB%E9%95%BF%E5%8F%98%E5%BC%82%E8%AF%81%E6%8D%AE%E6%AF%94%E8%BE%83)理解比对证据如何影响检测。

## 排错与自测

| 症状 | 优先检查 |
|---|---|
| 命中坐标错1 | 0/1基与Occ前缀是否包含当前位 |
| 搜索区间越来越怪 | l/r是否用旧区间同时更新，C是否严格小于 |
| 明明有子串却找不到 | 查询是否反向处理、字母表/大小写一致 |
| 大输入内存暴涨 | 是否保存全部后缀或旋转矩阵 |

**问：BWT字符顺序已变，原文是否丢失？** 答：有唯一终止符及对应约定时可逆。**问：r−l=2代表两种变异吗？** 答：只代表两个精确文本命中，不能作生物事件计数。**问：删除终止符再随便换排序能否仍用同公式？** 答：不能，唯一哨兵和排序定义属于算法条件。

本页Python代码在项目根目录执行，仅用标准库且只打印结果，不生成文件；Windows用`python`，macOS/Linux可用`python3`。学习顺序见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。


## 运行本页的完整小例

先下载项目并进入含README与examples的根目录，环境准备见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。用项目环境的Python运行以下命令（macOS/Linux路径示例；Windows用`.venv\Scripts\python.exe`）：

```bash
.venv/bin/python examples/methods_foundations_demo.py --case fm-index --output examples/output/methods-fm-index
```

输出summary.json保存本页计算与检查结果。目录非空时拒绝覆盖，再练习请更换`--output`目录；页面内短代码只打印结果。标准库case也可用系统Python直接运行。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
