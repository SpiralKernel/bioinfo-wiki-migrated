# 置换检验与Bootstrap实践

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26核对原始与官方资料；本页Python微型计算已执行；人工数据不代表专业软件或真实研究性能。

置换检验通过在零假设允许的方式下重排数据，生成检验统计量的参考分布；Bootstrap通过从观测单位有放回抽样，估计统计量的抽样不确定性。一个主要服务检验，一个主要服务估计精度，两者不能只因都“重复抽样”就混用。

前置为[假设检验与双侧单侧选择](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%81%87%E8%AE%BE%E6%A3%80%E9%AA%8C%E4%B8%8E%E5%8F%8C%E4%BE%A7%E5%8D%95%E4%BE%A7%E9%80%89%E6%8B%A9)和[置信区间与重采样](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%BD%AE%E4%BF%A1%E5%8C%BA%E9%97%B4%E4%B8%8E%E9%87%8D%E9%87%87%E6%A0%B7)。读完应能在六个数上枚举全部独立组标签分配，并解释为什么配对/家系不能直接乱洗。

## 先问哪些单位可以交换

独立两组置换通常要求H0下标签可交换，例如来自相同分布或随机分配机制。检验仅均值相等但方差不等时，简单原始均值差置换未必精确，应考虑合适的学生化统计量或模型。配对设计可在配对内交换条件；相关性检验可重配一组与另一组的对应；它们检验不同H0。[SciPy置换说明](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html)区分这些重排机制。

细胞嵌套在供体、样本嵌套家系或时间重复测量时，随意打乱所有行会破坏相关结构。应按设计使用供体/家系块、受限置换或适合模型的策略；不是“置换次数很多”就能补救错误单位。

## 六个数枚举20种分组

A=[1,2,3]，B=[4,5,6]，目标B−A=3。总共C(6,3)=20种保留每组3个对象的标签分配。只有A取最小三个或最大三个能达到绝对差至少3，因此绝对统计量双侧P=2/20=0.1。

> 本页例子在项目根目录的Python中运行；标准库例子无需安装包，含NumPy/SciPy的例子先按[R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83)建立项目环境。Windows命令使用`python`，macOS/Linux可使用`python3`；代码块可复制进该环境的Python交互会话，或保存为项目根目录下的临时脚本再用该环境解释器运行；不要保存到examples目录或任意外部目录后假定项目导入路径仍相同。除特别说明外，只向屏幕输出，不写文件。

```python
from itertools import combinations
values = [1.,2.,3.,4.,5.,6.]
observed = 3.
null = []
for Aidx in combinations(range(6),3):
    selected = set(Aidx)
    A = [x for i,x in enumerate(values) if i in selected]
    B = [x for i,x in enumerate(values) if i not in selected]
    null.append(sum(B)/3-sum(A)/3)
p = sum(abs(t)>=abs(observed)-1e-12 for t in null)/len(null)
print('allocations_exact_p', len(null), p)
assert len(null)==20 and p==.1
```

这是完整枚举，分母20包含原分配。Monte Carlo只随机抽B次时，常用(b+1)/(B+1)等包含观测的保守计算，b为至少同样极端次数，不能简单报告0/B为P=0。不同双侧定义可导致不同P，报告“绝对差双侧”或软件定义。

## Bootstrap是另一种抽样世界

对n个独立对象有放回抽n个，计算一次均值/差值，重复形成Bootstrap分布。若同一对象被抽中两次，它应作为两个抽样副本保留；先去重会改变分布。配对Bootstrap抽完整对，cluster Bootstrap抽整个独立簇。配对拆散就不再保留原相关。

百分位区间取重采样统计量的2.5%与97.5%分位数，直观但未必覆盖良好；basic和BCa采用不同校正。小样本、极端值、参数边界或离散退化下可能失败。见[Bootstrap官方文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html)。

```python
from random import Random
from statistics import mean
x, rng = [1.,2.,3.,4.], Random(20261026)
boot = [mean(rng.choices(x,k=len(x))) for _ in range(2000)]
print('bootstrap_draws', len(boot), 'mean_near', round(mean(boot),3))
print('range', min(boot), max(boot))
assert len(boot)==2000 and all(1<=v<=4 for v in boot)
```

均值接近2.5，所有重采样均值位于[1,4]；这只是人工抽样演示，未用这四点宣称可靠的95%区间。固定种子便于核验，却不保证不同随机库版本逐项一致。

## 怎样决定重复次数与报告

精确枚举可行时避免Monte Carlo误差；分组数很大时采用随机抽样，并记录种子、B、极端计数和双侧规则。小P附近需要更多抽样才能分辨，最小可报告尺度约为1/(B+1)。Bootstrap关注分位数稳定性和覆盖假设，不能只因重复一万次就称“精确置信区间”。

| 症状 | 优先核查 |
|---|---|
| 换种子P明显跨阈值 | Monte Carlo误差与次数，勿只留最小值 |
| 区间全部相同 | 原数据退化、统计量错误或样本太少 |
| Bootstrap方差极小 | 是否抽了技术重复或错误地去重 |
| 受限置换结果不同 | 参考零假设是否发生改变 |

**自测：把10个供体的10000个细胞逐细胞Bootstrap，会得到供体层面准确区间吗？** 答：通常不行，独立单位是供体时应保留该层级。**自测：Bootstrap生成的数据代表新的真实生物个体吗？** 答：不是，它是基于经验分布的重采样近似。更多已执行实例见[合成数据贯通统计推断](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%90%88%E6%88%90%E6%95%B0%E6%8D%AE%E8%B4%AF%E9%80%9A%E7%BB%9F%E8%AE%A1%E6%8E%A8%E6%96%AD)。


## 接着读什么

[假设检验与双侧单侧选择](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%81%87%E8%AE%BE%E6%A3%80%E9%AA%8C%E4%B8%8E%E5%8F%8C%E4%BE%A7%E5%8D%95%E4%BE%A7%E9%80%89%E6%8B%A9) → [置信区间与重采样](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%BD%AE%E4%BF%A1%E5%8C%BA%E9%97%B4%E4%B8%8E%E9%87%8D%E9%87%87%E6%A0%B7) → [合成数据贯通统计推断](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%90%88%E6%88%90%E6%95%B0%E6%8D%AE%E8%B4%AF%E9%80%9A%E7%BB%9F%E8%AE%A1%E6%8E%A8%E6%96%AD) → [R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83) → [方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。学习顺序与练习入口见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。


## 运行本页的完整小例

先下载项目并进入含README与examples的根目录，环境准备见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。用项目环境的Python运行以下命令（macOS/Linux路径示例；Windows用`.venv\Scripts\python.exe`）：

```bash
.venv/bin/python examples/methods_foundations_demo.py --case resampling --output examples/output/methods-resampling
```

输出summary.json保存本页计算与检查结果。目录非空时拒绝覆盖，再练习请更换`--output`目录；页面内短代码只打印结果。标准库case也可用系统Python直接运行。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
