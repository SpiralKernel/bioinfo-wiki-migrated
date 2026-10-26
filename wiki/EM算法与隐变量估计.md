# EM算法与隐变量估计

> 更新：2026-10-26；状态：可用；验证范围：官方/原始资料核对；本页标准库Python人工例已执行；未运行生产软件或真实生物数据。

期望最大化（expectation-maximization，EM）用于带有未观测变量的模型：若知道每个样本来自哪个成分，参数就容易估计；但成分恰好未知。EM交替估计软归属和更新参数。本页用两个成分的混合模型完成一次可复算迭代，避免把EM误解成“先强行分类再取均值”。

前置为[似然最大似然与MAP估计](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E4%BC%BC%E7%84%B6%E6%9C%80%E5%A4%A7%E4%BC%BC%E7%84%B6%E4%B8%8EMAP%E4%BC%B0%E8%AE%A1)、条件概率与加权均值。读完应能解释E步算什么、M步优化什么，以及似然不下降为什么不代表找到全局最优。

## 隐变量让似然出现求和

观测xi来自K个成分，隐变量zi标记来源。混合权重πk满足非负且和为1，成分密度为fk(xi|θk)，观测似然为$\prod_i\sum_k\pi_k f_k(x_i|\theta_k)$。log作用在求和外面，难以像已知分类一样直接分别求参数。

E步以当前参数计算责任度$r_{ik}=P(z_i=k|x_i)=\pi_kf_k(x_i)/\sum_h\pi_hf_h(x_i)$。每个样本对各成分责任度和为1。M步最大化按这些责任度加权的完整数据对数似然；高斯混合的一部分更新为$N_k=\sum_i r_{ik}$、$\pi_k=N_k/n$、$\mu_k=\sum_i r_{ik}x_i/N_k$。

“E”不是直接对原始数据取期望，而是针对当前参数下隐变量的条件分布。原始算法见[Dempster、Laird与Rubin 1977](https://doi.org/10.1111/j.2517-6161.1977.tb01600.x)；[scikit-learn混合模型说明](https://scikit-learn.org/stable/modules/mixture.html)给出高斯混合的实际建模条件与收敛问题。

## 两成分的一次更新

观测x=[0,1,4,5]，两个高斯成分方差都固定1，初始均值0与5、权重各0.5。靠近0的样本主要归属左成分，靠近5的主要归属右成分，但责任度并非硬0/1。E步之后，M步将两均值向各自两个点的中心移动到约0.5与4.5，权重保持0.5。

```python
from math import exp, log, pi
x, means, weights = [0.,1.,4.,5.], [0.,5.], [.5,.5]
def density(v,m):
    return exp(-.5*(v-m)**2)/(2*pi)**.5
def likelihood(ms,ws):
    return sum(log(sum(w*density(v,m) for w,m in zip(ws,ms))) for v in x)
responsibilities = []
for v in x:
    a = [w*density(v,m) for w,m in zip(weights,means)]
    responsibilities.append([z/sum(a) for z in a])
Nk = [sum(row[k] for row in responsibilities) for k in range(2)]
new_weights = [n/len(x) for n in Nk]
new_means = [sum(row[k]*v for row,v in zip(responsibilities,x))/Nk[k] for k in range(2)]
print('left_responsibilities',[round(row[0],6) for row in responsibilities])
print('new_means',[round(m,6) for m in new_means],'weights',new_weights)
print('loglik_before_after',round(likelihood(means,weights),6),round(likelihood(new_means,new_weights),6))
assert all(abs(sum(row)-1)<1e-12 for row in responsibilities)
assert likelihood(new_means,new_weights) >= likelihood(means,weights)-1e-12
```

检查点为左侧责任度接近1、1、0、0，两个新均值接近0.5和4.5；对数似然提高。这里方差固定，只更新均值与权重，是明确受约束的高斯混合EM，不是完整GaussianMixture拟合。

## 为什么需要初始化与多个起点

在适当条件下，每个精确EM更新不会降低观测似然，但目标可能有多个局部极值或退化点。两个均值都初始化在同一点，算法可能一直保持对称，无法自动分开成分；似然平台期也可能收敛很慢。换标签顺序产生相同模型，不应误判为不同生物群体。

完整高斯混合还估计协方差；某个成分追逐单个样本、方差趋0时可能产生奇异解，需要正则化、约束和诊断。K越大通常更能提高训练似然，但不保证更合理，需按任务采用BIC、独立验证和可解释性等依据。

## 与K-means、梯度下降的区别

K-means用硬成员分配并最小化平方距离；EM用成分模型的软责任度，成分尺度与权重参与概率。二者在特殊限制和极限下有关联，但不是同一算法。EM也不等于梯度下降：M步有时可解析更新，有时内部需数值优化，取决于模型。

基因表达量估计、祖源和隐状态模型中的隐变量具有不同含义；把混合模型成分直接命名为“细胞类型”或“真实祖先”需要额外注释与验证。专业模型的偏差、约束与可识别性不能由这个小例代替。

## 排错与自测

| 现象 | 先查什么 |
|---|---|
| 某轮似然明显下降 | E/M步参数是否混用、数值下溢、M步是否真正优化 |
| 所有样本同一成分 | 初始化、K、数据重叠、权重退化 |
| 责任度出现NaN | 所有密度下溢为0，需log-sum-exp |
| 达到最大迭代数 | 目标增量和参数变化，报告未收敛 |

**问：责任度0.9表示样本已被独立证实属于该生物类吗？** 答：不是，它是在当前模型和参数下的条件分配概率。**问：每轮似然提高就保证K选对了吗？** 答：不保证。**下一步**：[隐马尔可夫模型与Viterbi解码](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E9%9A%90%E9%A9%AC%E5%B0%94%E5%8F%AF%E5%A4%AB%E6%A8%A1%E5%9E%8B%E4%B8%8EViterbi%E8%A7%A3%E7%A0%81)说明当隐变量沿序列相关时，需要怎样的动态规划。

本页Python代码在项目根目录执行，仅用标准库且只打印结果，不生成文件；Windows用`python`，macOS/Linux可用`python3`。学习顺序见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。


## 运行本页的完整小例

先下载项目并进入含README与examples的根目录，环境准备见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。用项目环境的Python运行以下命令（macOS/Linux路径示例；Windows用`.venv\Scripts\python.exe`）：

```bash
.venv/bin/python examples/methods_foundations_demo.py --case em --output examples/output/methods-em
```

输出summary.json保存本页计算与检查结果。目录非空时拒绝覆盖，再练习请更换`--output`目录；页面内短代码只打印结果。标准库case也可用系统Python直接运行。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
