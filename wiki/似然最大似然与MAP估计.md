# 似然最大似然与MAP估计

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26核对原始与官方资料；本页Python微型计算已执行；人工数据不代表专业软件或真实研究性能。

似然问：“如果参数取某个值，观察到当前数据有多符合模型？”最大似然（MLE）选择使该量最大的参数；最大后验（MAP）还纳入先验。读完应能用抛硬币例算出二者差别，并理解似然不是参数本身的概率。

前置为[概率条件概率与贝叶斯定理](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%A6%82%E7%8E%87%E6%9D%A1%E4%BB%B6%E6%A6%82%E7%8E%87%E4%B8%8E%E8%B4%9D%E5%8F%B6%E6%96%AF%E5%AE%9A%E7%90%86)与[常见概率分布与生信场景](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%B8%B8%E8%A7%81%E6%A6%82%E7%8E%87%E5%88%86%E5%B8%83%E4%B8%8E%E7%94%9F%E4%BF%A1%E5%9C%BA%E6%99%AF)。实际优化可接[梯度下降与数值优化](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%A2%AF%E5%BA%A6%E4%B8%8B%E9%99%8D%E4%B8%8E%E6%95%B0%E5%80%BC%E4%BC%98%E5%8C%96)。

## 把观测固定，把参数当变量

若n次独立Bernoulli观测中有k次成功，成功率为p，忽略与p无关的组合常数，$L(p)=p^k(1-p)^{n-k}$。观测已经固定，p在0到1间变化。它不是关于p归一化的概率分布；给定p的观测概率与给定观测的参数概率不能倒置。

取对数得$\ell(p)=k\log p+(n-k)\log(1-p)$。对数单调递增，所以最大点不变，又能把乘法变加法，减少极小概率连乘下溢。内部极值由k/p−(n−k)/(1−p)=0得p̂=k/n；k=0或k=n时最大点在边界。

## 加入先验，为什么估计向中间移动

设p服从Beta(a,b)先验，密度正比$p^{a-1}(1-p)^{b-1}$，后验为Beta(k+a,n−k+b)。若后验两个形状参数均大于1，后验众数即MAP为$(k+a-1)/(n+a+b-2)$，后验均值却是$(k+a)/(n+a+b)$；两者不是同一量。边界参数不满足内部条件时需单独判断众数，不能机械代公式。

n=10、k=8：MLE=0.8。Beta(2,2)给MAP=9/12=0.75，后验均值10/14≈0.7143。这体现先验影响；样本很少或似然很平时，影响通常更明显。基础解释可参考[Harvard Stat110课程资料](https://stat110.hsites.harvard.edu/)。

> 本页例子在项目根目录的Python中运行；标准库例子无需安装包，含NumPy/SciPy的例子先按[R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83)建立项目环境。Windows命令使用`python`，macOS/Linux可使用`python3`；代码块可复制进该环境的Python交互会话，或保存为项目根目录下的临时脚本再用该环境解释器运行；不要保存到examples目录或任意外部目录后假定项目导入路径仍相同。除特别说明外，只向屏幕输出，不写文件。

```python
from math import log, log1p
n, k, a, b = 10, 8, 2, 2
mle = k/n
map_estimate = (k+a-1)/(n+a+b-2)
posterior_mean = (k+a)/(n+a+b)
def ll(p):
    return k*log(p)+(n-k)*log1p(-p)
print('MLE_MAP_mean', round(mle,6), round(map_estimate,6), round(posterior_mean,6))
print('likelihood_prefers_MLE', ll(mle) > ll(map_estimate))
assert abs(map_estimate-.75) < 1e-12
assert ll(mle) > ll(.5)
```

检查点为0.8、0.75、0.714286。似然在MLE处较高并不与MAP矛盾，因为MAP最大化的是“对数似然+对数先验”，目标函数不同。本例没有推断任何真实群体的等位频率。

## MAP与正则化怎样连接

在线性高斯误差模型中，加上以0为中心的高斯系数先验，负对数后验含平方系数惩罚，便与岭回归相联系；Laplace先验与L1惩罚有类似联系。比例因子取决于损失定义和方差参数。惩罚强度大意味着更强收缩，但不说明该先验适合所有性状，见[正则化过拟合与模型选择](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%AD%A3%E5%88%99%E5%8C%96%E8%BF%87%E6%8B%9F%E5%90%88%E4%B8%8E%E6%A8%A1%E5%9E%8B%E9%80%89%E6%8B%A9)与[贝叶斯回归与核方法用于GS](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%B4%9D%E5%8F%B6%E6%96%AF%E5%9B%9E%E5%BD%92%E4%B8%8E%E6%A0%B8%E6%96%B9%E6%B3%95%E7%94%A8%E4%BA%8EGS)。

MAP只给一个点；MCMC或其他后验近似则试图描述完整分布。报告MAP不等于已经报告后验不确定性，也不能把优化器的Hessian标准误不加说明地当作任何模型下的精确可信区间。

## 生信里的似然链条

基因型似然P(reads|GT)依据比对、质量和错误模型评价候选基因型；加入基因型先验才得到后验。混合模型估计方差分量、RNA计数模型拟合、祖源模型和共祖推断也各有特定似然。表格形状相同不代表其独立性和分布假设相同。

实际拟合需检查参数可识别性：不同参数能否产生相同预测？样本不足、共线、完全分离或边界方差都可能使极值和区间不稳定。模型错设时，最大化似然仍会得到一个“最符合错模型”的值，不能当作真实机制。

## 异常与自测

| 症状 | 可能原因 | 核查 |
|---|---|---|
| 概率连乘变0 | 浮点下溢 | 累加log概率、稳定log-sum-exp |
| 参数跑到边界 | 数据极端或约束/识别问题 | 看似然曲线与有效样本 |
| 换先验估计明显变 | 信息少或先验强 | 报告敏感性，勿隐藏 |
| 优化成功但SE巨大 | 平坦方向、共线 | 检查设计秩、Hessian和数据设计 |

**问：MLE=0.8表示p有80%概率等于0.8吗？** 答：不是。**问：MAP在参数重表达后一定简单变换吗？** 答：后验密度涉及Jacobian，众数通常不具备简单不变性；参数化和先验必须明确。MLE则对一一参数变换的极值对应有相应不变性。


## 接着读什么

[概率条件概率与贝叶斯定理](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%A6%82%E7%8E%87%E6%9D%A1%E4%BB%B6%E6%A6%82%E7%8E%87%E4%B8%8E%E8%B4%9D%E5%8F%B6%E6%96%AF%E5%AE%9A%E7%90%86) → [常见概率分布与生信场景](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%B8%B8%E8%A7%81%E6%A6%82%E7%8E%87%E5%88%86%E5%B8%83%E4%B8%8E%E7%94%9F%E4%BF%A1%E5%9C%BA%E6%99%AF) → [梯度下降与数值优化](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%A2%AF%E5%BA%A6%E4%B8%8B%E9%99%8D%E4%B8%8E%E6%95%B0%E5%80%BC%E4%BC%98%E5%8C%96) → [正则化过拟合与模型选择](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%AD%A3%E5%88%99%E5%8C%96%E8%BF%87%E6%8B%9F%E5%90%88%E4%B8%8E%E6%A8%A1%E5%9E%8B%E9%80%89%E6%8B%A9) → [贝叶斯回归与核方法用于GS](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%B4%9D%E5%8F%B6%E6%96%AF%E5%9B%9E%E5%BD%92%E4%B8%8E%E6%A0%B8%E6%96%B9%E6%B3%95%E7%94%A8%E4%BA%8EGS) → [R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83) → [GS中的REML混合模型方程与求解](https://github.com/SpiralKernel/bioinfo-wiki/wiki/GS%E4%B8%AD%E7%9A%84REML%E6%B7%B7%E5%90%88%E6%A8%A1%E5%9E%8B%E6%96%B9%E7%A8%8B%E4%B8%8E%E6%B1%82%E8%A7%A3) → [方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。学习顺序与练习入口见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。


## 运行本页的完整小例

先下载项目并进入含README与examples的根目录，环境准备见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。用项目环境的Python运行以下命令（macOS/Linux路径示例；Windows用`.venv\Scripts\python.exe`）：

```bash
.venv/bin/python examples/methods_foundations_demo.py --case mle-map --output examples/output/methods-mle-map
```

输出summary.json保存本页计算与检查结果。目录非空时拒绝覆盖，再练习请更换`--output`目录；页面内短代码只打印结果。标准库case也可用系统Python直接运行。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
