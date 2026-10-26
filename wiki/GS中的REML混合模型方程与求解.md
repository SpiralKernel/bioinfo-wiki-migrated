# GS中的REML混合模型方程与求解

> 更新：2026-10-26；状态：可用；验证范围：官方/原始来源核对；标准库教学算例已实跑；专业流程范围见正文。

遗传评估软件中的“估方差”“求育种值”“算不确定性”是三个相连但不同的任务。REML估计协方差参数，混合模型方程在给定这些参数时求BLUE/BLUP，求解器负责数值计算。本页以最小例解释自由度，再把它接回多随机效应模型；不把一次收敛日志当作生物学正确的证明。

前置为[线性模型与混合模型](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%BA%BF%E6%80%A7%E6%A8%A1%E5%9E%8B%E4%B8%8E%E6%B7%B7%E5%90%88%E6%A8%A1%E5%9E%8B)和矩阵乘法。输入包括表型、X固定效应、Z观测到随机效应的关联矩阵、关系结构、残差权重和ID映射。这里Z与[GBLUP与rrBLUP原理及实现](https://github.com/SpiralKernel/bioinfo-wiki/wiki/GBLUP%E4%B8%8ErrBLUP%E5%8E%9F%E7%90%86%E5%8F%8A%E5%AE%9E%E7%8E%B0)中的中心化标记矩阵记号不同，必须按上下文读维度。

## 为什么不是先减均值再普通ML

四个独立高斯观测1、2、3、4，未知共同均值的估计为2.5，残差平方和为5。最大似然残差方差是5/4=1.25；在仅截距模型下，REML得到5/(4−1)=1.666667，因为估计均值用掉了一个自由度。这里展示的是单一残差方差；复杂混合模型不是把每个方差都机械除以n−p。

REML可理解为基于消除了固定效应的误差对比建立似然。这样把估计固定效应带来的信息损失纳入协方差估计。比较含不同固定效应的模型时，REML似然对应的误差对比空间不同，不能不加处理地用其数值做普通似然比；常需用ML或相应方法。[REML原始论文](https://doi.org/10.1093/biomet/58.3.545)

## 模型、协方差和限制似然

设 $y=Xb+Zu+e$，$Var(u)=G\sigma_u^2$，$Var(e)=R\sigma_e^2$，两者独立，则 $V=ZGZ^T\sigma_u^2+R\sigma_e^2$。G、R在此是相对结构，不要又把绝对方差重复乘进去。消去固定效应后：

$$
P=V^{-1}-V^{-1}X(X^TV^{-1}X)^{-1}X^TV^{-1},
$$
$$
\ell_R= -\tfrac12[\log|V|+\log|X^TV^{-1}X|+y^TPy]+常数.
$$

参数改变V，也改变投影P；因此REML优化和一次普通回归不同。AI-REML、EM等更新策略在速度、初值敏感性和边界处理上不同。零方差可能是边界最优，也可能是随机效应彼此难以区分；两个几乎相同的关系结构很难靠小样本分别估计。

## 给定方差后的混合模型方程

若G、R可逆，令lambda=σe²/σu²，一种残差尺度下的方程是：

$$
\begin{bmatrix}X^TR^{-1}X&X^TR^{-1}Z\\
Z^TR^{-1}X&Z^TR^{-1}Z+\lambda G^{-1}\end{bmatrix}
\begin{bmatrix}\hat b\\\hat u\end{bmatrix}=
\begin{bmatrix}X^TR^{-1}y\\Z^TR^{-1}y\end{bmatrix}.
$$

它同时求固定和随机效应，右下角惩罚反映u的协方差。半正定G、不可估计固定效应或有约束模型需采用相应分解/广义逆处理；不能看到奇异就给所有矩阵随意加常数，也不能把显式求逆当成唯一实现方式。[rrBLUP模型与求解说明](https://rdrr.io/cran/rrBLUP/man/mixed.solve.html)

## 数值正确与统计正确分别检查

Cholesky等直接分解适合相应正定系统；PCG迭代求解需要预条件与残差检查。对同一方程，充分收敛的两种求解方式应接近；停止过早是数值误差。APY则引入对基因组逆的近似，核心个体和结构影响近似本身，不能视为另一个完全等价的精确求解器。

至少保存参数初值/最终值、边界标志、梯度或似然变化、迭代上限、求解残差、模型秩及失败样本。方差标准误大、不同初值给不同分解时，应先检查可识别性与设计，而不是只增加迭代次数。批次、场次与家系过度重合，也会让不同方差来源不可分。

[合成家系数据贯通GS](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%90%88%E6%88%90%E5%AE%B6%E7%B3%BB%E6%95%B0%E6%8D%AE%E8%B4%AF%E9%80%9AGS)使用训练内验证选lambda并直接解线性方程，不是REML拟合；专业软件默认REML时不应期待不同lambda下逐位相同。执行入口见[rrBLUP预测实践](https://github.com/SpiralKernel/bioinfo-wiki/wiki/rrBLUP%E9%A2%84%E6%B5%8B%E5%AE%9E%E8%B7%B5)、[HIBLUP基因组预测实践](https://github.com/SpiralKernel/bioinfo-wiki/wiki/HIBLUP%E5%9F%BA%E5%9B%A0%E7%BB%84%E9%A2%84%E6%B5%8B%E5%AE%9E%E8%B7%B5)，生产资源还取决于稀疏结构和随机效应个数。

## 动手：运行后先对照什么

在下载项目的根目录执行；只需Python 3标准库，不安装专业统计包：

```bash
python3 examples/methods_genetics_demo.py --case reml
```

JSON应为`residual_SS=5`、`ML_variance=1.25`、`REML_variance≈1.666667`。将四个残差逐个平方相加，确认分母为什么相差1。

运行的是独立高斯截距模型的解析特例；没有优化两个方差组分、AI-REML、PCG或APY。完整混合模型公式用于读专业输出，不代表本脚本实现了它们。

命令只向终端打印JSON。若要保存所有22个小例，可用`python3 examples/methods_genetics_demo.py --output examples/output/methods-genetics-practice1`；目录必须新建或为空，已有结果不会被覆盖。源代码见[标准库方法算例](https://github.com/SpiralKernel/bioinfo-wiki/blob/main/examples/methods_genetics_demo.py)。

## 排错：先定位原因再改参数

| 症状 | 可能原因 | 下一步检查 |
|---|---|---|
| 某方差在0边界 | 信息不足、真实小方差或结构冗余 | 检查随机效应是否可区分及边界记录 |
| 迭代结束仍结果飘 | 停止阈值、初值或病态系统 | 比较梯度、似然轨迹和方程残差 |
| 模型间REML数值不宜直接比 | 固定效应设计不同 | 先确认比较目标，使用合适的ML/检验设计 |

## 自测与答案

1. **PCG收敛是否证明模型假设正确？** 不证明；它只说明指定线性系统被充分求解。

2. **为什么REML并不等于对所有方差统一用n−p作分母？** 简单自由度公式只适用于特定平衡/独立模型，复杂V需要整体限制似然优化。

## 原有实践与延伸入口

[rrBLUP手册](https://cran.r-project.org/web/packages/rrBLUP/rrBLUP.pdf)；[HIBLUP教程](https://www.hiblup.com/tutorials)；[BLUPF90说明](https://nce.ads.uga.edu/wiki/doku.php?id=readme.blupf90)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
