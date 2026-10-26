# GBLUP与rrBLUP原理及实现

> 更新：2026-10-26；状态：可用；验证范围：官方/原始来源核对；标准库教学算例已实跑；专业流程范围见正文。

GBLUP按个体之间的关系预测遗传值，RR-BLUP按许多标记的小效应预测。两者看起来使用不同规模的矩阵，却在一组明确的线性高斯假设下给出相同预测。本页要能检查这种等价的条件、把候选个体接到训练基准，并区分遗传预测与完整表型预测。

## 输入和中心化先于求解

设训练基因型M为n×m剂量矩阵，列j始终计数同一种等位基因；p_j是选定训练或基础群体的频率。定义 $Z_{ij}=M_{ij}-2p_j$，$s=\sum_j2p_j(1-p_j)$，一种常用关系矩阵为 $K=ZZ^T/s$。K也常记G，此处不用Z表示随机效应的关联矩阵，以免混淆。频率基准及缩放会影响遗传方差单位。[VanRaden原始方法](https://pubmed.ncbi.nlm.nih.gov/18946147/)

缺失填补、剔除列和中心化规则要连同标记ID保存。候选Z*必须沿用相同频率、计数方向、列顺序和缺失规则，不能在每批候选重新中心化。若使用允许候选基因型参与的参考设计，应明确其部署条件，不与训练限定的评估混称。

## 个体模型及预测

写 $y=Xb+u+e$，其中y是n维训练表型，X是固定效应设计，$u\sim N(0,K\sigma_u^2)$，$e\sim N(0,I\sigma_e^2)$。令 $\lambda_G=\sigma_e^2/\sigma_u^2$、$A=K+\lambda_G I$。先求解线性系统而非直接求逆：

$$
\hat b=(X^TA^{-1}X)^{-1}X^TA^{-1}y,\quad
\alpha=A^{-1}(y-X\hat b),\quad \hat u_*=K_{*,t}\alpha.
$$

K*,t是候选对训练的交叉关系块。预测完整表型还需加X* b；育种值与表型不是同一输出。正lambda使K即使只有半正定也可得到可解系统，X仍须列满秩；多环境、多性状或异方差残差需要扩展模型。

## 标记模型及尺度关系

若 $u=Za$，且 $a_j$独立同分布为零均值高斯，则 $Var(u)=ZZ^T\sigma_a^2$。要与上面K匹配，必须有 $\sigma_a^2=\sigma_u^2/s$，所以标记惩罚 $\lambda_a=s\lambda_G$。对应方程为：

$$
\hat a=(Z^TZ+s\lambda_GI)^{-1}Z^T(y-X\hat b),\quad \hat u_*=Z_*\hat a.
$$

同一数据、固定效应、中心化、权重与方差尺度下才等价。两个软件都叫lambda却一个表示方差比、一个还除以n，直接填同一数通常不是同一模型。rrBLUP的求解器同时提供ML/REML估计与预测接口，但教学脚本固定lambda的过程并没有执行REML。[rrBLUP原作者接口说明](https://rdrr.io/cran/rrBLUP/man/mixed.solve.html)

## 两个个体的可核算例

设硬分型M两行为(0,0)、(2,1)，训练频率为(0.5,0.25)，则中心化矩阵Z两行为(−1,−0.5)、(1,0.5)，s=0.875，表型取−1、1。Z中的半整数来自减去2p，不是观察到了半份基因拷贝。K对角为10/7、非对角为−10/7。取lambda_G=1、截距0，对比方向上的特征值20/7被收缩为(20/7)/(20/7+1)=20/27。因此两条遗传预测为±0.740741。

标记路径应使用lambda_a=0.875，可得a约(0.592593,0.296296)，代回Z也得到同一结果。这个例子直接暴露“关系矩阵缩放变了但惩罚没变”的错误；不适合估计真实群体的遗传力或准确度。

## 模型大小、参数与交付

当m远大于n时，个体方程可能更省；当n远大于m时，标记方程可能更合适。生产程序还会利用矩阵分解、迭代和稀疏性，不能以本例消元器推断百万样本资源需求。lambda可由REML或训练内验证选择，两者优化目标不同；最终比较必须固定外层划分。

[合成家系数据贯通GS](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%90%88%E6%88%90%E5%AE%B6%E7%B3%BB%E6%95%B0%E6%8D%AE%E8%B4%AF%E9%80%9AGS)及[既有NumPy实现](https://github.com/SpiralKernel/bioinfo-wiki/blob/main/examples/genomic_prediction_demo.py)已经验证个体/标记预测等价到浮点误差级。交付还需训练频率、保留标记、固定效应编码、参数、候选交叉块和验证设计；只有一列GEBV不能复现下一批预测。

## 动手：运行后先对照什么

在下载项目的根目录执行；只需Python 3标准库，不安装专业统计包：

```bash
python3 examples/methods_genetics_demo.py --case gblup
```

输出`scale=0.875`（公式中的s）、`lambda_marker=0.875`，个体与标记路径均约`[-0.740741,0.740741]`，两者差小于1e-12。手工将a代回两行Z核对符号。

这是已知方差比、截距为0的二个体算例；没有REML、PEV或真实数据性能评估。小型消元器仅为教学，不用于大型遗传评估。

命令只向终端打印JSON。若要保存所有22个小例，可用`python3 examples/methods_genetics_demo.py --output examples/output/methods-genetics-practice1`；目录必须新建或为空，已有结果不会被覆盖。源代码见[标准库方法算例](https://github.com/SpiralKernel/bioinfo-wiki/blob/main/examples/methods_genetics_demo.py)。

## 排错：先定位原因再改参数

| 症状 | 可能原因 | 下一步检查 |
|---|---|---|
| 两条路径不相同 | lambda缩放、中心或固定效应不同 | 比较s、lambda_a/s、X和候选标记顺序 |
| 候选批次均值漂移 | 候选重新计算了频率基准 | 将训练转换器原样应用到候选 |
| 关系矩阵有负非对角 | 相对基准下的中心化关系 | 先查定义；不因负值就截0 |

## 自测与答案

1. **s减半而lambda_G固定时，标记lambda应如何变化？** 也减半；否则标记与个体模型不再匹配。

2. **GEBV为何不能直接当预测表型？** 完整表型还含固定效应及可能的其他随机效应，且实际表型含不可预测残差。

## 原有实践与延伸入口

[genomic_prediction_demo.py](https://github.com/SpiralKernel/bioinfo-wiki/blob/main/examples/genomic_prediction_demo.py)；[VanRaden方法](https://pubmed.ncbi.nlm.nih.gov/18946147/)；[rrBLUP预测实践](https://github.com/SpiralKernel/bioinfo-wiki/wiki/rrBLUP%E9%A2%84%E6%B5%8B%E5%AE%9E%E8%B7%B5)；[GS训练验证与数据泄漏](https://github.com/SpiralKernel/bioinfo-wiki/wiki/GS%E8%AE%AD%E7%BB%83%E9%AA%8C%E8%AF%81%E4%B8%8E%E6%95%B0%E6%8D%AE%E6%B3%84%E6%BC%8F)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
