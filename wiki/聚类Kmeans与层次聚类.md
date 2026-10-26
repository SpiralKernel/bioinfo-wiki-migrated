# 聚类Kmeans与层次聚类

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26核对原始与官方资料；本页Python微型计算已执行；人工数据不代表专业软件或真实研究性能。

聚类在没有已知类别标签时，按特征表示和相似规则把对象分组。K-means寻找围绕中心的紧凑簇，层次聚类产生逐步合并的树。读完应能手工完成一次K-means迭代，并说明“得到两个簇”为什么不等于发现两个真实生物群体。

前置为[距离标准化与相似性](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%B7%9D%E7%A6%BB%E6%A0%87%E5%87%86%E5%8C%96%E4%B8%8E%E7%9B%B8%E4%BC%BC%E6%80%A7)。本页的二维点和一维手算均为人工数据，结果只检验目标函数与算法步骤。

## K-means在优化什么

给n个p维点xi与K个中心μk，最小化簇内平方欧氏距离和：$J=\sum_i\|x_i-\mu_{c_i}\|^2$，ci为每点所属簇。Lloyd算法交替两步：固定中心，把每点分到最近中心；固定成员，用成员均值更新中心。每一步不会增加该目标，但通常只保证达到局部解或稳定分配，不保证全局最优。初始化、K和停止条件均影响结果。[scikit-learn聚类指南](https://scikit-learn.org/stable/modules/clustering.html)比较这些假设与其他聚类方法。

为何中心是均值？在成员固定时，一维平方偏差和对μ求导得2∑(μ−xi)=0，故μ=均值；各维独立成立。若换成绝对距离，最优代表是中位数，已不是相同的K-means目标。

## 四个点做一次更新

点为0、1、9、10，K=2，初始中心0与10。分配得到{0,1}与{9,10}；更新后中心0.5和9.5，平方误差为4×0.25=1。再次分配不变则收敛。若把K设4，每点自己成簇，误差降到0，但这不说明四簇一定更合理。

> 本页例子在项目根目录的Python中运行；标准库例子无需安装包，含NumPy/SciPy的例子先按[R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83)建立项目环境。Windows命令使用`python`，macOS/Linux可使用`python3`；代码块可复制进该环境的Python交互会话，或保存为项目根目录下的临时脚本再用该环境解释器运行；不要保存到examples目录或任意外部目录后假定项目导入路径仍相同。除特别说明外，只向屏幕输出，不写文件。

```python
points, centers = [0.,1.,9.,10.], [0.,10.]
labels = [min(range(2), key=lambda k:(x-centers[k])**2) for x in points]
updated = [sum(x for x,c in zip(points,labels) if c==k)/labels.count(k)
           for k in range(2)]
sse = sum((x-updated[c])**2 for x,c in zip(points,labels))
print('labels', labels, 'centers', updated, 'SSE', sse)
assert labels == [0,0,1,1] and sse == 1
from examples.basic_algorithms_demo import kmeans
import numpy as np
C, labels2, history = kmeans(np.array(points)[:,None],2)
print('final_SSE', history[-1], 'nonincreasing', bool((np.diff(history)<=1e-10).all()))
assert np.isclose(history[-1],1)
```

第一个部分标准库即可，第二部分调用依赖NumPy的项目小实现。簇编号0和1可以交换，应比较中心集合、成员关系或目标值，不能把标签编号当稳定生物名称。空簇的处理也属于实现选择；本项目保留空簇旧中心，不冒充所有软件行为。

## 层次聚类为什么依赖linkage

凝聚式层次聚类从每点一簇开始，每次合并最近的两簇。single linkage看最近跨簇点，容易链式连接；complete看最远点，更偏紧凑簇；average看平均距离；Ward按合并带来的簇内平方误差增量选择，要求与欧氏几何相容。相同距离输入换linkage会得到不同树，并非某一棵必然错误。

树状图高度表达该算法定义的合并距离/代价，不能直接当成进化时间。叶片左右旋转通常不改变聚类树，图上视觉相邻也不一定代表最近。截树高度和K应记录，不能只保存彩色热图。

## 怎样评价与选择参数

K-means偏好相对凸、近似球状且尺度可比的簇；弯月形、密度差异大或离群点很多时可能不适合。K越大，训练SSE通常越低；“肘部”是启发判断，不是显著性检验。轮廓系数、重采样稳定性与生物注释帮助评估，但同一数据上多次试参再挑最漂亮分组存在选择偏差。

单细胞邻居图聚类常使用图社群目标，不应直接解释为K-means；bulk表达热图又常先按基因缩放，需确认正在聚类样本还是基因。样本数很大时完整两两距离和层次树可能需要O(n²)资源，不能直接沿用小热图做法。

## 异常诊断与自测

| 症状 | 检查 | 可采取的动作 |
|---|---|---|
| 每次结果都变 | 随机种子、初始化、弱分离 | 多次起点，报告稳定性 |
| 一点单独成簇 | 离群、尺度、K过大 | 查QC与影响，勿直接删除 |
| 聚类完全按批次 | 预处理与混杂 | 先诊断批次，不重命名为细胞类型 |
| SSE增加 | 记录时点、分配/更新实现 | 固定小例逐步重算 |

**问：SSE从2降到1说明生物分群更真实吗？** 答：仅说明更好优化了此目标。**问：两次标签恰好反过来是否不稳定？** 答：标签置换不影响成员分组，需先匹配标签再比较。


## 接着读什么

[距离标准化与相似性](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%B7%9D%E7%A6%BB%E6%A0%87%E5%87%86%E5%8C%96%E4%B8%8E%E7%9B%B8%E4%BC%BC%E6%80%A7) → [R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83) → [梯度下降与数值优化](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%A2%AF%E5%BA%A6%E4%B8%8B%E9%99%8D%E4%B8%8E%E6%95%B0%E5%80%BC%E4%BC%98%E5%8C%96) → [基础算法可复现练习](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E5%9F%BA%E7%A1%80%E7%AE%97%E6%B3%95%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%BB%83%E4%B9%A0) → [方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。学习顺序与练习入口见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。


## 运行本页的完整小例

先下载项目并进入含README与examples的根目录，环境准备见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。用项目环境的Python运行以下命令（macOS/Linux路径示例；Windows用`.venv\Scripts\python.exe`）：

```bash
.venv/bin/python examples/methods_foundations_demo.py --case clustering --output examples/output/methods-clustering
```

输出summary.json保存本页计算与检查结果。目录非空时拒绝覆盖，再练习请更换`--output`目录；页面内短代码只打印结果。标准库case也可用系统Python直接运行。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
