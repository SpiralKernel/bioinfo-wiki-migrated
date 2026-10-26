# 分类指标ROC与PR曲线

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26核对原始与官方资料；本页Python微型计算已执行；人工数据不代表专业软件或真实研究性能。

分类评价需要把真实标签、预测分数与决策阈值分开。ROC观察阳性与阴性排序，PR更直接显示预测阳性的可信比例和阳性找回率；二者都不能单独证明模型在新群体可用。读完应能从四个样本手算混淆矩阵、AP和ROC AUC，并解释类别比例的影响。

前置为[近邻朴素贝叶斯与决策树入门](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%BF%91%E9%82%BB%E6%9C%B4%E7%B4%A0%E8%B4%9D%E5%8F%B6%E6%96%AF%E4%B8%8E%E5%86%B3%E7%AD%96%E6%A0%91%E5%85%A5%E9%97%A8)和[数据泄漏与交叉验证](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%95%B0%E6%8D%AE%E6%B3%84%E6%BC%8F%E4%B8%8E%E4%BA%A4%E5%8F%89%E9%AA%8C%E8%AF%81)。先说明哪个标签是阳性，测试样本是否独立，以及分数越大是否越偏向阳性。

## 单个阈值与整条曲线

阈值t将分数≥t判阳性。真阳TP、假阳FP、真阴TN、假阴FN得到：precision=TP/(TP+FP)，recall=TP/(TP+FN)，specificity=TN/(TN+FP)，FPR=1−specificity。F1是precision与recall的调和平均，强调二者平衡却不计TN；accuracy=(TP+TN)/N在极不平衡时可能误导。

逐渐降低阈值可形成ROC的(FPR,TPR)和PR的(recall,precision)。ROC AUC可解释为随机阳性比分数随机阴性更高的概率，平分计半；PR面积有不同积分约定，average precision（AP）按recall增量加权precision，不等同于对PR点任意梯形积分。[scikit-learn评价指南](https://scikit-learn.org/stable/modules/model_evaluation.html)说明这些定义与零分母处理。

## 四个样本完整走一遍

按分数降序：0.9阳、0.8阴、0.7阳、0.1阴。阈值0.75时，TP=1、FP=1、FN=1、TN=1，precision=recall=0.5。阳性0.9胜过两个阴性，阳性0.7只胜过0.1，所以AUC=3/4=0.75。AP在第1和第3个位置累计recall，每次增加1/2，故AP=0.5×1+0.5×(2/3)=5/6。

> 本页例子在项目根目录的Python中运行；标准库例子无需安装包，含NumPy/SciPy的例子先按[R包Python包与可复现环境](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E5%8C%85Python%E5%8C%85%E4%B8%8E%E5%8F%AF%E5%A4%8D%E7%8E%B0%E7%8E%AF%E5%A2%83)建立项目环境。Windows命令使用`python`，macOS/Linux可使用`python3`；代码块可复制进该环境的Python交互会话，或保存为项目根目录下的临时脚本再用该环境解释器运行；不要保存到examples目录或任意外部目录后假定项目导入路径仍相同。除特别说明外，只向屏幕输出，不写文件。

```python
pairs = [(0.9,1),(0.8,0),(0.7,1),(0.1,0)]
t = .75
TP = sum(y==1 and s>=t for s,y in pairs)
FP = sum(y==0 and s>=t for s,y in pairs)
FN = sum(y==1 and s<t for s,y in pairs)
TN = sum(y==0 and s<t for s,y in pairs)
pos = [s for s,y in pairs if y==1]
neg = [s for s,y in pairs if y==0]
auc = sum((a>b)+.5*(a==b) for a in pos for b in neg)/(len(pos)*len(neg))
tp, ap = 0, 0.
for rank,(score,y) in enumerate(pairs,1):
    if y:
        tp += 1
        ap += (tp/rank)/len(pos)
print('TP_FP_FN_TN',TP,FP,FN,TN,'AUC',auc,'AP',round(ap,6))
assert (TP,FP,FN,TN)==(1,1,1,1) and auc==.75
assert abs(ap-5/6)<1e-12
```

此AP小实现假定已降序且无并列；真实分数并列应按同一阈值整批处理，使用明确文档的实现。没有阳性或没有阴性时，有的指标不可定义，不能用0替代后再平均。

## 不平衡为什么影响precision

假设灵敏度0.9、特异度0.9，若1000人只有10个阳性，预期TP=9，FP=99，precision=9/108≈0.0833。若阳性比例50%，同样灵敏度特异度给precision=0.9。模型的条件错误率没变，应用人群基率已经改变，所以PR和实际决策成本不可脱离部署场景。

ROC AUC也不是不平衡下万能指标：大量容易阴性可能掩盖关键低FPR区域，且整体AUC不说明所需阈值是否达到要求。报告阈值处指标、有效分母和置信区间，必要时评价特定FPR或特定召回率区域。

## 分数排序与概率校准

AUC高只说明排序较好，不保证分数0.8的对象有80%阳性。校准曲线或Brier分数研究概率是否匹配频率。阈值应由训练/验证数据和错误代价确定，最终测试集不可反复调阈值；交叉验证的折也需按供体、家系、批次或时间设计。

多分类要说明one-vs-rest/one-vs-one及micro、macro或加权平均；macro使小类同等重要，micro汇总所有决策，数值回答不同问题。不同研究若阳性定义、难度与采样不同，AUC不宜直接排名。

## 排错与自测

| 症状 | 检查 |
|---|---|
| AUC远低于0.5 | 分数方向、阳性编码、样本对齐 |
| accuracy极高但recall=0 | 多数类比例与始终判阴基线 |
| AP与另一个工具不同 | 并列、插值/积分约定、权重 |
| 测试结果异常完美 | 泄漏、同供体重复、测试集调参 |

**问：AUC=0.75等于75%样本分类正确吗？** 答：不是，它是排序性质。**问：AP高就可直接用0.5为概率阈值吗？** 答：不可以，需校准和决策设计。**问：四个测试样本足够稳定估计性能吗？** 答：不足，本例只验证算术。


## 接着读什么

[组学验证的层级与证据](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%BB%84%E5%AD%A6%E9%AA%8C%E8%AF%81%E7%9A%84%E5%B1%82%E7%BA%A7%E4%B8%8E%E8%AF%81%E6%8D%AE)。学习顺序与练习入口见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。


## 运行本页的完整小例

先下载项目并进入含README与examples的根目录，环境准备见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。用项目环境的Python运行以下命令（macOS/Linux路径示例；Windows用`.venv\Scripts\python.exe`）：

```bash
.venv/bin/python examples/methods_foundations_demo.py --case classification --output examples/output/methods-classification
```

输出summary.json保存本页计算与检查结果。目录非空时拒绝覆盖，再练习请更换`--output`目录；页面内短代码只打印结果。标准库case也可用系统Python直接运行。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
