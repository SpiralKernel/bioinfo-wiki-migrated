# 蛋白质鉴定FDR与差异丰度

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26官方/原始来源核对；标准库教学演示实际执行；专业验证范围见正文。

## 两道统计关口

“这个谱图对应哪条肽段”和“这个蛋白在条件间是否改变”使用不同数据与零假设。鉴定层FDR不能替代差异层FDR。学完应能画出谱图→PSM→肽段→蛋白组→样本丰度的层级，并解释共享肽段及未检出。前置读[P值效应量与多重检验](https://github.com/SpiralKernel/bioinfo-wiki/wiki/P%E5%80%BC%E6%95%88%E5%BA%94%E9%87%8F%E4%B8%8E%E5%A4%9A%E9%87%8D%E6%A3%80%E9%AA%8C)、[缺失值与零值](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%BC%BA%E5%A4%B1%E5%80%BC%E4%B8%8E%E9%9B%B6%E5%80%BC)及[生物学重复与技术重复](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%94%9F%E7%89%A9%E5%AD%A6%E9%87%8D%E5%A4%8D%E4%B8%8E%E6%8A%80%E6%9C%AF%E9%87%8D%E5%A4%8D)。

## 鉴定错误如何估计

PSM是谱图与候选肽段的匹配。搜索目标序列和构造的诱饵decoy序列时，诱饵命中可以为错误匹配提供估计线索；前提包括目标/诱饵搜索空间、竞争规则与错误得分行为。[Elias与Gygi原始研究](https://pubmed.ncbi.nlm.nih.gov/17327847/)讨论了target-decoy的假设和设计。它不表示某一个被接受肽段“恰有1%的概率是错的”，也不自动把PSM层阈值转成蛋白层FDR。

在一个假定目标/诱饵规模和错误竞争可比的教学约定下，某阈值上100个目标命中、2个诱饵命中，可用D/T=0.02直观估计错误比例量级。真实流程可能采用不同竞争、加一修正或层级控制，不能未经核对把2/100写成所有软件统一算法，也不能把这个toy当作已控制FDR的证明。

蛋白推断又有共享肽段：同一序列可能对应多个蛋白或同工型。蛋白组是按证据可区分性组织的集合，不代表每个候选同工型均被独立测到。增加数据库、放宽修饰或搜索范围会改变候选空间与错误控制，不是“搜索越全越好”而无需代价。

## 从肽段强度到条件效应

设$y_{psj}$是蛋白$p$的肽段/特征$s$在样本$j$的log2强度，可用示意模型 $y_{psj}=\mu_p+\eta_{ps}+\beta_pT_j+u_{pj}+\epsilon_{psj}$。$\eta$代表肽段响应差异，$T$是条件，$u$为生物样本层变化；模型选择取决于DDA/DIA、标记方式和实验层级，不能把同一供体多个谱图当独立供体。[MSstats文档](https://bioconductor.org/packages/MSstats/)提供对应工作流入口，本页不提供统一搜索/定量参数。

若两条相应肽段的B/A强度比都为2，log2比值均1，简单中位数汇总得到蛋白log2FC=1。若其中一条是共享肽段或受修饰影响，就不该未经规则直接混入同一蛋白丰度。若某组未检出，也不等于真实强度0：检测限、随机采样、离子化和缺失机制需另外判断。

## 差异层的多重校正另做

四个蛋白原始P按序为0.001、0.01、0.04、0.2，BH先算$mP_{(i)}/i$再从尾向前取累计最小并限于1，得到0.004、0.02、约0.05333、0.2。鉴定的2%错误估计不能代替这一步；反过来差异padj小也不能确认搜索身份。

## 输出、参数与排错

同时交付鉴定层级、数据库/搜索设定、共享肽段规则、定量单位、缺失/归一化、样本设计、效应与区间/P/FDR。每一层过滤都保留数量与理由。

| 症状 | 原因候选 | 检查 |
|---|---|---|
| 蛋白数量随搜索扩张大增 | 候选与错误空间变大 | 层级FDR、数据库和诱饵规则 |
| 补零后出现极大倍数 | 未检出被当真零 | 检测限、原强度与缺失机制 |
| 同一蛋白各肽段方向相反 | 共享、修饰或特征异常 | 肽段归属、特异性和原谱图 |
| P极小但只有一供体 | 技术记录被当重复 | 生物学单位与模型自由度 |

## 自测与答案

1. PSM通过1%FDR是否证明蛋白也1%？答：层级不同，需分别控制。
2. 蛋白丰度增加是否等于活性增强？答：修饰、定位与反应状态另需证据。
3. 共享肽段能确认所有同工型吗？答：不能超过证据可区分范围。

## 在电脑上核对这一小步

先[下载项目ZIP](https://github.com/SpiralKernel/bioinfo-wiki/archive/refs/heads/main.zip)并解压，进入同时含README.md和examples的项目根目录。不熟悉终端和当前目录时读[R与Python路径文件与项目组织](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E4%B8%8EPython%E8%B7%AF%E5%BE%84%E6%96%87%E4%BB%B6%E4%B8%8E%E9%A1%B9%E7%9B%AE%E7%BB%84%E7%BB%87)。以下在**系统终端**执行，需Python 3.10+，只用标准库，无需pip安装；Windows可将python3替换为已安装的python或py。

```bash
python3 examples/methods_omics_demo.py --case proteomics
```

输入是脚本内置人工数字，不需下载专业数据。成功会打印JSON，并生成`examples/output/methods-omics/proteomics/summary.json`；`cases.proteomics`下的字段对应本页演算。非空输出目录会拒绝覆盖，第二次练习可在命令后加`--output examples/output/methods-omics/proteomics-practice2`。不要为了重跑删除不明来源的数据。

**检查点：**illustrative_decoy_over_target=0.02；BH=[0.004,0.02,约0.0533333,0.2]；protein_median_log2_ratio=1。

阅读源码时先找到同名函数，逐一对照输入、运算与输出；要修改数字，先复制整个解压项目并运行副本，再手算新的预期。脚本含若干固定结果断言，修改教学条件时也需重新核对相应断言；断言失败不自动说明专业方法失效。本轮已实际运行全部22个标准库模式，Python 3.14.8；本页范围为诱饵比例示意、肽段倍数与BH算术；未运行质谱搜索、蛋白推断或MSstats。演示没有为真实样本产生专业P值或性能结论。


## 来源与下一步

本轮查阅日期为2026-10-26；正文原始来源支持对应技术事实，手算与代码是本项目人工教学设计。方法页解释模型，安装及大数据运行交给相应工具/工作流页。教学组织参考[GTN](https://training.galaxyproject.org/training-material/)的步骤与检查点、[Bioinformatics Workbook](https://bioinformaticsworkbook.org/)的问题驱动练习；未复制其教程正文。

保留的相关入口：[蛋白质组与代谢组验证设计](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%9B%8B%E7%99%BD%E8%B4%A8%E7%BB%84%E4%B8%8E%E4%BB%A3%E8%B0%A2%E7%BB%84%E9%AA%8C%E8%AF%81%E8%AE%BE%E8%AE%A1)、[target-decoy原始论文](https://www.nature.com/articles/nmeth1019)、[MSstats官方入口](https://bioconductor.org/packages/MSstats/)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
