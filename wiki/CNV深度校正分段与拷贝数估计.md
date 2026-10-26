# CNV深度校正分段与拷贝数估计

> 更新：2026-10-26；状态：可用；验证范围：2026-10-26官方/原始来源核对；标准库教学演示实际执行；专业验证范围见正文。

## 先说明从深度能知道多少

CNV关注拷贝剂量，read depth提供相对信号，但不是整数CN的直接读数。学完应能从bin深度算log2比值、解释分段目标及纯度影响，并分清区域剂量与精确断点。前置读[测序深度覆盖度与质量分数](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%B5%8B%E5%BA%8F%E6%B7%B1%E5%BA%A6%E8%A6%86%E7%9B%96%E5%BA%A6%E4%B8%8E%E8%B4%A8%E9%87%8F%E5%88%86%E6%95%B0)、[结构变异CNV与PAV](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E7%BB%93%E6%9E%84%E5%8F%98%E5%BC%82CNV%E4%B8%8EPAV)、[拷贝数基线非整倍体与LOH](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%8B%B7%E8%B4%9D%E6%95%B0%E5%9F%BA%E7%BA%BF%E9%9D%9E%E6%95%B4%E5%80%8D%E4%BD%93%E4%B8%8ELOH)。

## 从可比窗口生成信号

将可分析区域分成bin或捕获靶区，得到样本覆盖$d_i$；匹配正常参考或模型基线$b_i$后，产生相对比值。简单示意为 $r_i=\log_2[(d_i/s)/b_i]$，$s$为样本尺度。真实工作还要处理GC、可比对性、异常区域、重复和测量权重；对极低/零覆盖不能机械取log再当高置信缺失。[CNVkit流程](https://cnvkit.readthedocs.io/en/stable/pipeline.html)将覆盖、参考、校正与分段分开。

WGS、捕获与扩增子可分析区域和偏差不同，不能复用同一baseline。性染色体、区域倍性、非整倍体和正常对照质量都会改变中性水平。匹配参考是一个估计任务，不是随便选一个样本深度作分母。

## 手算纯净样本与混合样本

基线深度100，六个bin深度100、100、50、50、150、150，log2比值为0、0、−1、−1、约0.585、0.585。在纯二倍体、正确尺度和理想测量下，$CN=2\times2^r$，对应2、2、1、1、3、3。[CNVkit拷贝数说明](https://cnvkit.readthedocs.io/en/stable/calling.html)明确这种换算依赖纯度与倍性。

若一半细胞CN=1、一半正常CN=2，混合平均为1.5份，相对二倍体比值0.75，log2≈−0.415，而非−1。更一般的简化关系为 $R=[\rho C+(1-\rho)2]/2$，$\rho$为异常细胞比例、$C$为其CN；它未包含整体倍性、亚克隆和归一化复杂性，不能凭一个log2值同时精确反推所有量。

## 分段是有代价的合并

噪声bin需要组织成区域，可用概念目标

$$
\min_{segments}\sum_s\sum_{i\in s}w_i(r_i-\theta_s)^2+\lambda M.
$$

$\theta_s$是第$s$段的信号水平，$w_i$是可靠性权重，$M$为段数，$\lambda$惩罚过多分段。它是分段思想的示意，不声称CNVkit的所有算法正好优化此式。先形成候选边界、评估相邻信号差，再按误差与复杂度合并/拆分，最后赋予剂量解释。$\lambda$过大可漏短事件，过小可把噪声切成很多段；bin更小也会降低每窗计数与稳定性。

## 输入输出与检查

输入应有相同参考的排序索引BAM、样本/对照、目标/可访问区域及倍性设计。输出保留bin信号、校正值、分段边界、CN/质量、基线和可调用区间。深度不能独立确认平衡倒位、精确连接或拷贝中性LOH；可结合等位平衡、分裂reads、组装和其他证据。

| 症状 | 可能原因 | 先查 |
|---|---|---|
| 高GC区域周期性波动 | 覆盖偏差 | GC/可比对性与残差 |
| 分段很多且很短 | 噪声或过度切分 | bin计数、权重与参数敏感性 |
| 整体水平漂移 | 基线/纯度/倍性错误 | 正常参考和染色体分布 |
| CN=1被标整基因不存在 | 剂量与PAV混淆 | 单位、另一拷贝和基因覆盖范围 |

## 自测与答案

1. log2=-0.415是否一定CN=1.5？答：只是相对混合信号，整数状态与纯度需模型。
2. 深度零能否直接CN=0？答：还可能不可比对或没测到，应保留未知。
3. 分段边界就是断点吗？答：受bin和方法分辨率限制，需连接证据定位。

## 在电脑上核对这一小步

先[下载项目ZIP](https://github.com/SpiralKernel/bioinfo-wiki/archive/refs/heads/main.zip)并解压，进入同时含README.md和examples的项目根目录。不熟悉终端和当前目录时读[R与Python路径文件与项目组织](https://github.com/SpiralKernel/bioinfo-wiki/wiki/R%E4%B8%8EPython%E8%B7%AF%E5%BE%84%E6%96%87%E4%BB%B6%E4%B8%8E%E9%A1%B9%E7%9B%AE%E7%BB%84%E7%BB%87)。以下在**系统终端**执行，需Python 3.10+，只用标准库，无需pip安装；Windows可将python3替换为已安装的python或py。

```bash
python3 examples/methods_omics_demo.py --case cnv
```

输入是脚本内置人工数字，不需下载专业数据。成功会打印JSON，并生成`examples/output/methods-omics/cnv/summary.json`；`cases.cnv`下的字段对应本页演算。非空输出目录会拒绝覆盖，第二次练习可在命令后加`--output examples/output/methods-omics/cnv-practice2`。不要为了重跑删除不明来源的数据。

**检查点：**log2_ratios为[0,0,-1,-1,约0.584963,约0.584963]，ideal_CN=[2,2,1,1,3,3]；混合CN1信号约-0.415037。

阅读源码时先找到同名函数，逐一对照输入、运算与输出；要修改数字，先复制整个解压项目并运行副本，再手算新的预期。脚本含若干固定结果断言，修改教学条件时也需重新核对相应断言；断言失败不自动说明专业方法失效。本轮已实际运行全部22个标准库模式，Python 3.14.8；本页范围为已知基线与纯度换算；未执行GC校正、分段拟合或真实绝对CN推断。演示没有为真实样本产生专业P值或性能结论。


## 来源与下一步

本轮查阅日期为2026-10-26；正文原始来源支持对应技术事实，手算与代码是本项目人工教学设计。方法页解释模型，安装及大数据运行交给相应工具/工作流页。教学组织参考[GTN](https://training.galaxyproject.org/training-material/)的步骤与检查点、[Bioinformatics Workbook](https://bioinformaticsworkbook.org/)的问题驱动练习；未复制其教程正文。

保留的相关入口：[CNV检测与结果提取](https://github.com/SpiralKernel/bioinfo-wiki/wiki/CNV%E6%A3%80%E6%B5%8B%E4%B8%8E%E7%BB%93%E6%9E%9C%E6%8F%90%E5%8F%96)、[CNVkit流程](https://cnvkit.readthedocs.io/en/stable/pipeline.html)、[CNVkit拷贝数解释](https://cnvkit.readthedocs.io/en/stable/calling.html)。

学习顺序、前置与练习清单见[方法与算法系统学习指南](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E6%96%B9%E6%B3%95%E4%B8%8E%E7%AE%97%E6%B3%95%E7%B3%BB%E7%BB%9F%E5%AD%A6%E4%B9%A0%E6%8C%87%E5%8D%97)。

公式里的矩阵、下标或对数不熟时，可回查[读懂方法公式的最小数学基础](https://github.com/SpiralKernel/bioinfo-wiki/wiki/%E8%AF%BB%E6%87%82%E6%96%B9%E6%B3%95%E5%85%AC%E5%BC%8F%E7%9A%84%E6%9C%80%E5%B0%8F%E6%95%B0%E5%AD%A6%E5%9F%BA%E7%A1%80)。
