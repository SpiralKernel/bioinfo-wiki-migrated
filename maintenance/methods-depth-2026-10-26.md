# 方法与算法系统深化（2026-10-26）

本轮回应入门读者反馈：原有方法条目多为定义与提示，缺少能够跟随的推理、数值过程和结果检查。逐篇深化原有66篇方法，新增EM、HMM、BWT/FM索引和de Bruijn图4篇方法；再新增系统指南和公式阅读桥接页。全库正文318→324，方法66→70，管理327个Wiki Markdown文件。新增标题少于本轮修订，重点是加深现有解释。

## 借鉴哪些知识库，具体改变什么

| 参考 | 本轮借鉴的组织方式 | 在本项目的落实 |
|---|---|---|
| [Galaxy Training Network教程写法](https://training.galaxyproject.org/training-material/topics/contributing/tutorials/create-new-tutorial/tutorial.html) | 明确目标、可管理的小输入、操作与问题/解答 | 每页有学习目标、数值例、精确case与检查点、自测答案 |
| [Bioinformatics Workbook](https://bioinformaticsworkbook.org/about.html) | 从真实问题组织解释与练习，再连接背景 | 先说明解决的问题与前置，后接已有工作流及真实案例 |
| [SeqAn成对比对教程](https://seqan.readthedocs.io/en/main/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html) | 说明前置、评分、输出与练习变体 | 全局/局部比对补矩阵、回溯、评分复核和改变条件的题目 |
| [Truvari Wiki bench](https://github.com/ACEnglish/truvari/wiki/bench) | 解释匹配参数的作用和评估对象 | 先定义事件/区域/阈值与TP两侧，再做一对一匹配反例 |

只借鉴教学结构和可复核习惯；正文、数值例与图独立编写，没有复制外部教程或图表。数学原理与技术事实分别在条目关键位置引用原论文或官方资料，查阅日期为2026-10-26；被验证码/权限阻挡的原论文没有冒称已读全文，保留为进一步阅读并用可读取的作者文档核对操作和公式。

## 知识结构与深度

[系统指南](../wiki/方法与算法系统学习指南.md)给十二个共同基础检查点，以及计算/序列、统计推断、几何/隐变量/预测、群体/GWAS/GS、变异/组装、组学、测量技术七条路线。每条路线同时写下一步工作流与完成标志。首页与侧边栏提供直接入口，避免只剩一长串互不相连的标题。

[数学桥接页](../wiki/读懂方法公式的最小数学基础.md)用同一张三样本表解释求和、下标、向量矩阵、转置、乘法、秩、对数、条件概率及常见模型记号；涉及更深推导时链接到MIT/Harvard课程与对应方法。该页是阅读帮助，不是完整数学教材。

原有66篇方法均增加具体问题、模型或算法步骤、符号/单位解释、数值例、可复制计算、参数影响、假设、排错、自测与答案。复杂专业方法按层次解释，人工例只核对其中可计算的一步；软件安装和生产调用继续由工具/工作流承接。逐页改动及case映射见[method-pages.json](../assets/methods-depth-20261026/method-pages.json)；正文长度变化只作审计，不作为质量评分。

写作规范升至v0.3，方法模板与CONTRIBUTING同步补齐深度要求、代码/正文一致性、数值重核与验证范围。原稿仍在Obsidian，项目wiki为生成发布版，没有建立另一份独立维护正文。

## 实际运行了哪些例子

| runner/内容 | 数量 | 实际验证的内容 |
|---|---:|---|
| [methods_foundations_demo.py](../examples/methods_foundations_demo.py) | 26个case、27个页面Python块 | 查找/图/比对、距离/PCA/聚类/优化、基础推断/GLM算术、EM/HMM/FM/已知图路径 |
| [methods_genetics_demo.py](../examples/methods_genetics_demo.py) | 22个case | FST/LD/ROH、祖源概率、BLUP/REML特例、历史尺度、SV匹配/状态、泄漏/正则化 |
| [methods_omics_demo.py](../examples/methods_omics_demo.py) | 22个case | 扩增/信号模型、计数/交互、供体/图/整合、富集/组成、CNV/区间/稀有变异分数 |
| 数学桥接页 | 1个短例 | 三样本设计矩阵维度、线性求解、方差与对数 |

70个runner算例均在最终脚本版本的新输出目录运行，数学短例另行执行。三个`--list`均与case集合相符；再次指定非空结果目录均拒绝运行，summary SHA保持不变。[CLI检查](../assets/methods-depth-20261026/cli-checks.json)。这些不计作70个新的独立专业软件测试。

基础例实际环境为Python3.11.17、NumPy2.4.6、SciPy1.17.1、Matplotlib3.11.2；标准库两套例使用Python3.14.8。复用既有隔离Linux环境，无新增全局安装；未逐平台重新安装依赖。[环境记录](../assets/methods-depth-20261026/environment.json)。读者依赖仍按requirements-teaching.txt准备，边界与运行目录逐页说明。

公开结果分别为[基础例](../assets/methods-depth-20261026/foundations-results.json)、[遗传/变异例](../assets/methods-depth-20261026/genetics-results.json)、[组学/测量例](../assets/methods-depth-20261026/omics-results.json)与[数学短例](../assets/methods-depth-20261026/math-primer-output.txt)。输入均为人工数字或字符串，无真实个体资料。脚本和结果SHA见[provenance.json](../assets/methods-depth-20261026/provenance.json)。

### 可复算的检查点

- 全局比对ACGT/AGT在匹配+2、错配−1、缺口−2下为4分；逐列重算并恢复输入。
- 秩一PCA解释率[1,0]、重建误差约0；符号翻转不改变轴。
- Welch与配对例分别t≈1.632993和4.898979，由采样设计决定用哪一个。
- 完整枚举六数分组20种，绝对均值差双侧置换P=0.1；不是Monte Carlo P=0。
- 混杂反例总体斜率约−1.303226，调整组别后为1；不由此证明因果。
- EM固定方差例责任度每行和1、更新后对数似然提高；没有完整高斯混合软件验证。
- HMM最佳路径[0,0,1]联合概率0.05184，8条路径总概率0.10452；最优路径与求和分开。
- FM检索ana在banana的0基命中1、3，与直接字符串检索核对。
- GBLUP两种路径固定同一尺度后预测±20/27；REML例仅为截距高斯模型自由度特例。
- SV一候选对两真值的一对一反例TP=1、FN=1，保留群体缺失与低GQ=0参考输出的不同。

[基础代码块SHA清单](../assets/methods-depth-20261026/foundations-cases.json)保证页面短例与runner相同；[逐块实际结果](../assets/methods-depth-20261026/inline-code-results.json)保存27块输出。源代码与正文双处维护，更新时须同步并复跑受影响case，约定已写进CONTRIBUTING。

## 两位子代理交叉纠错

两位代理各负责22篇深化及独立标准库算例，随后交叉审阅：读者代理复核主代理22篇、四篇新增方法、数学桥接与总指南；技术代理逐篇复核组学/测序22篇、全部相关算例和四篇新增方法。主代理核查群体/GS/变异内容与输出，执行整套与逐块验证，再整合修订。

| 提出的问题 | 核对与修正 |
|---|---|
| 读者：临时脚本放哪里，为什么from examples导入失败？ | 明确保存项目根目录，使用相同环境解释器；总指南区分首次建环境与已有环境，补全PowerShell及两个标准库实际case |
| 读者：log2是log₂还是ln(2)，PCA符号和数组轴不懂？ | GLM系数明确ln(2)；PCA解释n、Σ、转置、@与axis=0，局部比对补状态下标和评分符号 |
| 读者：如何确认检验输出，三套summary为什么不同？ | 补Welch/配对/Fisher的实际P检查点；指南列各runner的JSON路径、默认写出与只打印规则 |
| 读者：HMM代码复制整条路径也能直接声称最优复杂度吗？ | 区分优化回溯指针递推与教学代码的额外路径复制成本；BWT哨兵另加代码标记避免误渲染 |
| 主代理：GBLUP检查点中的s是实际JSON键吗？ | 改为scale（公式中的s），逐页核对case与键；给出原始剂量和中心化值以免把半整数误当硬GT |
| 技术：模块度的总边权是否错了一倍？ | 直接定义k_i=Σ_j A_ij、2m=Σ_ij A_ij，无向边在邻接矩阵计两次，消除含糊定义 |
| 技术：clr值为负时能否计算所列Bray–Curtis？ | 限定非负丰度且分母非零；明确clr不直接套该式，区分Aitchison/欧氏几何 |
| 技术：toy是否冒充仪器或专业统计结果？ | 逐页保留模型假设、专业未测与历史实测来源；HMM/EM/FM/图拼接均注明实现边界 |

上述修订均已落实，Python片段未因文字澄清而变化。审阅关注解释、公式、数值与操作对应，不把两个代理的判断替代真实样本验证。可复核的审阅摘要与来源清单见[review-notes.json](../assets/methods-depth-20261026/review-notes.json)。

## 图与维护检查

[绘图脚本](../examples/plot_methods_foundations.py)用同一人工输入生成比对矩阵/回溯与PCA中心化/投影两张图。已视觉检查坐标、箭头、数值和文字，无外部图片转载。[图与材料说明](../assets/methods-depth-20261026/README.md)。HMM页另有状态与观测关系图，区分模型结构和某次解码。

既有54项回归测试全部通过、无跳过，涵盖先前的教学算法、统计性质、专业微型案例、导出与安全同步；[记录](../assets/methods-depth-20261026/regression-tests.txt)。本轮教学断言、文档代码核对和CLI保护另行报告，不与54项混加。

导出校验正文属性、别名、内部链接与代码块，并检查327个受管理页面的同步碰撞。Obsidian写回前逐页核对初始SHA、备份旧稿；如有并发编辑则停止。Wiki发布仍保护远程人工修改。下载的专业软件源码、二进制、环境与私有备份不进入公开仓库；公开只保留教学脚本与汇总。

## 未验证范围与后续

本轮没有因为扩写而重新执行全部GWAS/GS或新增真实组学流程。没有执行FarmCPU/BLINK、专业REML/MCMC、ADMIXTURE/NGSadmix、PSMC/MSMC/SMC++、SKAT、真实组装、真实仪器/湿实验或各组学完整上游；对应页明确限定教学算术与概念范围。原有已实跑的公开GWAS、人工GS、DESeq2、Sniffles/TRGT等案例保留各自验证证据。

下一轮按学生实际练习反馈补前置与难点，并选真实数据闭环，尤其系统发育、实际组装、真实CNV/TE和组学原始数据。不能把本轮目录覆盖、数值断言或图形效果当成真实研究性能证明。

## 发布时核对仓库迁移

当前主仓库和Wiki的origin均已配置为`SpiralKernel/bioinfo-wiki-migrated`。发布时通过GitHub API确认旧名与新名指向同一公开仓库且Wiki启用；旧Wiki地址返回新名页面，旧raw图片地址仍返回200。发布脚本原先固定旧名并会拒绝当前Wiki远程，已改为从主仓库当前HTTPS origin推导Wiki地址，保留目标远程一致性与人工修改保护。历史正文中的旧链接继续通过已验证的GitHub重定向访问。
