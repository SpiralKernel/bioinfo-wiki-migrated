# 中文生信 Wiki

面向入门学生，从群体分析、GWAS和GS出发，连接术语解释、算法、组学工作流和实践。

**阅读入口：[GitHub Wiki](https://github.com/SpiralKernel/bioinfo-wiki/wiki)** · **[关键词索引](https://github.com/SpiralKernel/bioinfo-wiki/wiki/关键词索引)**

## 从哪里开始

- [方法与算法系统学习指南](wiki/方法与算法系统学习指南.md)：70篇方法按七条路线学习，含推导、算例、可执行代码、排错和自测。
- [读懂方法公式的最小数学基础](wiki/读懂方法公式的最小数学基础.md)：先用一张表理解矩阵、对数和概率符号。

- [变异注释与功能解释学习路线](wiki/变异注释与功能解释学习路线.md)：基因模型、转录本、后果提取和候选证据。
- [生信缺口清单与维护优先级](wiki/生信缺口清单与维护优先级.md)：尚需补齐的真实数据和专业软件任务。
- [全类型变异检测总览与证据地图](wiki/全类型变异检测总览与证据地图.md)：小变异、SV、CNV/PAV、STR/VNTR、TE和组装分支。
- [人工长读长贯通SV调用与评估](wiki/人工长读长贯通SV调用与评估.md)：真实执行minimap2、Sniffles/SNF/force、Truvari两事件核对。
- [官方TRGT小数据贯通重复分型](wiki/官方TRGT小数据贯通重复分型.md)：固定PacBio提供例，实际重复分型。
- [变异检测分型与提取学习路线](wiki/变异检测分型与提取学习路线.md)：SNV/Indel、SV、CNV/PAV、STR、TE及生殖系/体细胞场景。
- [人工VCF贯通变异分类与提取](wiki/人工VCF贯通变异分类与提取.md)：16条人工记录、真实bcftools操作与可核对输出。

- [R与Python生信编程入门路线](wiki/R与Python生信编程入门路线.md)：下载项目、运行脚本、语法、表格、绘图、环境与排错。
- [bulkRNA-seq分析与实践路线](wiki/bulkRNA-seq分析与实践路线.md)：原始计数、配对设计、DESeq2与诊断。
- [公开airway数据贯通DESeq2](wiki/公开airway数据贯通DESeq2.md)：8个真实公开计数样本、原生R实际运行。
- [组学分析方法学习路线](wiki/组学分析方法学习路线.md)：RNA、单细胞、表观组、蛋白、代谢、微生物与多组学方法。
- [群体遗传学分析学习路线](wiki/群体遗传学分析学习路线.md)：多样性、FST、祖源、近交、LD、选择和历史。

- [测序技术与PCR学习路线](wiki/测序技术与PCR学习路线.md)：PCR、Sanger、二代、长读长、建库与质控。
- [人工区间贯通测序指标](wiki/人工区间贯通测序指标.md)：可运行的深度、双端重叠、重复、Phred和N50算例。

- [统计学与组学验证学习路线](wiki/统计学与组学验证学习路线.md)：概率、估计、基础检验及各组学验证。
- [生信基础算法学习路线](wiki/生信基础算法学习路线.md)：数据结构、查找、图、动态规划、聚类与优化。
- [组学验证的层级与证据](wiki/组学验证的层级与证据.md)：计算、统计、技术、独立样本和功能证据。

- [中文生信学习地图](wiki/中文生信学习地图.md)：前置知识与阅读入口。
- [群体分析与GWAS学习路线](wiki/群体分析与GWAS学习路线.md)：质控、PCA、模型与关联解释。
- [GWAS软件选择与比较](wiki/GWAS软件选择与比较.md)：PLINK、GEMMA、GCTA、REGENIE、SAIGE、BOLT及育种常用工具。
- [公开小鼠数据贯通多软件GWAS](wiki/公开小鼠数据贯通多软件GWAS.md)：真实公开数据、四个软件、七种配置及结果诊断图。
- [GS算法软件与实践路线](wiki/GS算法软件与实践路线.md)：算法、六种软件、验证设计与候选预测。
- [合成家系数据贯通GS](wiki/合成家系数据贯通GS.md)：300个人工个体、GBLUP与岭回归等价、嵌套验证和60个候选预测。
- [组学分析总览与选型](wiki/组学分析总览与选型.md)：按研究问题选择测量与流程。
- [合成数据贯通PCA与GWAS](wiki/合成数据贯通PCA与GWAS.md)：本项目实际运行的练习。
- [BLAST原理与种子延伸](wiki/BLAST原理与种子延伸.md)：人工字符串教学实现。

各篇标明更新时间与验证范围。组学实践包含官方数据、检查点和详细教程入口；本机没有运行的流程已明确标注。

## 运行教学例子

先[下载项目ZIP](https://github.com/SpiralKernel/bioinfo-wiki/archive/refs/heads/main.zip)并解压，在含README.md与examples的根目录运行。Windows、macOS/Linux和RStudio的具体进入目录步骤见[路径与项目组织](wiki/R与Python路径文件与项目组织.md)。

第一张表无需第三方包。macOS/Linux用`python3 examples/programming_basics.py`，Windows PowerShell用`python examples/programming_basics.py`；R控制台在项目根目录执行`source("examples/programming_basics.R")`。两张乱序人工表均得到8行、1个缺失、A有效n=3/均值2、B有效n=4/均值5.5。重复执行覆盖`examples/output/programming/`对应语言的同名文件。

需要第三方依赖时统一用项目环境；以下为macOS/Linux，Windows完整对应命令见[R包Python包与可复现环境](wiki/R包Python包与可复现环境.md)：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-teaching.txt
.venv/bin/python examples/programming_basics.py --plot
.venv/bin/python examples/population_genetics_demo.py
```

群体案例内置8条完整单倍型×12位置，输出到`examples/output/population-genetics/`；导出的TSV为输出副本。教学指标不是专业WC/Hudson、ROH或历史推断实现。Python与群体案例已运行，base R完整表格/模型/绘图脚本通过WebR R4.6.0执行；2026-10-15又在原生R4.4.3复跑基础R脚本，并实际执行airway/DESeq2/apeglm案例；其他专业包仍按各页验证范围。结果、两位审阅者的反馈与修改见[本轮维护记录](maintenance/omics-population-programming-2026-10-14.md)。

本轮变异案例需要bcftools1.22与+fill-tags插件，基础Python代码只用标准库。先按[bcftools环境页](wiki/bcftools变异操作与表达式.md)建立环境；以下从项目根目录运行：

```bash
conda run -n variants python examples/variant_extraction_demo.py
python3 examples/pav_matrix_demo.py
# 配图另需已安装Matplotlib的项目Python环境：
.venv/bin/python examples/plot_variant_extraction.py
```

固定人工输入在`examples/data/variants/`，输出分别为`examples/output/variants/`与`examples/output/pav-matrix/`，重复执行覆盖同名输出。VCF案例实际运行了标准化、类型/区域/样本提取、GT掩蔽与频率重算，还用40条人工SAM比对调用一个已知杂合SNV；未执行真实read mapper或专业SV/CNV/STR/TE调用器。PAV脚本只整理已经人工指定的状态，保留未知为NA，不发现基因缺失。详情见[变异首轮记录](maintenance/variants-round1-2026-10-14.md)。

新增注释、真实mapper调用链和工作流案例：Java/SnpEff/bedtools环境准备见[注释案例](wiki/人工基因贯通区间注释与后果提取.md)，minimap2/samtools/bcftools环境见[FASTQ案例](wiki/从FASTQ到比对检查与变异调用.md)。确认工具在PATH并设置SNPEFF_JAR为核心jar的绝对路径后，在项目根目录运行：

```bash
python3 examples/annotation_demo.py
python3 examples/genome_alignment_demo.py
# 绘图另需项目Matplotlib环境：
.venv/bin/python examples/plot_annotation.py
```

输出分别在examples/output/annotation、alignment-chain、annotation-plots，重复执行覆盖同名文件。固定输入是1600bp人工参考、三条转录本、12条VCF及40条人工FASTQ；注释VCF与reads调用结果是两个不同任务。实际验证3条CDS/蛋白零错误、负链与多ALT关联，以及预设杂合SNV。Snakemake9.27.0的安装、首次运行与缺失表恢复见[工作流案例](wiki/最小注释工作流与Snakemake恢复.md)，验证器仅使用新的独立目录。

本轮VEP、版本转换、多倍体/细胞器专业分析未执行。结果、工具归档SHA和范围见[注释与工作流记录](maintenance/annotation-workflow-round1-2026-10-14.md)。

真实RNA计数案例使用固定airway包，已有Conda时在项目根目录运行：

```bash
conda env create -f examples/environments/rnaseq.yml
conda run -n wiki-rnaseq Rscript examples/test_rnaseq_inputs.R
conda run -n wiki-rnaseq Rscript examples/rnaseq_airway_demo.R
```

本机Linux原生R4.4.3、DESeq2 1.46.0、airway1.26.0与apeglm1.28.0实际安装运行。8个样本配对设计，保留22,369个基因；普通Wald FDR<0.05为4,030个，绝对log2FC超过1的专门检验为243个。完整结果在examples/output/rnaseq-airway，非空输出目录会拒绝运行；传新的目录可另作练习。重画已有结果用`conda run -n wiki-rnaseq Rscript examples/plot_rnaseq_airway.R`，只更新图。

输入来自已有公开计数，未重跑FASTQ、STAR、featureCounts、Salmon/tximport或富集。基础R输入检查无需第三方包；完整检查用RSCRIPT_BIN指定原生Rscript。环境、输出和边界见[本轮记录](maintenance/rnaseq-round1-2026-10-15.md)。

本轮专业调用微型案例在Linux x86_64隔离环境实际执行，在项目根目录运行：

```bash
conda env create -f examples/environments/variant-calling.yml
conda run -n wiki-variant-calling python examples/variant_calling_demo.py
conda run -n wiki-variant-calling python examples/variant_calling_demo.py --case small
# TRGT先按对应工具页取得官方二进制，并设置TRGT_BIN为绝对路径：
conda run -n wiki-variant-calling python examples/trgt_official_demo.py
```

SV/FreeBayes输入完全人工，TRGT输入为固定官方fixture并核对SHA。输出分别为examples/output/calling-sv、calling-small、trgt-official；非空目录拒绝覆盖，重复练习用--output换新目录。Sniffles2.8.1、Truvari5.5.0、FreeBayes1.3.10、TRGT5.1.0均已执行；Sniffles force格式/继承样本表头的问题保留raw与受限.compat副本，GQ=0和群体缺失如实报告，不外推真实检测性能。

图另需Matplotlib，可在环境内安装后执行`conda run -n wiki-variant-calling python examples/plot_variant_calling.py`。完整回归还需requirements-teaching.txt中的科学计算包，以及各集成测试对应工具；没有时明确跳过。范围、两位子代理问题及修订见[本轮记录](maintenance/variant-detection-round2-2026-10-15.md)。

Python 3.9+ 可运行序列比对教学例子；它不依赖第三方包，不是完整BLAST实现：

```bash
python3 examples/alignment_demo.py
```

群体练习另需[PLINK 2](https://www.cog-genomics.org/plink/2.0/)，使用固定种子生成120个人工样本和600个位点：

```bash
bash examples/run_population_demo.sh
# 若不在PATH：
PLINK2_BIN=/完整路径/plink2 bash examples/run_population_demo.sh
```

输出在 `examples/output/population/`，重复运行覆盖其中同名结果。详见[实测记录](maintenance/validation-2026-10-13.md)。二进制、练习输出和原稿私有备份不提交。

完整公开小鼠案例的工具依赖、运行方式与验证范围见[第二轮GWAS记录](maintenance/gwas-round2-2026-10-13.md)。代码为 `examples/gwas_software_demo.py` 和 `examples/plot_gwas_comparison.py`；公开配图与汇总指标在 `assets/gwas-round2/`，个体数据仍不纳入仓库。

GS教学例子依赖NumPy，配图另需Matplotlib：

```bash
.venv/bin/python examples/genomic_prediction_demo.py
.venv/bin/python examples/plot_genomic_prediction.py
```

输出在 `examples/output/gs/`，重复运行覆盖同名文件；建议使用项目虚拟环境。实际结果和验证范围见[GS第一轮记录](maintenance/gs-round1-2026-10-13.md)。专业GS软件的模板尚未运行，教学实现不包含REML、MCMC或ssGBLUP。

统计与基础算法练习建议Python 3.10+，依赖NumPy、SciPy 1.15+（rng接口）和Matplotlib：

```bash
.venv/bin/python examples/statistics_demo.py
.venv/bin/python examples/basic_algorithms_demo.py
```

[统计案例](wiki/合成数据贯通统计推断.md)实际运行了均值、配对、分类表、区间、置换、多重校正与供体层级反例；[算法练习](wiki/基础算法可复现练习.md)包含六类小型实现。两者均为人工数据，原始输出分别在 `examples/output/statistics/` 和 `examples/output/basic-algorithms/`，重复运行覆盖同名文件。版本与边界见[本轮维护记录](maintenance/statistics-round1-2026-10-13.md)。分组学验证页提供设计指南，未执行湿实验。

测序指标算例只需Python 3.9+；配图另需Matplotlib：

```bash
python3 examples/sequencing_metrics_demo.py
python3 examples/plot_sequencing_metrics.py
```

输出在`examples/output/sequencing/`，重复运行覆盖同名文件。区间和分子身份是人工真值，没有运行比对或去重工具。PCR、建库、测序仪及外部质控软件页为原理和判读指南，尚未执行实验。范围与结果见[本轮测序记录](maintenance/sequencing-round1-2026-10-13.md)。

## 方法与算法的配套练习

原有66篇方法已逐篇深化，新增EM、HMM、BWT/FM索引和de Bruijn图；系统指南给出前置、共同基础与七条专题路线。教学结构参考GTN、Bioinformatics Workbook、SeqAn与Truvari Wiki，正文和人工例独立撰写。三套runner共70个方法算例，部分基础例需要NumPy/SciPy，其余只用标准库：

```bash
# 按前文建立.venv并安装requirements-teaching.txt，在项目根目录运行：
.venv/bin/python examples/methods_foundations_demo.py --list
.venv/bin/python examples/methods_foundations_demo.py --case pca --output examples/output/methods-pca
python3 examples/methods_genetics_demo.py --list
python3 examples/methods_genetics_demo.py --case gblup --output examples/output/methods-gblup
python3 examples/methods_omics_demo.py --list
python3 examples/methods_omics_demo.py --case pseudobulk --output examples/output/methods-pseudobulk
# 两张自绘图：
.venv/bin/python examples/plot_methods_foundations.py
```

Windows PowerShell使用`.venv\Scripts\python.exe`替换`.venv/bin/python`，系统Python命令也可使用已安装的`python`。运行目录需含README与examples。保存结果时拒绝覆盖非空目录，重复练习用新的--output；genetics不传--output只打印结果，omics默认按case写入其输出目录。不要把教学算术视为执行了REML、SKAT、真实组装或测序实验，专业范围逐页声明。结果、图、环境与双代理审阅见[本轮维护记录](maintenance/methods-depth-2026-10-26.md)。

## 在Obsidian编写并发布

原稿在 Obsidian 的 `知识库/中文生信Wiki`；本仓库 `wiki/` 是生成的发布版。写作属性和维护规则见 [CONTRIBUTING.md](CONTRIBUTING.md)。

```bash
python3 scripts/export_wiki.py --source '/你的Obsidian库/知识库/中文生信Wiki' --refresh-index
# 完整集成检查另需bcftools及+fill-tags、Java/SnpEff、bedtools、minimap2和samtools；缺少工具时相应项目明确跳过：
conda run -n variants .venv/bin/python -m unittest discover -s tests -v
bash scripts/publish_wiki.sh
```

导出只读取六类目录：学习路线、术语与概念、方法与算法、工具、数据库、分析工作流。草稿和待编写页面默认不导出；正文属性会转换为公开的状态、更新时间及验证说明，首页、侧边栏和关键词索引自动生成。默认不修改原稿；`--refresh-index` 只刷新原稿中的自动索引。模板、写作入口、私有备份和其他Obsidian笔记不导出。

内部链接支持标题、相对笔记路径、Wiki根路径和无歧义别名；代码块与行内代码中的双括号保持原样。断链、重名、缺少必要属性、未闭合代码块会停止导出。Obsidian嵌入和标题锚点暂不支持；附件不自动复制，添加图片前需自行安排公开资源及检查标准Markdown路径。

发布脚本使用Git与已登录的GitHub CLI，在旁边的独立 `bioinfo-wiki.wiki` 仓库拉取、校验、同步并推送。仅同步清单中的文件；删除或改名只清理此前受管理的旧页，保留其他页面。若受管理页面曾在网页上修改，发布会停止，需先把差异合回Obsidian原稿。首次从旧版升级也会核对旧文件校验值。初次启用GitHub Wiki时需要先在网页创建第一页，本项目已完成。

## 项目结构

```text
examples/                  人工数据与算法教学代码
scripts/export_wiki.py     校验链接、导出文章、生成导航和索引
scripts/sync_wiki_files.py  校验清单并保护远程人工修改
scripts/publish_wiki.sh     发布到Wiki独立仓库
tests/                     教学算法、字段语义、导出与安全同步检查
maintenance/               验证记录（私有备份目录除外）
wiki/                      可检查和提交的发布文件
```
