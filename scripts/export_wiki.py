#!/usr/bin/env python3
"""Export authored notes, validate links and generate navigation with stdlib only."""
import argparse
import ast
from datetime import date
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile
from urllib.parse import quote

WIKI_URL = 'https://github.com/SpiralKernel/bioinfo-wiki/wiki'
FOLDERS = ('学习路线', '术语与概念', '方法与算法', '工具', '数据库', '分析工作流')
INDEX = '关键词索引'
MANIFEST = '.wiki-export.json'
LINK = re.compile(r'(!?)\[\[([^\]]+)\]\]')
RESERVED = {'Home', '_Sidebar', '_Footer'}


def page_url(name):
    return f'{WIKI_URL}/{quote(name, safe="")}'


def scalar(value):
    value = value.strip()
    if not value:
        return ''
    for parser in (json.loads, ast.literal_eval):
        try:
            return parser(value)
        except (ValueError, SyntaxError):
            pass
    return value


def read_note(text):
    """Read scalar properties and JSON/quoted inline or block alias lists."""
    lines = text.splitlines(keepends=True)
    if not lines or lines[0].strip() != '---':
        return {}, text
    end = next((i for i in range(1, len(lines)) if lines[i].strip() == '---'), None)
    if end is None:
        raise ValueError('YAML 属性缺少结束分隔符')
    properties = {}
    active = None
    for line in lines[1:end]:
        item = re.match(r'^\s+-\s+(.*)$', line)
        if item and active:
            if not isinstance(properties.get(active), list):
                raise ValueError(f'{active}: 不支持的列表属性格式')
            properties[active].append(scalar(item[1]))
            continue
        match = re.match(r'^([\w-]+):\s*(.*?)\s*$', line)
        if not match:
            continue
        key, value = match.groups()
        active = key if not value else None
        properties[key] = [] if active and key in {'aliases', 'tags'} else scalar(value)
    aliases = properties.get('aliases', [])
    if isinstance(aliases, str):
        aliases = [aliases] if aliases else []
    if not isinstance(aliases, list) or any(not isinstance(a, str) for a in aliases):
        raise ValueError('aliases 需要字符串或字符串列表，请使用模板中的属性格式')
    properties['aliases'] = aliases
    return properties, ''.join(lines[end + 1:]).lstrip('\n')


def strip_properties(text):
    return read_note(text)[1]


def transform_prose(text, transform):
    output = []
    fence = None
    for line in text.splitlines(keepends=True):
        marker = re.match(r'^\s{0,3}(`{3,}|~{3,})', line)
        if fence is not None:
            output.append(line)
            if re.match(rf'^\s{{0,3}}{re.escape(fence[0])}{{{len(fence)},}}\s*$', line):
                fence = None
            continue
        if marker:
            fence = marker[1]
            output.append(line)
            continue
        # Paired backtick spans can contain shorter runs of backticks.
        parts = re.split(r'(`+)(.*?)\1', line)
        rendered = []
        for i in range(0, len(parts), 3):
            rendered.append(transform(parts[i]))
            if i + 2 < len(parts):
                rendered.append(parts[i+1] + parts[i+2] + parts[i+1])
        output.append(''.join(rendered))
    if fence:
        raise ValueError('代码块未关闭')
    return ''.join(output)


def render_note(text, path, targets):
    def replace(match):
        if match[1]:
            raise ValueError(f'{path.name}: ![[...]] 嵌入请先转换为标准 Markdown')
        target, separator, label = match[2].replace('\\|', '|').partition('|')
        page, anchor, _ = target.partition('#')
        page = page.removesuffix('.md')
        if anchor:
            raise ValueError(f'{path.name}: 请使用稳定的页面链接替代标题锚点：{target}')
        if page not in targets:
            raise ValueError(f'{path.name}: 内部链接不在发布范围或别名有歧义：{target}')
        display = label if separator else page.rsplit('/', 1)[-1]
        return f'[{display}]({page_url(targets[page])})'
    return transform_prose(strip_properties(text), lambda line: LINK.sub(replace, line))


def table_cell(value):
    return str(value).replace('|', '\\|').replace('\n', ' ')


def build_index(records, obsidian=False):
    lines = ['# 关键词索引', '',
             '按标题、英文名或缩写查找。可以使用浏览器或 Obsidian 的查找功能；别名只是检索入口，不改变概念定义。', '',
             '| 条目 | 常见名称与缩写 | 类型 | 一句话说明 |', '|---|---|---|---|']
    for item in sorted(records, key=lambda r: r['title'].casefold()):
        title = item['title']
        link = f'[[{title}\\|{title}]]' if obsidian else f'[{title}]({page_url(title)})'
        aliases = '、'.join(item['aliases']) or '—'
        lines.append(f'| {link} | {table_cell(aliases)} | {item["kind"]} | {table_cell(item["summary"])} |')
    lines += ['', '本页由导出脚本生成；新增词条后重新生成即可，正文仍在各自的主要页面维护。', '']
    return '\n'.join(lines)


def export(source, output, refresh_index=False):
    if not source.is_dir():
        raise ValueError(f'原稿目录不存在：{source}')
    if source == output or source in output.parents or output in source.parents:
        raise ValueError('原稿与发布目录必须分开，且不能互相包含')
    notes = sorted(p for folder in FOLDERS for p in (source/folder).rglob('*.md')
                   if p.stem != INDEX)
    records = []
    targets = {'知识库/中文生信Wiki/00-开始写作': 'Home', '00-开始写作': 'Home', INDEX: INDEX,
               f'学习路线/{INDEX}': INDEX, f'知识库/中文生信Wiki/学习路线/{INDEX}': INDEX}
    names = set()
    alias_owners = {}
    for note in notes:
        if not note.resolve().is_relative_to(source):
            raise ValueError(f'发布范围外的文件链接：{note}')
        if note.stem in names or note.stem in RESERVED:
            raise ValueError(f'重复或保留名称：{note.stem}')
        if re.search(r'[\\/:*?"<>|]', note.stem):
            raise ValueError(f'Wiki 页面名含不支持字符：{note.stem}')
        names.add(note.stem)
        metadata, body = read_note(note.read_text())
        if metadata.get('status') in {'草稿', '待编写'}:
            names.remove(note.stem)
            continue
        if metadata.get('status') not in {'可用', '待核验', '待更新'}:
            raise ValueError(f'{note.name}: status 需使用草稿、可用、待核验或待更新')
        if not body.startswith('# ') or not isinstance(metadata.get('summary'), str) or not metadata['summary'].strip():
            raise ValueError(f'{note.name}: 缺少正文标题或 summary 属性')
        if not metadata.get('updated') or not metadata.get('validation'):
            raise ValueError(f'{note.name}: 请填写 updated 与 validation 属性')
        relative = note.relative_to(source).with_suffix('').as_posix()
        for key in (note.stem, relative, f'知识库/中文生信Wiki/{relative}'):
            targets[key] = note.stem
        item = {'title': note.stem, 'source': relative+'.md', 'folder': note.relative_to(source).parts[0],
                'kind': metadata.get('article_kind', ''), 'status': metadata['status'], 'aliases': metadata.get('aliases', []),
                'summary': metadata['summary'], 'validation': metadata.get('validation', ''),
                'updated': metadata.get('updated', ''), 'body': note.read_text(), 'path': note}
        records.append(item)
        for alias in item['aliases']:
            alias_owners.setdefault(alias, set()).add(note.stem)
    for alias, owners in alias_owners.items():
        if alias not in targets and len(owners) == 1:
            targets[alias] = next(iter(owners))
    rendered = {}
    for item in records:
        text = render_note(item['body'], item['path'], targets)
        # Keep execution status visible after removing private/editorial properties.
        first, rest = text.split('\n', 1)
        stamp = f'> 更新：{item["updated"]}；状态：{item["status"]}；验证范围：{item["validation"]}。'
        rendered[item['title']+'.md'] = first+'\n\n'+stamp+'\n'+rest
    rendered[INDEX+'.md'] = build_index(records)
    home = ['# 中文生信 Wiki', '',
            '持续完善的中文生物信息学知识库，系统连接概念与理论、方法与算法、工具与数据库、分析工作流及可复核案例。', '',
            '## 专题导航与检索', '',
            f'- [中文生信学习地图]({page_url("中文生信学习地图")})：了解前置知识和阅读顺序。',
            f'- [方法与算法系统学习指南]({page_url("方法与算法系统学习指南")})：按七条路线理解原理、手算、动手核对和排错。',
            f'- [全类型变异检测总览与证据地图]({page_url("全类型变异检测总览与证据地图")})：按输入、平台、倍性和类型选择检测分支。',
            f'- [变异检测分型与提取学习路线]({page_url("变异检测分型与提取学习路线")})：小变异、SV、CNV/PAV、重复及按类型/区域/样本提取。',
            f'- [变异注释与功能解释学习路线]({page_url("变异注释与功能解释学习路线")})：基因模型、后果提取和候选证据。',
            f'- [生信缺口清单与维护优先级]({page_url("生信缺口清单与维护优先级")})：查看仍需完善的真实数据与专业软件闭环。',
            f'- [R与Python生信编程入门路线]({page_url("R与Python生信编程入门路线")})：从运行脚本、样本表和缺失开始。',
            f'- [bulkRNA-seq分析与实践路线]({page_url("bulkRNA-seq分析与实践路线")})：从真实计数学习配对、DESeq2和结果诊断。',
            f'- [组学分析方法学习路线]({page_url("组学分析方法学习路线")})：逐类理解模型、算法与证据。',
            f'- [群体遗传学分析学习路线]({page_url("群体遗传学分析学习路线")})：多样性、分化、祖源、近交、选择与历史。',
            f'- [测序技术与PCR学习路线]({page_url("测序技术与PCR学习路线")})：PCR、各代测序、建库和质控。',
            f'- [统计学与组学验证学习路线]({page_url("统计学与组学验证学习路线")})：从基础检验进入跨组学验证。',
            f'- [生信基础算法学习路线]({page_url("生信基础算法学习路线")})：查找、图、动态规划、聚类与优化。',
            f'- [群体分析与GWAS学习路线]({page_url("群体分析与GWAS学习路线")})：沿熟悉的主线进入分析。',
            f'- [GWAS软件选择与比较]({page_url("GWAS软件选择与比较")})：按场景选模型，进入各软件详解。',
            f'- [GS算法软件与实践路线]({page_url("GS算法软件与实践路线")})：算法、六种软件和候选预测实践。',
            f'- [组学分析总览与选型]({page_url("组学分析总览与选型")})：按测量对象选择组学流程。',
            f'- [关键词索引]({page_url(INDEX)})：搜索中文名、英文名和缩写。', '',
            '## 实践与验证案例', '',
            f'- [人工长读长贯通SV调用与评估]({page_url("人工长读长贯通SV调用与评估")})：实际Sniffles/Truvari、群体缺失与格式问题。',
            f'- [官方TRGT小数据贯通重复分型]({page_url("官方TRGT小数据贯通重复分型")})：固定官方微型输入、长度和motif字段。',
            f'- [人工VCF贯通变异分类与提取]({page_url("人工VCF贯通变异分类与提取")})：bcftools分类、标准化、子集、掩蔽与人工小变异调用。',
            f'- [人工基因贯通区间注释与后果提取]({page_url("人工基因贯通区间注释与后果提取")})：SnpEff序列检查、bedtools和多ALT后果。',
            f'- [从FASTQ到比对检查与变异调用]({page_url("从FASTQ到比对检查与变异调用")})：minimap2、samtools和bcftools人工实测。',
            f'- [最小注释工作流与Snakemake恢复]({page_url("最小注释工作流与Snakemake恢复")})：首次运行、无变化重跑与缺失交付表恢复。',
            f'- [同一份样本表贯通R与Python]({page_url("同一份样本表贯通R与Python")})：八个样本、按ID合并、保留缺失和核对均值。',
            f'- [公开airway数据贯通DESeq2]({page_url("公开airway数据贯通DESeq2")})：原生R、真实公开计数、apeglm和效应阈值。',
            f'- [人工数据贯通群体遗传统计]({page_url("人工数据贯通群体遗传统计")})：多样性、描述性分化、LD与中性漂变。',
            f'- [人工区间贯通测序指标]({page_url("人工区间贯通测序指标")})：重复、双端重叠、覆盖、Phred和N50。',
            f'- [合成数据贯通统计推断]({page_url("合成数据贯通统计推断")})：基础检验、重采样、FDR与供体反例。',
            f'- [基础算法可复现练习]({page_url("基础算法可复现练习")})：六类人工算法与结果核验。',
            f'- [合成数据贯通PCA与GWAS]({page_url("合成数据贯通PCA与GWAS")})：120 样本、600 位点的实测示例。',
            f'- [公开小鼠数据贯通多软件GWAS]({page_url("公开小鼠数据贯通多软件GWAS")})：四个软件、七种配置的公开数据案例。',
            f'- [合成家系数据贯通GS]({page_url("合成家系数据贯通GS")})：GBLUP、家系验证与60个候选预测。',
            f'- [BLAST原理与种子延伸]({page_url("BLAST原理与种子延伸")})：用人工字符串理解启发式与精确比对。', '',
            '其他流程分别提供官方公开数据和实践入口，并注明文档核对与实际运行的区别。', '']
    # Curated routes are optional when exporting a smaller authored collection.
    home = [line for line in home if not line.startswith('- [') or
            any(page_url(name) in line for name in names | {INDEX})]
    sidebar = [f'[首页]({page_url("Home")}) · [关键词索引]({page_url(INDEX)})', '']
    if "方法与算法系统学习指南" in names:
        sidebar += [f'[方法系统学习]({page_url("方法与算法系统学习指南")})', ""]
    for folder in FOLDERS:
        subset = [r for r in records if r['folder'] == folder]
        home += [f'## {folder}（{len(subset)}）', '']
        sidebar += [f'### {folder}', '']
        for item in subset:
            link = f'[{item["title"]}]({page_url(item["title"])})'
            home.append(f'- {link}：{item["summary"]}')
            sidebar.append(f'- {link}')
        home.append(''); sidebar.append('')
    home += ['## 编写与维护', '',
             '[项目、代码与发布说明](https://github.com/SpiralKernel/bioinfo-wiki)。原稿在 Obsidian 维护，发布版自动生成导航与关键词索引。', '']
    rendered['Home.md'] = '\n'.join(home)
    rendered['_Sidebar.md'] = '\n'.join(sidebar)
    payload = {'schema': 1, 'generated_files': sorted(rendered), 'pages': [
        {k:v for k,v in item.items() if k not in {'body', 'path'}} for item in records]}
    payload['file_hashes'] = {name: hashlib.sha256(text.encode()).hexdigest()
                              for name, text in rendered.items()}
    for item in payload['pages']:
        item['sha256'] = payload['file_hashes'][item['title']+'.md']
    # Complete validation before updating either output or the generated source index.
    if (output/MANIFEST).is_symlink():
        raise ValueError('导出清单不能是符号链接')
    for name in rendered:
        if (output/name).is_symlink():
            raise ValueError(f'导出目标不能是符号链接：{name}')
    previous = json.loads((output/MANIFEST).read_text()) if (output/MANIFEST).exists() else {}
    for name in previous.get('generated_files', []):
        if Path(name).name != name or not name.endswith('.md'):
            raise ValueError('旧导出清单包含不安全的文件名')
    with tempfile.TemporaryDirectory() as tmp:
        staging = Path(tmp)
        for name, text in rendered.items():
            (staging/name).write_text(text)
        (staging/MANIFEST).write_text(json.dumps(payload, ensure_ascii=False, indent=2)+'\n')
        output.mkdir(parents=True, exist_ok=True)
        for name in previous.get('generated_files', []):
            if name not in rendered:
                (output/name).unlink(missing_ok=True)
        for file in staging.iterdir():
            shutil.copy2(file, output/file.name)
    if refresh_index:
        target = source/'学习路线'/f'{INDEX}.md'
        target.parent.mkdir(exist_ok=True)
        properties = ('---\ntags: ["知识库", "领域/生物信息学"]\ntype: 索引\narticle_kind: 学习路线\n'
                      'status: 可用\naliases: ["术语索引", "关键词检索"]\n'
                      f'updated: {date.today().isoformat()}\n---\n\n')
        target.write_text(properties+build_index(records, obsidian=True))
    print(f'已导出 {len(records)} 篇正文、关键词索引、首页和侧边栏。')
    return payload


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1]/'wiki')
    parser.add_argument('--refresh-index', action='store_true', help='同时刷新 Obsidian 中自动生成的关键词索引')
    args = parser.parse_args()
    try:
        export(args.source.resolve(), args.output.resolve(), args.refresh_index)
    except (ValueError, OSError) as error:
        parser.exit(1, f'导出失败：{error}\n')
