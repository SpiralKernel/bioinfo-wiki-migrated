#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
wiki_checkout="$(dirname -- "$project_dir")/bioinfo-wiki.wiki"
main_remote="$(git -C "$project_dir" remote get-url origin)"
if [[ "$main_remote" != https://github.com/*/*.git ]]; then
    echo '主仓库 origin 需为 GitHub HTTPS .git 地址，请先核对远程配置。' >&2
    exit 1
fi
project_url="${main_remote%.git}"
wiki_remote="${project_url}.wiki.git"
wiki_url="${project_url}/wiki"
git_command=(git)
if command -v gh >/dev/null 2>&1; then
    git_command+=(-c 'credential.helper=!gh auth git-credential')
fi

if [[ ! -f "$project_dir/wiki/Home.md" ]]; then
    echo '请先运行 export_wiki.py 导出文章。' >&2
    exit 1
fi

if [[ ! -e "$wiki_checkout" ]]; then
    if ! "${git_command[@]}" clone "$wiki_remote" "$wiki_checkout"; then
        echo "请先打开 $wiki_url 并创建保存第一页，然后重试。" >&2
        exit 1
    fi
fi

if [[ ! -d "$wiki_checkout/.git" ]]; then
    echo "目标目录不是独立的 Wiki Git 仓库：$wiki_checkout" >&2
    exit 1
fi
if [[ "$("${git_command[@]}" -C "$wiki_checkout" remote get-url origin)" != "$wiki_remote" ]]; then
    echo '目标目录的远程地址与本项目不符，已停止。' >&2
    exit 1
fi
if [[ -n "$("${git_command[@]}" -C "$wiki_checkout" status --porcelain)" ]]; then
    echo 'Wiki 工作区有未提交修改，请先处理后再发布。' >&2
    exit 1
fi

"${git_command[@]}" -C "$wiki_checkout" pull --ff-only
python3 "$project_dir/scripts/sync_wiki_files.py" \
    --source "$project_dir/wiki" --destination "$wiki_checkout"
"${git_command[@]}" -C "$wiki_checkout" add -- .
if "${git_command[@]}" -C "$wiki_checkout" diff --cached --quiet; then
    "${git_command[@]}" -C "$wiki_checkout" push
    echo 'Wiki 文件没有变化，已同步待推送的提交。'
    exit 0
fi
"${git_command[@]}" -C "$wiki_checkout" commit -m '更新中文生信Wiki'
"${git_command[@]}" -C "$wiki_checkout" push
echo "已发布：$wiki_url"
