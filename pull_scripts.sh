#!/bin/bash

CLUSTER="scuhpc"
REMOTE_DIR="/mnt/may1nov1/u5023/pal/"
LOCAL_DIR="$HOME/work/NGS_analysis/"

mkdir -p "$LOCAL_DIR"
echo "===== 传输开始于: $(date) ====="

echo "正在扫描远程.sh文件..."
ssh "$CLUSTER" "find $REMOTE_DIR -type f -name '*.sh'" > filelist.tmp 2>/dev/null

if [ ! -s "filelist.tmp" ]; then
    echo "错误: 未找到任何.sh文件"
    rm -f filelist.tmp
    exit 1
fi

TOTAL_FILES=$(wc -l < filelist.tmp)
echo "找到 $TOTAL_FILES 个.sh文件"

COUNTER=0
SUCCESS=0
FAILED=0

while read -r remote_file; do
    ((COUNTER++))
    local_path="$LOCAL_DIR${remote_file#$REMOTE_DIR}"
    mkdir -p "$(dirname "$local_path")"
    echo "[$COUNTER/$TOTAL_FILES] 正在传输: $remote_file"
    if scp -q "$CLUSTER:$remote_file" "$local_path" 2>/dev/null; then
        echo "传输成功: ${remote_file#$REMOTE_DIR}"
        ((SUCCESS++))
    else
        echo "传输失败: ${remote_file#$REMOTE_DIR}"
        ((FAILED++))
    fi
done < filelist.tmp

rm -f filelist.tmp

echo "===== 传输完成于: $(date) ====="
echo "总计: $TOTAL_FILES 个文件"
echo "成功: $SUCCESS 个"
echo "失败: $FAILED 个"
echo "文件已保存到: $LOCAL_DIR"