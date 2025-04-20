#!/bin/bash
# ATAC-seq flagstat 统计脚本（最终修正版）

# 设置报告目录
REPORT_DIR="/mnt/may1nov1/u5023/pal/atac-seq/report"

# 输出表头
printf "%-20s %-12s %-12s %-12s %-12s %-12s\n" \
       "Sample" "Total" "Mapped(%)" "Proper(%)" "Singletons(%)" "DiffChr(%)"
echo "========================================================================"

# 遍历所有 flagstat 文件
for file in ${REPORT_DIR}/*.flagstat.txt; do
    # 提取样本名
    sample=$(basename "$file" .flagstat.txt)
    
    # 解析关键统计值
    total=$(awk 'NR==1{print $1}' "$file")
    mapped_pct=$(grep "mapped (" "$file" | head -n 1 | grep -oP '\(\K[0-9.]+(?=%)')
    proper_pct=$(grep "properly paired" "$file" | grep -oP '\(\K[0-9.]+(?=%)')
    singleton_pct=$(grep "singletons (" "$file" | grep -oP '\(\K[0-9.]+(?=%)')
    diff_chr_pct=$(grep "different chr (mapQ>=5)" "$file" | grep -oP '\(\K[0-9.]+(?=%)')

    # 处理可能缺失的值
    [ -z "$mapped_pct" ] && mapped_pct="N/A"
    [ -z "$diff_chr_pct" ] && diff_chr_pct="N/A"
    
    # 输出格式化结果
    printf "%-20s %-12s %-12s %-12s %-12s %-12s\n" \
           "$sample" "$total" "$mapped_pct%" "$proper_pct%" "$singleton_pct%" "$diff_chr_pct%"
done
