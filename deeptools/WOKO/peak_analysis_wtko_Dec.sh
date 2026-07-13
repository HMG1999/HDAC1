#!/usr/bin/env bash
set -euo pipefail

# =========================
# 输入文件
# =========================
BW1="/mnt/My_disk/Li_linghang/liumeng_HDCAC1/dataset/output/proj_d19_wtko/GroupBigWigs/Sample_Clusters2/new/D19WT1_Dec.bw"
BW2="/mnt/My_disk/Li_linghang/liumeng_HDCAC1/dataset/output/proj_d19_wtko/GroupBigWigs/Sample_Clusters2/new/D19KO_Dec.bw"

BED1="D19WT1_Dec_peaks.bed"
BED2="D19KO_Dec_peaks.bed"

TSS_SRC="/mnt/My_disk/Li_linghang/liumeng_HDCAC1/dataset/output/proj_d19_wtko/GroupBigWigs/Sample_Clusters2/mm10_TSS.bed"
TSS_LOCAL="mm10_TSS.bed"
BEDTOOLS="/home/lilinghang/anaconda3/envs/deeptools/bin/bedtools"

# 过滤掉 signal_sum 最低 20% peak
LOW_QUANTILE=0.20

LFC_CUTOFF=1
PSEUDO=0.01

# =========================
# 1) union peaks
# =========================
cat "$BED1" "$BED2" \
  | cut -f1-3 \
  | sort -k1,1V -k2,2n \
  | "${BEDTOOLS}" merge -i - > union_peaks.bed

# =========================
# 2) 提取 bigwig 信号
# =========================
multiBigwigSummary BED-file \
  --bwfiles "$BW1" "$BW2" \
  --BED union_peaks.bed \
  --outFileName peak_signal.npz \
  --outRawCounts peak_signal.tab

# =========================
# 3) 自动阈值 + 分组
# =========================
python - <<PY
import pandas as pd
import numpy as np

low_q = float("${LOW_QUANTILE}")
lfc_cut = float("${LFC_CUTOFF}")
pseudo = float("${PSEUDO}")

df = pd.read_csv("peak_signal.tab", sep="\t")
c1, c2 = df.columns[3], df.columns[4]

signal_sum = df[c1] + df[c2]
thr = np.quantile(signal_sum, low_q)
print(f"[INFO] Auto threshold (quantile={low_q}): {thr:.6f}")

keep = signal_sum >= thr
sub = df.loc[keep].copy()

sub["log2FC"] = np.log2((sub[c2] + pseudo) / (sub[c1] + pseudo))
sub["Group"] = np.where(sub["log2FC"] >= lfc_cut, "Gained",
                 np.where(sub["log2FC"] <= -lfc_cut, "Lost", "Unchanged"))

sub.to_csv("peak_signal_grouped.tsv", sep="\t", index=False)

vc = sub["Group"].value_counts()
print("[INFO] Kept peaks:", len(sub))
print("[INFO] Filtered peaks:", len(df)-len(sub))
print("[INFO] Group counts:")
print(vc.to_string())
PY

# =========================
# 4) 按类别拆分 BED
# =========================
awk 'NR>1 && $NF=="Gained"    {print $1,$2,$3}' OFS="\t" peak_signal_grouped.tsv > gained.bed
awk 'NR>1 && $NF=="Lost"      {print $1,$2,$3}' OFS="\t" peak_signal_grouped.tsv > lost.bed
awk 'NR>1 && $NF=="Unchanged" {print $1,$2,$3}' OFS="\t" peak_signal_grouped.tsv > unchanged.bed

echo "[INFO] line counts:"
wc -l union_peaks.bed gained.bed lost.bed unchanged.bed

# =========================
# 5) TSS profile (reference-point)
# =========================
cp "$TSS_SRC" "$TSS_LOCAL"

computeMatrix reference-point \
  -S "$BW1" "$BW2" \
  -R "$TSS_LOCAL" \
  -a 3000 -b 3000 \
  --referencePoint center \
  -o matrix_tss_refpoint.gz \
  -p 20 \
  --missingDataAsZero \
  --skipZeros

plotProfile -m matrix_tss_refpoint.gz \
  --perGroup \
  --plotTitle "ATAC-seq signal around TSS" \
  -out atac_tss_profile.png

# =========================
# 6) gained/lost/unchanged heatmap
# =========================
computeMatrix reference-point \
  -S "$BW1" "$BW2" \
  -R gained.bed lost.bed unchanged.bed \
  -a 3000 -b 3000 \
  -o d19_wt_ko_matrix.gz \
  --referencePoint center \
  -p 30 \
  --missingDataAsZero \
  --skipZeros

plotHeatmap \
  -m d19_wt_ko_matrix.gz \
  --regionsLabel "Gained" "Lost" "Unchanged" \
  --outFileName heatmap_D19_WT_vs_KO.pdf \
  --colorList "white,blue,red" \
  --zMin 0 --zMax auto

# =========================
# 7) scale-regions TSS heatmap
# =========================
computeMatrix scale-regions \
  -S "$BW1" "$BW2" \
  -R "$TSS_LOCAL" \
  --beforeRegionStartLength 3000 \
  --regionBodyLength 5000 \
  --afterRegionStartLength 3000 \
  --skipZeros \
  -o matrix_tss_scale.gz

plotHeatmap \
  -m matrix_tss_scale.gz \
  -out ExampleHeatmap1.png \
  --colorList "white,red"

echo "[INFO] Done."
