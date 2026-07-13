#!/usr/bin/env bash

# =========================
# 比较方向：log2FC = log2(BW2/BW1)
#   BW1/BAM1/BED1 = 对照/分母（reference）→ WT
#   BW2/BAM2/BED2 = 处理/分子（treatment）→ KO
#   Gained = KO相比WT上调, Lost = KO相比WT下调
# =========================
set -euo pipefail

# =========================
# 输入文件（D16_D19_H3K27Ac: D16 vs D19）
# =========================
BASE="/mnt/My_disk/Li_linghang/liumeng_HDCAC1/Cut_Tag/D16_D19_H3K27Ac"
BAM1="$BASE/D16_H3K27acAligned.sortedByCoord.out.bam"
BAM2="$BASE/D19_H3K27acAligned.sortedByCoord.out.bam"

BW1="D16_H3K27Ac.bw"
BW2="D19_H3K27Ac.bw"

BED1="D16_H3K27ac.bed"
BED2="D19_H3K27ac.bed"

OUTDIR="/mnt/My_disk/Li_linghang/liumeng_HDCAC1/Cut_Tag/cuttag_pipeline_1/H3K27Ac"
mkdir -p "$OUTDIR"/{beds,bigwig,matrices,heatmaps,tables}

# =========================
# 过滤参数
# =========================
LOW_QUANTILE=0.20
LFC_CUTOFFS="1"
PSEUDO=0.01
THREADS=20

# deepTools conda
DT="/home/lilinghang/anaconda3/envs/deeptools/bin"
ST="/home/lilinghang/software/samtools-1.21/bin/samtools"

# =========================
# 1) xls → BED
# =========================
echo "[INFO] xls -> bed"
for xls in "$BASE"/*_peaks.xls; do
  sample="$(basename "$xls" _peaks.xls)"
  bed="$OUTDIR/beds/${sample}.bed"
  awk -F'\t' '!/^#/ && $1 != "chr" && NF >= 3 {
    print $1"\t"$2"\t"$3"\t"$9"\t"$8
  }' "$xls" > "$bed"
  echo "[INFO]  $(basename $xls) -> $bed  ($(wc -l < "$bed") peaks)"
done

# =========================
# 2) union peaks
# =========================
echo "[INFO] merge peaks"
cat "$OUTDIR/beds/$BED1" "$OUTDIR/beds/$BED2" \
  | cut -f1-3 \
  | sort -k1,1V -k2,2n \
  | "$DT/bedtools" merge -i - > "$OUTDIR/beds/union_peaks.bed"
echo "[INFO] union peaks: $(wc -l < "$OUTDIR/beds/union_peaks.bed")"

# =========================
# 3) BAM index（如缺失）
# =========================
echo "[INFO] check BAM index"
for bam in "$BAM1" "$BAM2"; do
  if [[ ! -s "$bam.bai" ]]; then
    echo "[INFO] indexing $bam"
    "$ST" index -@ "$THREADS" -o "$bam.bai" "$bam"
  fi
done

# =========================
# 4) BAM → BigWig
# =========================
echo "[INFO] bamCoverage"
"$DT/bamCoverage" \
  --bam "$BAM1" \
  --outFileName "$OUTDIR/bigwig/$BW1" \
  --outFileFormat bigwig \
  --binSize 10 \
  --normalizeUsing RPGC \
  --effectiveGenomeSize 1870000000 \
  --extendReads 200 \
  --numberOfProcessors "$THREADS"

"$DT/bamCoverage" \
  --bam "$BAM2" \
  --outFileName "$OUTDIR/bigwig/$BW2" \
  --outFileFormat bigwig \
  --binSize 10 \
  --normalizeUsing RPGC \
  --effectiveGenomeSize 1870000000 \
  --extendReads 200 \
  --numberOfProcessors "$THREADS"

# =========================
# 5) 提取 union peaks 信号
# =========================
echo "[INFO] multiBigwigSummary"
"$DT/multiBigwigSummary" BED-file \
  --bwfiles "$OUTDIR/bigwig/$BW1" "$OUTDIR/bigwig/$BW2" \
  --BED "$OUTDIR/beds/union_peaks.bed" \
  --outFileName "$OUTDIR/matrices/peak_signal.npz" \
  --outRawCounts "$OUTDIR/tables/peak_signal.tab" \
  --numberOfProcessors "$THREADS"

# =========================
# 6) 整体信号热图
# =========================
echo "[INFO] overall heatmap"
"$DT/computeMatrix" reference-point \
  -S "$OUTDIR/bigwig/$BW1" "$OUTDIR/bigwig/$BW2" \
  -R "$OUTDIR/beds/union_peaks.bed" \
  -a 3000 -b 3000 \
  --referencePoint center \
  -o "$OUTDIR/matrices/overall_union_peaks.matrix.gz" \
  --missingDataAsZero --skipZeros \
  -p "$THREADS"

"$DT/plotHeatmap" \
  -m "$OUTDIR/matrices/overall_union_peaks.matrix.gz" \
  --outFileName "$OUTDIR/heatmaps/overall_union_peaks.white_to_153169.pdf" \
  --colorList "white,#3376CD" \
  --plotTitle "D16_D19_H3K27Ac: union peaks signal" \
  --regionsLabel "Union peaks" \
  --zMin 0 --zMax auto \
  --interpolationMethod bilinear \
  --sortRegions descend \
  --sortUsing mean \
  --missingDataColor white \
  --whatToShow "plot, heatmap and colorbar" \
  --heatmapWidth 4 --heatmapHeight 30 --dpi 300

# =========================
# 7) log2FC bigWig（KO/WT）
# =========================
echo "[INFO] bigwigCompare"
"$DT/bigwigCompare" \
  --bigwig1 "$OUTDIR/bigwig/$BW2" \
  --bigwig2 "$OUTDIR/bigwig/$BW1" \
  --operation log2 \
  --pseudocount "$PSEUDO" \
  --binSize 10 \
  --numberOfProcessors "$THREADS" \
  --outFileFormat bigwig \
  --outFileName "$OUTDIR/bigwig/D19_vs_D16.log2fc.pseudo${PSEUDO}.bw"

# =========================
# 8) 分类 + 热图
# =========================
echo "[INFO] classify and plot: cutoff=$LFC_CUTOFFS"

for cutoff in $LFC_CUTOFFS; do

  tag=$(python -c "x=float('$cutoff'); print(f'{x:.6g}'.replace('.','p'))")
  echo "[INFO]  cutoff=$cutoff (tag=$tag)"

  # 8a) 分类
  python - "$OUTDIR/tables/peak_signal.tab" "$OUTDIR/tables/peak_signal_grouped_${tag}.tsv" "$OUTDIR/beds" "$cutoff" "$LOW_QUANTILE" "$PSEUDO" <<'PY'
import sys, math, csv
from pathlib import Path
infile, out_tsv, beds_dir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
lfc_cut, low_q, pseudo = float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
tag = (f'{lfc_cut:.6g}').replace('.','p')
with infile.open(newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    rows = [r for r in reader if r]
values = [(float(r[3])+float(r[4]), float(r[3]), float(r[4]), r) for r in rows]
s = sorted(v[0] for v in values)
pos = low_q * (len(s)-1)
lo = int(pos); hi = min(lo+1, len(s)-1); frac = pos - lo
thr = s[lo]*(1-frac) + s[hi]*frac if len(s) > 1 else s[0]
cnt = {'Gained':0, 'Lost':0, 'Non sig':0}
bed_out = {}
for g in cnt:
    bn = beds_dir / (g.lower().replace(' ','_') + '_lfc_' + tag + '.bed')
    bed_out[g] = bn.open('w')
kept = []
try:
    for sm,d,n,r in values:
        if sm < thr: continue
        l2 = math.log2((n+pseudo)/(d+pseudo))
        g = 'Gained' if l2 >= lfc_cut else ('Lost' if l2 <= -lfc_cut else 'Non sig')
        cnt[g] += 1; bed_out[g].write('\t'.join(r[:3])+'\n')
        kept.append(r+[f'{sm:.12g}',f'{l2:.12g}',g])
finally:
    for h in bed_out.values(): h.close()
with out_tsv.open('w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(header+['signal_sum','log2FC','Group']); w.writerows(kept)
with (beds_dir.parent/'tables'/f'summary_lfc_{tag}.txt').open('w') as f:
    f.write(f'low_quantile\t{low_q}\nkept_peaks\t{len(kept)}\nfiltered_peaks\t{len(rows)-len(kept)}\n')
    for k in cnt: f.write(f'{k}\t{cnt[k]}\n')
PY

  echo "[INFO]   Gained:$(wc -l < $OUTDIR/beds/gained_lfc_${tag}.bed)  Lost:$(wc -l < $OUTDIR/beds/lost_lfc_${tag}.bed)  Non_sig:$(wc -l < $OUTDIR/beds/non_sig_lfc_${tag}.bed)"

  # 8b) signal groups heatmap
  "$DT/computeMatrix" reference-point \
    -S "$OUTDIR/bigwig/$BW1" "$OUTDIR/bigwig/$BW2" \
    -R "$OUTDIR/beds/gained_lfc_${tag}.bed" "$OUTDIR/beds/lost_lfc_${tag}.bed" "$OUTDIR/beds/non_sig_lfc_${tag}.bed" \
    -a 3000 -b 3000 --referencePoint center \
    -o "$OUTDIR/matrices/signal_groups_lfc_${tag}.matrix.gz" \
    --missingDataAsZero --skipZeros -p "$THREADS"

  "$DT/plotHeatmap" \
    -m "$OUTDIR/matrices/signal_groups_lfc_${tag}.matrix.gz" \
    --regionsLabel "Gained" "Lost" "Non sig" \
    --outFileName "$OUTDIR/heatmaps/signal_groups_lfc_${tag}.white_to_153169.pdf" \
    --colorList "white,#3376CD" \
    --plotTitle "D16_D19_H3K27Ac: signal groups |log2FC| >= ${cutoff}" \
    --zMin 0 --zMax auto --interpolationMethod bilinear \
    --sortRegions descend --sortUsing mean \
    --missingDataColor white --whatToShow "plot, heatmap and colorbar" \
    --heatmapWidth 4 --heatmapHeight 30 --dpi 300

  # 8c) 自定义：Gained vs Lost 红白热图
  echo "[INFO]  custom Gained vs Lost red heatmap"
  if [[ -s "$OUTDIR/beds/gained_lfc_${tag}.bed" && -s "$OUTDIR/beds/lost_lfc_${tag}.bed" ]]; then
    "$DT/computeMatrix" reference-point \
      -S "$OUTDIR/bigwig/$BW1" "$OUTDIR/bigwig/$BW2" \
      -R "$OUTDIR/beds/gained_lfc_${tag}.bed" "$OUTDIR/beds/lost_lfc_${tag}.bed" \
      -a 3000 -b 3000 --referencePoint center \
      -o "$OUTDIR/matrices/signal_gainlost_lfc_${tag}.matrix.gz" \
      --missingDataAsZero --skipZeros -p "$THREADS"
    "$DT/plotHeatmap" \
      -m "$OUTDIR/matrices/signal_gainlost_lfc_${tag}.matrix.gz" \
      --regionsLabel "Gained" "Lost" \
      --outFileName "$OUTDIR/heatmaps/signal_gainlost_lfc_${tag}.white_to_red.pdf" \
      --colorList "white,red" \
      --plotTitle "D16_D19_H3K27Ac: Gained vs Lost |log2FC| >= ${cutoff} (custom)" \
      --zMin 0 --zMax auto --interpolationMethod bilinear \
      --sortRegions descend --sortUsing mean \
      --missingDataColor white --whatToShow "plot, heatmap and colorbar" \
      --heatmapWidth 4 --heatmapHeight 24 --dpi 300
    echo "[INFO]  custom heatmap done"
  else
    echo "[WARN]  gained or lost BED empty, skip custom heatmap"
  fi

done

echo "[INFO] ALL DONE"
