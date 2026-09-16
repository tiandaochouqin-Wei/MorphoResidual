#!/bin/bash
set -e
CANCER=$1; PDC=$2; S=/public/home/fjhui/ZW/scripts
echo ">>> [$CANCER] aliquot"; python $S/build_aliquot_xwalk.py $PDC $CANCER
echo ">>> [$CANCER] slide_map"; python $S/build_slide_map.py $CANCER
PROT=$(ls /public/home/fjhui/ZW/$CANCER/omics/protein/*.tmt1*.tsv | head -1); echo ">>> 蛋白: $PROT"
echo ">>> [$CANCER] residual_analysis"
MORPHO_ROOT=/public/home/fjhui/ZW/$CANCER MORPHO_RNA_MANIFEST=$S/manifest_rna_tumor_$CANCER.tsv \
MORPHO_ALIQUOT_XWALK=$S/aliquot_to_case_tumor_$CANCER.tsv MORPHO_SLIDE_MAP=$S/slide_type_map_$CANCER.tsv \
MORPHO_PROTEIN_TSV=$PROT MORPHO_WSI_EMB_DIR=/public/home/fjhui/ZW/$CANCER/WSI/emb MORPHO_N_WORKERS=8 \
python -u $S/residual_analysis.py
