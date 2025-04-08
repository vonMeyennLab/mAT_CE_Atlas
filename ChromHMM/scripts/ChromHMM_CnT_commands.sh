##################################################
Concatenated design
##################################################

### AdipoCRE

# Binarize Bam
sbatch --time=4:00:00 --ntasks=4 --mem-per-cpu=6000 --wrap="java -mx20000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
BinarizeBam \
-paired /cluster/work/nme/software/apps/ChromHMM/1.22/CHROMSIZES/mm10_nochr.txt \
/scratch/../ChromHMM_design_concatenated_AdipoCre.txt \
/scratch/../concatenated_AdipoCre/binary_files" \
/

# Learn Model
sbatch --time=24:00:00 --ntasks=20 --mem-per-cpu=5000 --wrap="seq 1 30 | parallel java -mx4000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
LearnModel \
-s 123 -nobrowser -noimage -nobed -noenrich -noautoopen \
/scratch/../concatenated_AdipoCre/binary_files/ \
/scratch/../concatenated_AdipoCre/model/ \
{} \
GRCm38" \
/

# Compare Models
sbatch --time=4:00:00 --ntasks=10 --mem-per-cpu=2000 --wrap="java -mx20000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
CompareModels \
/scratch/../concatenated_AdipoCre/model/emissions_30.txt \
/scratch/../concatenated_AdipoCre/model/ \
/scratch/../concatenated_AdipoCre/compare_models/concatenated_AdipoCre_1_to_30" \
/

# Make Segmentation
sbatch --time=24:00:00 --ntasks=10 --mem-per-cpu=6000 --wrap="ls /scratch/../concatenated_AdipoCre/model/model_9.txt | \
parallel java -mx5000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
MakeSegmentation \
{} \
/scratch/../concatenated_AdipoCre/binary_files \
/scratch/../concatenated_AdipoCre/segmentation" \
/

# Genomic Features Enrichment
sbatch --time=4:00:00 --ntasks=10 -R --mem-per-cpu=2000 --wrap="java -mx4000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
OverlapEnrichment \
/scratch/../concatenated_AdipoCre/segmentation/CE_DT_9_segments.bed \
/scratch/../ChromHMM/annotatr_mm/ \
/scratch/../concatenated_AdipoCre/annotatr_enrichment/CE_DT_genomicFeature_enrichment" \
/

# ENCODE cCRE Enrichment
sbatch --time=4:00:00 --ntasks=10 -R --mem-per-cpu=2000 --wrap="java -mx4000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
OverlapEnrichment \
/scratch/../concatenated_AdipoCre/segmentation/CE_DT_9_segments.bed \
/scratch/../ChromHMM/ENCODE_mm_collapsed/ \
/scratch/../concatenated_AdipoCre/cCRE_enrichment/CE_DT_cCRE_enrichment" \
/


####################################################################################################

### Ucp1ERCre

# Binarize Bam
sbatch --time=4:00:00 --ntasks=4 --mem-per-cpu=6000 --wrap="java -mx20000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
BinarizeBam \
-paired /cluster/work/nme/software/apps/ChromHMM/1.22/CHROMSIZES/mm10_nochr.txt \
/scratch/../ChromHMM_design_concatenated_Ucp1ERCre.txt \
/scratch/../concatenated_Ucp1ERCre/binary_files" \
/

# Learn Model
sbatch --time=24:00:00 --ntasks=20 --mem-per-cpu=5000 --wrap="seq 1 30 | parallel java -mx4000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
LearnModel \
-s 123 -nobrowser -noimage -nobed -noenrich -noautoopen \
/scratch/../concatenated_Ucp1ERCre/binary_files/ \
/scratch/../concatenated_Ucp1ERCre/model/ \
{} \
GRCm38" \
/

# Compare Models
sbatch --time=4:00:00 --ntasks=10 --mem-per-cpu=2000 --wrap="java -mx20000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
CompareModels \
/scratch/../concatenated_Ucp1ERCre/model/emissions_30.txt \
/scratch/../concatenated_Ucp1ERCre/model/ \
/scratch/../concatenated_Ucp1ERCre/compare_models/concatenated_Ucp1ERCre_1_to_30" \
/

# Make Segmentation
sbatch --time=24:00:00 --ntasks=10 --mem-per-cpu=6000 --wrap="ls /scratch/../concatenated_Ucp1ERCre/model/model_9.txt | \
parallel java -mx5000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
MakeSegmentation \
{} \
/scratch/../concatenated_Ucp1ERCre/binary_files \
/scratch/../concatenated_Ucp1ERCre/segmentation" \
/

# Genomic Features Enrichment
sbatch --time=4:00:00 --ntasks=10 -R --mem-per-cpu=2000 --wrap="java -mx4000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
OverlapEnrichment \
/scratch/../concatenated_Ucp1ERCre/segmentation/CETN_9_segments.bed \
/scratch/../ChromHMM/annotatr_mm/ \
/scratch/../concatenated_Ucp1ERCre/annotatr_enrichment/CETN_genomicFeature_enrichment" \
/

# ENCODE cCRE Enrichment
sbatch --time=4:00:00 --ntasks=10 -R --mem-per-cpu=2000 --wrap="java -mx4000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
OverlapEnrichment \
/scratch/../concatenated_Ucp1ERCre/segmentation/CETN_9_segments.bed \
/scratch/../ChromHMM/ENCODE_mm_collapsed/ \
/scratch/../concatenated_Ucp1ERCre/cCRE_enrichment/CETN_cCRE_enrichment" \
/