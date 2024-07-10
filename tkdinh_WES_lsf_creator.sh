#!/bin/bash
sample_list=('scA_4001-01_WES' 'scA_4004-01_WES' 'scA_4005-01_WES' 'scA_4006-01_WES' 'scA_4008-01_WES')

for sample in "${sample_list[@]}";
do
  echo "#BSUB -J ${sample}" >> ${sample}.lsf
	echo "#BSUB -W 15:00" >> ${sample}.lsf
	echo "#BSUB -o /rsrch6/scratch/lym_myl_rsch/tkdinh/WES/log" >> ${sample}.lsf
        echo "#BSUB -e /rsrch6/scratch/lym_myl_rsch/tkdinh/WES/stderr" >> ${sample}.lsf
        echo "#BSUB -cwd /rsrch6/scratch/lym_myl_rsch/tkdinh/WES" >> ${sample}.lsf
	echo "#BSUB -q medium" >> ${sample}.lsf
	echo "#BSUB -n 12" >> ${sample}.lsf
	echo "#BSUB -M 96G" >> ${sample}.lsf
	echo "#BSUB -R rusage[mem=96G]" >> ${sample}.lsf
	echo "#BSUB -u tkdinh@mdanderson.org" >> ${sample}.lsf
	echo "module load annovar/2019.10.24" >> ${sample}.lsf
	echo "module load R" >> ${sample}.lsf

	echo "python /rsrch6/scratch/lym_myl_rsch/tkdinh/WES/scripts/run_preprocessing.py \\" >> ${sample}.lsf
	echo "--threads 12 \\" >> ${sample}.lsf
	echo "--sequence \"WES\" \\" >> ${sample}.lsf
	echo "--input \"/rsrch6/scratch/lym_myl_rsch/tkdinh/WES/input\" \\" >> ${sample}.lsf
	echo "--sample \"$sample\" \\" >> ${sample}.lsf
	echo "--input1 \"/rsrch4/home/lym_myl_rsch/Green_Lab_NGS2/DLBCL_scAtlas/Raw_data/MedGenomeP2006445/FASTQ/${sample}_R1.fastq.gz\" \\" >> ${sample}.lsf
	echo "--input2 \"/rsrch4/home/lym_myl_rsch/Green_Lab_NGS2/DLBCL_scAtlas/Raw_data/MedGenomeP2006445/FASTQ/${sample}_R2.fastq.gz\" \\" >> ${sample}.lsf
	echo "--wes_type \"Agilent\" \\" >> ${sample}.lsf
	echo "--wes_model \"human\" \\"  >> ${sample}.lsf
	echo "--output \"/rsrch6/scratch/lym_myl_rsch/tkdinh/WES/output\" \\" >> ${sample}.lsf
done
