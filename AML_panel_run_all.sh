# Run all analysis code for R1 130nt + R2 21nt paired-end fastq file

# Please change R1 and R2 fastq file names

# Triming
python adapter_trim_f60_v2.py --FastqR1 aml7_S9_L001_R1_001.fastq --FastqR2 aml7_S9_L001_R2_001.fastq --OutputName trim_Lib_7.fastq

# # build index, generates .bt2 files, don't need to run again
# bowtie2-build BDA_Leukemia_v20220717.fasta BDA_Leukemia

for i in 7 #1 2 3 4 5 6 7 8 9 10 11
do
	echo Alignment...
	bowtie2 -x BDA_Leukemia -U trim_Lib_$i.fastq -S alignlib$i.sam

	echo Sorting...
	python sort_index_v2.py --InputSam alignlib$i.sam --OutputBam alignlib$i.bam

	echo VariantCall...
	python UMI_counter3_Vote_0_dyAmplicon_20201109.py --AlignmentFile alignlib$i.bam --FastaFile BDA_Leukemia_v20220717.fasta --EnrichFile BDA_Leukemia_EnR_v20220717.txt --UMIlen 15
done

