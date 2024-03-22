# ODZ-NOB
Code for analyzing oxygen deficient zone nitrite oxidizing bacteria


### Assembly and Binning with KBase
#### KBase Assembly: 
Each sample was assembled individually using MEGAHIT v1.2.9. The “advanced” parameters including minimum multiplicity, minimum k-mer size, maximum k-mer size, and increment k-mer size were left as the default. Parameter preset was set to meta-sensitive, and minimum contig length was set to 1000. 

#### KBase Binning:
Three binning programs were used on the assembly from each sample.
For CONCOCT v1.1, the read mapping tool was Bowtie2 very-sensitive, the minimum contig length was 1000, contig split size was 10,000, and the remaining parameters were left as default including contig split overlap, kmer length, maximum number of clusters and iterations for VGMM, and percent of total PCA used.
For MaxBin2 v2.2.4, the probability threshold was 0.8, the marker set was 107 Bacterial Marker Gene Set, and the minimum contain length was 1000.
For MetaBAT2 v1.7, the only controllable parameter on KBase is the minimum contain length (set to 1500, the minimum allowed), and all other parameters are forced to default.
To compare the three binning results for each sample and select the best possible MAGs, DAS-Tool v1.1.2 with diamond as the gene identification tool, a score threshold of 0.5, a duplicate penalty of 0.6, and a megabit penalty of 0.5.

#### KBase initial processing:
CheckM v1.0.18 was used to check the quality of each bin using default settings.

GTDB-Tk classify was used to identify the closest known taxonomy with the default minimum alignment percent of 10. 

FastANI v0.1.3 was used to calculate the average nucleotide identity of NOB MAGs and previously published NOB genomes and SAGs. Default settings were used.

A first look for nxrB genes was performed on KBase using the “Search with HMMs of Environmental Bioelement Families” v1 and the available Nitrite oxidation-NxrB HMM. Each result was checked manually using NCBI blastn searches. Further nxrB genes were identified using BLASTn v2.12.0 on KBase using known nxrB genes previously identified from NOB, an E-value threshold of 0.001, and a minimum sequence identity threshold of 50%.

### Mapping Relative Abundance of MAGs with Bowtie2
Short reads from OMZ metagenome samples produced in this study and published metagenomes from OMZ regions and the Tara Oceans database were mapped to ODZ NOB MAGs created in this study and published NOB genomes and SAGs using bowtie2. 

    #!/bin/bash

    #SBATCH --mem=500
    #SBATCH -t 14:00:00
    #SBATCH --array=1-2
    #SBATCH --mail-type=begin
    #SBATCH --mail-type=fail
    #SBATCH --mail-type=end
    #SBATCH --mail-user=sf3033@princeton.edu

    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p array_set1.txt)

    bowtie2 --end-to-end --very-sensitive -x all_NOB_one -1 /projects/WARD/Sam_Fortin/JGI_2019_seq/${SAMPLE}_QC_reads_R1.fastq -2   /projects/WARD/Sam_Fortin/JGI_2019_seq/${SAMPLE}_QC_reads_R2.fastq>${SAMPLE}_all_NOB_bowtie2.sam

    samtools view -S -b ${SAMPLE}_all_NOB_bowtie2.sam>${SAMPLE}_all_NOB.bam
    samtools sort ${SAMPLE}_all_NOB.bam -o ${SAMPLE}_all_NOB_sort.bam
    samtools index ${SAMPLE}_all_NOB_sort.bam
    samtools idxstats ${SAMPLE}_all_NOB_sort.bam | tee ${SAMPLE}_all_NOB_one.txt
