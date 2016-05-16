#!/bin/sh
#
#$ -N GM12878
#$ -q KA
#$ -pe smp 8
#$ -cwd

#julia -p 8 --depwarn=no  ~/.julia/v0.4/IDPASE/src/prep_runs.jl -a /Shared/Au/bdeonovic/work/haplotype/gold_standard8/data/SRR1153470.hisat-grch37.unique.renamed.psl /Shared/Au/bdeonovic/work/haplotype/gold_standard8/data/sny.gmap-grch37-best.psl /Shared/Au/bdeonovic/work/haplotype/gold_standard8/data/ccs75.gmap-grch37-best.psl -g /Shared/Au/bdeonovic/work/haplotype/gold_standard8/data/refFlat_20140202.uniquenamed-20151222.newfmt.grch37.gpd -v /Shared/Au/yunhao/Project/Haplotype/Plot/Overlap/Samtools-1KG-Illumina.overlap.vcf -q /Shared/Au/jason/Data/UtahTrio/GM12878/IDP-ASE/Data/SRR1153470.hisat-hg19.unique.renamed.fastq /Shared/Au/data/SNYDER_Sept2015/data/fastq/SRR1163655.fastq /Shared/Au/jason/Data/UtahTrio/GM12878/IDP-ASE/Data/ccs75.fq -d /Shared/Au/bdeonovic/work/haplotype/gold_standard8/temp -c 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X -f 1 1 1 -o /Shared/Au/bdeonovic/work/haplotype/gold_standard8/test_in/out -i -s -r /Shared/Au/bdeonovic/work/haplotype/gold_standard8/test_out2/results/real_idp.txt -l 101
julia -p 8 --depwarn=no  ~/.julia/v0.4/IDPASE/src/prep_runs.jl -a /Shared/Au/bdeonovic/work/haplotype/gold_standard8/data/SRR1153470.hisat-grch37.unique.renamed.psl /Shared/Au/bdeonovic/work/haplotype/gold_standard8/data/sny.gmap-grch37-best.psl /Shared/Au/bdeonovic/work/haplotype/gold_standard8/data/ccs75.gmap-grch37-best.psl -g /Shared/Au/bdeonovic/work/haplotype/gold_standard8/data/idp.unique.renamed.newfmt.corrected.gpd -v /Shared/Au/bdeonovic/work/haplotype/gold_standard8/data/gm12878_biallelic.idp.exonic.10depth.sort.header.corrected.vcf -q /Shared/Au/jason/Data/UtahTrio/GM12878/IDP-ASE/Data/SRR1153470.hisat-hg19.unique.renamed.fastq /Shared/Au/data/SNYDER_Sept2015/data/fastq/SRR1163655.fastq /Shared/Au/jason/Data/UtahTrio/GM12878/IDP-ASE/Data/ccs75.fq -d /Shared/Au/bdeonovic/work/haplotype/gold_standard8/temp -c 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X -f 1 1 1 -o /Shared/Au/bdeonovic/work/haplotype/gold_standard8/test_in/out -i -r /Shared/Au/bdeonovic/work/haplotype/gold_standard8/test_out3/results/real_idp.txt -l 101


