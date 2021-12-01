# CNVL
This repository holds the somatic copy number variant caller computation used in the research and generation of the following manuscript:

Farshidfar, et. al.  Integrative Molecular and Clinical Profiling of Acral Melanoma Identifies LZTR1 as a Key Tumor Promoter and Therapeutic Target.  https://www.biorxiv.org/content/10.1101/2021.04.20.440286v1

This caller was adapted from the computation used in [1] in 2013, developed by Siming Zhao and Murim Choi, and the original implementation of the p01CNV_segmentation code was developed by Siming Zhao in 2013. The computation has also been used in [2], [3], [4] and [5].

### Installation/Setup

The software consists of a main bash script CNVL, along with additional code files that were run in the following environment:
* g++ 10.2.0 compiler
* python 2.7.11, with the numpy package installed
* R 3.2.3, with the DNAcopy library installed
* samtools 1.8

It is likely that similar versions of each of these tools will be compatible, but that has not been tested.

The following steps will setup the software for use on a system:
1. Clone this git repository.
2. Cd into the CNVL repository and run the commands "g++ -o bamMetrics -static -O2 bamMetrics.cpp" and "chmod +x CNVL".
3. Ensure that python, R and samtools are all on the PATH environment variable.

### Running the software

CNVL was used to perform somatic calling of the exome data for the manuscript, and it takes the aligned bam files
for a paired tumor and normal sample, along with a bed file of target regions.  For the manuscript, the exome data
was aligned to the hs37d5 human reference using GATK 3 best practices, and the RefGene coding regions were used as
the target regions (the bed file for these regions are included in the repository).  Exome kit target regions can
also be used as the targets for the software.

The software expects that the tumor and normal bam files were generated using GATK best practices, and that the human
reference fasta file used to generate the reference is available.  Also, the "tumor purity" of the tumor sample
must be estimated prior to using the software, as it is one of the arguments to the CNVL command.

The CNVL command is the following:
```
CNVL ref.fasta targets.bed tumor.bam normal.bam purity outputPrefix
```
where "ref.fasta" is the fasta of the human reference, "targets.bed" is the bed file of the target regions,
"tumor.bam" and "normal.bam" are the BAM files for the tumor and normal samples, "purity" is the tumor purity
estimate (and should be a fraction between 0.0 and 1.0), and "outputPrefix" is the filename prefix given to
the CNVL output and intermediate files.

#### Output and Intermediate Files

Upon completion of the software, the main output file is "prefix.calls.txt", containing the CNV calls made by
the software.  The file is a tab-delimited file of the CNV calls, such as with this example:
```
chr     start   end     length  # markers       copy_ratio      copy_count      gainloss
chr1    60001   104300000       104240000       2096            0.702           1       loss
chr3    360001  88220000        87860000        1329            0.695           1       loss
```

The software also generates intermediate files used by the computation, that can be used for more detailed inspection of the calls:
* prefix_tumorMetrics.txt and prefix_normalMetrics.txt - Basic "exome" metrics for the tumor and normal samples
* prefix_tumorCov.txt and prefix_normalCov.txt - Per-target-region read depth information for the tumor and normal samples
* prefix_covRatio.txt - Normalized read depth ratios across the genome
* prefix_CBS_calling.txt - Raw results from the DNAcopy CBS computation
* prefix_cnvfull.txt - Final results for each of the CBS identified regions, including regions identified as copy-neutral

#### References

[1]  Zhao S, et al.  Landscape of somatic single-nucleotide and copy-number mutations in uterine serous carcinoma.  Proc. Natl. Acad. Sci. U. S. A. 110, 2916–2921 (2013).

[2] Zhao S, et al.  Mutational landscape of uterine and ovarian carcinosarcomas implicates histone genes in epithelial-mesenchymal transition.  Proc Natl Acad Sci U S A. 2016 Oct 25;113(43):12238-12243.

[3] Bi M, et al.  Genomic characterization of sarcomatoid transformation in clear cell renal cell carcinoma.  Proc Natl Acad Sci U S A. 2016 Feb 23;113(8):2170-5.
 
[4] Zhao S, et al.  Mutational landscape of uterine and ovarian carcinosarcomas implicates histone genes in epithelial-mesenchymal transition.  Proc Natl Acad Sci U S A. 2016 Oct 25;113(43):12238-12243.
 
[5] Choi J, et al.  Integrated mutational landscape analysis of uterine leiomyosarcomas.  Proc Natl Acad Sci U S A. 2021 Apr 13;118(15):e2025182118.
