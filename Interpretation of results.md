# Interpretation of results


In this tutorial, it is our goal to aid clinitians in how they can interpret the results produced by NanoCAH.


## The output files

After you have performed the preperatory steps, like producing the unfiltered alignment file, and have run this alignment through the NanoCAH software; you should have produced the following files;

* A filtered and phased alignment file(_NanoCAH_filtered_alignment_sorted.bam)
* An assortment of aligned assemblies (_asm_Hap1_sorted.bam, _asm_Hap2_sorted.bam, asm_all_sorted.bam)
* A variant file containing small variants (_SNV.vcf.gz)
* A variant file containing larger variants (_SV.vcf.gz)


The variant files can be run through software like VarSeq to extrapolate relevant pathogenic variants.


## Interpretation of the structure

In most commonly used human reference sequences, like hg38 or T2T-chm13 there is one copy of the active CYP21A2 and one copy of the inactive pseudogene CYP21A1P. Unfortunately, as is explained in the introduction, high genomic instability of the region can cause many types of structural variants. To get a better understanding of the structural variants of the two haplotypes for a particular sample the alignment files should be evaluated using IGV. 

In this following section we will show you real-world examples of results obtained from clinically diagnosed CAH patients, showing outputs from NanoCAH and comparing them to the results obtained using the classical method of simple alignment filtering and phasing. Here we will show how using NanoCAH can be especially powerful in cases where pathogenetic variants are caused by or come in conjuction with large structural variants. 


#### A deletion including a part fo CYP21A2P, C4B and CYP21A2

In the case of a large deletion in one or more of the haplotypes, you will see a gap in coverage. When evaluating the alignment results obtained without NanoCAH (top picture), you see this aforementioned gap in coverage. If variants are called on this type of alignment, you can often run into problems where pathogenic variants are not called due to low coverage. If you on the otherhand look at the result after using NanoCAH (bottom picture), there is high coverage on both pseudogene and active gene, which avoids this problem. The way this is done, is that each side of the same haplotype, is split into a left and a right side, when there is a large deletion or insertion. NanoCAH then retains the reads that map to both sides of the same haplotype. Each side is named HapX_left and HapX_right, or Hap1_left and Hap1_right in this case. You should consider the two sides as overlapping each other, but only one side represents the correct alignment. In this case we see a lot of SNPs on the active gene side of the haplotype (Hap1_right). This indicates that the correct alignment is more likely on Hap1_left, and the deletion is therefore likely overlapping with the active gene. Both the assembly and SV variant call files output by NanoCAH show that the deletion has breakpoints in both the pseudogene and active gene, explaining the pathogenic variant of Hap1 in this case. Hap2 shows a normal structural variant with one pseudogene and one active gene, but the active gene contains a pathogenic SNV (not shown in detail). Therefore, it is in this case possible to explain the pathogenic variants of each haplotype and therefore give an exact genetic diagnosis to this patient.

![Case with deletion, without NanoCAH](./PLOTS/deletion_case_results_without_NanoCAH.png)

![Case with deletion, with NanoCAH](./PLOTS/deletion_case_results_with_NanoCAH.png)

#### An insertion

In the case of a large insetion in one or more of the haplotypes, you will often see inconsistencies in the SNPs of the haplotypes and an increased coverage around the pseudogene and active gene. The insertion/duplication is therefore most of the time missed by variant callers. When using NanoCAH, as mentioned for a deletion, there will be a left and right side of the haplotype, but it will be clear that an insertion is present. Looking at Hap1 in the example shown above the both the left and right side have reads mapping to the active gene and the pseudogene. In this case it seems that there are two pseudogenes and one active gene present. The active gene in Hap1 is found to have a pathogenic variant (SNV not shown in details), which explains the pathology of Hap1. Insertions/duplications can also appear with what has been called a "within"-block, as is seen for Hap2 in this case. When there is a "within"-block, then there are reads which we know are connected to a haplotype, in this case Hap2, but could not be placed during processing. This can happen, as is the case here, if there are two copies of the pseudogene that are very homologous. In the case of Hap2 there are three copies of the pseudogene and no active gene.

![Case with insertion, without NanoCAH](./PLOTS/insertion_case_results_without_NanoCAH.png)

![Case with insertion, with NanoCAH](./PLOTS/insertion_case_results_with_NanoCAH.png)

#### A chimera

The last example is a chimera. In this patient, one of the haplotypes has one copy of the pseudogene and a pseudogene/gene chimera. When using the standard method of filtering and phasing, the problem you see is that the alignment algorithm favors read alignment to the pseudogene, which makes the read coverage around the pseudogene high, but the read coverage around the active gene is low – maybe resulting in a deletion call. Like in the case of a true deletion, the low coverage will often prevent variant calling and any pathogenic variant call can therefore be missed. When using NanoCAH, the result shows a haplotype with continuous coverage and it is now clear that Hap2 has one pseudogene and one pseudogene/gene chimera, which is inactive. In this patient, Hap1 has a pathogenic SNV in the active gene (not shown in detail).

![Case with chimera, without NanoCAH](./PLOTS/chimera_case_results_without_NanoCAH.png)

![Case with chimera, with NanoCAH](./PLOTS/chimera_case_results_with_NanoCAH.png)



