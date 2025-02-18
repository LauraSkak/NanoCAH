# NanoCAH
A local rephaser build to deal with difficult to phase regions with segmental duplications using Oxford Nanopore long-read data.

### Introduction

This software was developed to solve the notoriously difficult gene region around CYP21A2, which cause congenital adrenal hyperplasia (CAH). The CYP21A2 gene is difficult to resolve due to the homologous pseudogene CYP21A1P located just upstream of the active CYP21A2 gene. The close proximity and high homology between the two often make read mapping to the region ambiguous. Additionally, high genomic instability in the area causes large insertions, deletions and chimeric genes. All these factors lead to variant calls from this region often being unreliable.
To solve this problem Oxford Nanopore adaptive sampling of chromosome 6 for each sample and the resulting data is run through the NanoCAH algorithm. By employing this specialized bioinformatic algorithm, it is not possible to accurately distinguish between reads that map to the active CYP21A2 gene and the reads that map to the pseudogene, determine the presence of large deletion, insertions or chimeric genes and phase reads, thereby making variant calling in the region much more reliable.

### The algorithm

The following is a more technical description of how the NanoCAH algorithm works.

The input parameters:

* maxdiff - The maximal SNP difference, used for merging blocks. (DEFAULT: 2)
* minimum_overlap - The minimal identical overlap between blocks, used for low variance merging. (DEFAULT: 3)
* minimum_count - The minimum SNP representation count for a SNP to remain usable. (DEFAULT: 2)
* minimum_depth - The minimum read count at a particular SNP position, is used for block splitting. (DEFAULT: 3)
* minimum_coverage - The minimum read count in one block, is used in read selection. (DEFAULT: 5)

##### Creating data dictionaries

All reads are first loaded into a read dictionary, which will contain the base at each position which have found at least one read that have a SNP that differ from the reference. A variant dictionary is created containing information from each position with base variantion between reads.

The variant dictionary is then filtered for low coverage SNPs. This means is a base is found in less than the minimum count reads, then the base is removed. If the variant now only contains reads that match the reference, then the variant is deem unusable and is removed from the variant dictionary.

... missing add a figure showing this initial removal of unusable variants ...

Each read is then added to a haploblock dictionary. If a read overlaps with at least one block with the minimum overlap, then the read is added to the largest block. If a read overlaps with no existing blocks, then a new block is created.

If a read contains no usable variants, it is classed as unphasable. 

Because a read that contain less SNPs than the minimum overlap will never overlap with an existing block; a round of no variance merging is performed where these one read blocks are added to the largest block with an identity of 100%. 

... missing add a picture showing this initial merging ...


##### Performing low variance merging

Once all phasable reads have been assigned to a block; the blocks are merged with increasing variance, iteration from 0 count difference to the maximal difference.
For example, if two blocks have a SNP difference of 3, meaning they have different bases at 3 positions, and at each of these 3 positions one block has a depth less than the minimum count, then the two blocks can be merged. For each mismatched position the reads from the block with the lowest coverage are temperarily removed.
If one of these 3 positions have a depth larger than the minimum count in both blocks, then that is defined as a high confidence SNP, thereby preventing the two blocks from being merged. 

During this low variance merging process addition variants are deemed unusable. This happens, if the temperary removal of read SNPs causes the base count at that position to be covered with less than the minimum count reads. If so, the base is removed. If this causes the variant to only contain a base that matches the reference and/or an indel, then the variant is deemed as unusable and removed.

... missing figure showing merging of two blocks with 3 differences ...

Rounds of ambiguous block cleanups moves reads that map to more than one block to the block with the highest coverage. Each block is then split, if two sections are covered by less than the minimum depth reads. After block splitting the blocks are merged again. This causes necessary shuffling in regions with low variant counts. It also classes reads that during low variance merging got its SNPs removed.

Lastly the reads that map to more than one block is marked as ungrouped and blocks are split, if there are neighboring SNPs with no reads in common.

##### Read selection

After low variance merging each block should only contain reads that uniquely map to one block. Now we utilize that some reads map more than once to select the best mapping and identify blocks which are connected as an insertion, deletion or chimera. This is done by evaluating each block and determining if the secondary reads mapped to that block should be kept or removed. Here there are 2 scenarios;

* keeping both the primary and secondary reads with multiple mapping will cause another block's coverage to fall below the minimum coverage, and therefore that other block should be removed. 

* Two blocks have more than the minimum coverage reads in common and is therefore likely connected by a large deletion or insertion, or is a chimera.

In the case where it is impossible to know which mapping is the correct one, all read alignments will be classed as unchosen.

##### Assemble haplotypes

The last step is stitching blocks together to make haplotypes. If two blocks comes from the same haplotype but contains a large deletion or insertion, then there will be a left and right grouping of that haplotype.

... missing figure showing how blocks are stitched together ...

### FAQ

Since this algorithm is very new and is specialized to work on data available to me. I'm very willing to assist with any problems that may arise when the algorithm meets new data. I am also very open to any suggestions for improvements.


### The author (As of Febuary 2025) 

I, Laura Skak Rasmussen, am a recent graduate of bioinformatics from Aarhus University and now work as a PhD student at MOMA in Aarhus, Denmark working on exploring the epigenetics of individuals with sex chromosome aneuploidies.
