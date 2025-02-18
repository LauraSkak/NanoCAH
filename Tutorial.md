# NanoCAH tutorial

In this page i will go through how you prepare your data for the NanoCAH algorithm and how you run the algorithm.

### Step 0: Create enviroment

```bash

conda create -n NanoCAH -c bioconda -c conda-forge pysam argparse sys samtools bwa

```

... missing ...

### Step 1: Nanopore long-read sequencing

It is recommended that you use data which have an average read depth of at least 15 to 30 reads. Samples with a read depth less than 15 reads in some cases can still yield good results, but it is less likely than at higher depth.

### Step 2: Data preparation

If the data received is in fast5 or pod5 format the data should first be basecalled using the [dorado](https://github.com/nanoporetech/dorado) basecaller. For this, refer to the part 1 in the [Modification count workflow](https://github.com/LauraSkak/Modification-count-workflow). Modification calling is not necessary in this case.

If the data received is in the modbam format, then the data can be converted to a fastq format.

... A instructive command is missing here ...

If the data received is in the fastq format, then the data is ready for the next step.

### Step 3: Alignment

The NanoCAH needs aligned data, so the next step is to align med reads to our reference.

Input data:

* reference - Path to genome reference used. hg38 is recommended since CAH variants are annotated with hg38 positions.
* threads - The thread count, which will depend on your setup. If possible 32 or 64 threads are recommended.
* infile - Path to fastq file with all basecalled reads.
* outfile - Desired path to alignment outfile.

```bash

minimap2 -ax map-ont -Y {reference} -t {threads} {infile} | samtools view -bS | samtools sort > {outfile}

```
Tested with minimap2 v. 2.28-r1209.

The data should remain unfiltered because the NanoCAH algorithm utilizes reads with multiple mapping, which would likely otherwise be removed.

### Step 4: The NanoCAH algorithm

* reference - The genome reference as the one used for alignment.
* range - The range which will undergo local rephasing. It is recommended to chose range that span over the entire region of segmental duplication with an additional 1M basepairs on each flank. If hg38 as a reference and the target is the CYP21A2 gene, the recommended range is chr6:31900000-32200000. 
* infile - The path to the alignment file produced in step 3.
* Extra flags are added if you want to change input parameters or output file paths from the default.

```bash

python NanoCAH.py \
        -i {infile} \
        -r {reference} \
        -s {range} \
        {extra flags}

```

The NanoCAH will output a phased and filtered alignment file and a fasta file containing the assembly for each haplotype.

### Step 5: Inspecting result in IGV (optional)

To inspect the resulting filtered alignment, the alignment file needs to be sorted and indexed.

```bash

samtools sort {alignment outfile} > {sorted alignment outfile}
samtools index {sorted alignment outfile}

```

To inspect the haplotype assemblies, the fasta file should be mapped, sorted and indexed.

```bash

bwa mem {reference} {assembly outfile} > {assembly alignment outfile}

samtools view -bS {assembly alignment outfile} | samtools sort > {sorted assembly alignment outfile}
samtools index {sorted assembly alignment outfile}

```

Tested with samtools v. 1.6 and bwa v. ... missing ... .