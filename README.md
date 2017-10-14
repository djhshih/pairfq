# pairfq

[![travis-ci](https://travis-ci.org/djhshih/pairfq.svg?branch=master)](https://travis-ci.org/djhshih/pairfq)
[![codecov](https://codecov.io/gh/djhshih/pairfq/branch/master/graph/badge.svg)](https://codecov.io/gh/djhshih/pairfq)

Pair first-read and second-second fastq files together, separating unpaired reads.

```{bash}
usage: pairfq <r1.fq> <r2.fq> <out.fq> [unpaired.fq]
```

Assume that each sequence and quality score entry `r1.fq` and `r2.fq` is single-line.
Assume that `r1.fq` and `r2.fq` are both sorted naturally.
(e.g. the source BAM file was sorted by `sambamba sort -N` or `samtools sort -n`).

The output file `out.fq` will contain interleaved read pairs with 
unpaired reads removed, and this file should be suitable for alignment with `bwa mem -p`.

Unpaired reads will be written to `unpaired.fq`, if available.

## Example workflow

Suppose we have reads aligned to hg19 and we want to re-align the the reads
to hg38. Starting with the input file `sample.bam`, we should ensure that
the reads are *not* sorted by coordinates, which could affect downstream
alignment due to biased insert size estimation.

We choose here to sort the reads by read name (in natural order), and
extract the sorted reads.

```{bash}
samtools sort -n -o qsorted.bam sample.bam
samtools fastq -0 marked-unpaired.fq -1 r1.fq -2 r2.fq qsorted.bam
```

If all unpaired reads in `sample.bam` are properly marked, the above command
should work as unexpected, producing `r1.fq` and `r2.fq` that
are already paired (and `pairfq` is unnecessary).
Otherwise, the below commands are necessary to ensure that `r1.fq` and `r2.fq`
are paired correctly.

```{bash}
pairfq r1.fq r2.fq out.fq unmarked-unpaired.fq
```

Now, we can re-align the read pairs to another reference genome by

```{bash}
bwa mem -p hg38.fa out.fq | samtools view -b - > sample_hg38.bam
```

