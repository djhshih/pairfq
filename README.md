# pairfq

Pair first-read and second-second fastq files together, separating unpaired reads.

```{bash}
usage: pairfq <r1.fq> <r2.fq> <out.fq> [unpaired.fq]
```

Assume that each sequence and quality score entry `r1.fq` and `r2.fq` is single-line.
Assume that `r1.fq` and `r2.fq` are both sorted naturally.
(e.g. the source BAM file was sorted by `sambamba sort -N` or `samtools sort -n`).

Unpaired reads will be written to `unpaired.fq`, if available.

