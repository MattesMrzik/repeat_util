## .gz files
created using GNU Gzip
```bash
gzip -k <file>
```
## .bam files
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562G1AlnRep1.bam
```
```
samtools view  -h  wgEncodeUwRepliSeqK562G1AlnRep1.bam | head -n 30 > small.sam
```

Changing some positions and reference names to check if reading them with the cpp code actually works. Also changing
some sequences to contain repeats

```
samtools view -bS small.sam -o small.bam
```