# repeat_util

## consecutive_kmers.py
### input 
You must provide one of the following:
1. `--sequence` and pass a string on the command line
2. `--fasta` a path to fasta file
3. `--fastq` a path to fastq file

Optional: `-k` for a specific k-mer size to look for. Default uses 3,4, and 5.

### output
For a sequence all consecutive k-mer repetitions are found in all reading frames. 

For example, for sequence `ACTGATGATGGTGATGATGCGAA` and reading frame 0 `ACT(GAT)_2 GGT(GAT)_2 GCGAA` will be the output.
