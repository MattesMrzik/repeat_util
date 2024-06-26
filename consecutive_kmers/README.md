# consecutive_kmers
This is a project that simply detects consecutive k-mer repeats in sequences (sequencing reads).
For example, for sequence `ACTGATGATGGTGATGATGCGAA` and reading frame 0 (although all reading frames are considered) `ACT(GAT)_2 GGT(GAT)_2 GCGAA`
will be the output.

It contains:
1. `brainstorming` directory where new ideas are explored
2. the `cpp_code` directory that contains the c++ implementation ![CI Status](https://github.com/MattesMrzik/repeat_util/actions/workflows/c-cpp.yml/badge.svg)
3. a `jupyter notebook` in the directory `eval` to analyze the output of the c++ program

## cpp_code
can scan fasta, fastq, gzip fasta and fastq files (gzip -k file), and bam
### prerequisites:
- install [HSTlib](https://github.com/samtools/htslib/tree/develop) and make sure the library path is accessible to the program
- libseqan2-dev

### compile
```bash
make
```
or to set a new MAX_SEQ_LEN
```bash
make clean
make MAX_SEQ_LEN=200
```
### compile test with
```
make test
```
### Acknowledgments
- uses the header-only test framework [doctest](https://github.com/doctest/doctest)

## Open questions
- What about `(CAC)_2 CAGCAACAGCAA(CAG)_15 `? Do we also want to recognize the `(CAGCAA)_2` in the middle?
- Use different a score. Perhaps `CAG`-score together with `acc`. Maybe also a score for purity?
- Do we want to consider interruptions that are not multiplicative of `k`?


## TODO
make input -d option for dir to scan for fasta, fastq, and bam files
make output default dir to root dir ie the dir where it was executed from
