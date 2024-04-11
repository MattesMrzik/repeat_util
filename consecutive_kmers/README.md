# consecutive_kmers
This is a project that simply detects consecutive k-mer repeats in sequences (sequencing reads).
For example, for sequence `ACTGATGATGGTGATGATGCGAA` and reading frame 0 (although all reading frames are considered) `ACT(GAT)_2 GGT(GAT)_2 GCGAA`
will be the output.

It contains:
1. `brainstorming` directory where new ideas are explored
2. the `cpp_code` directory that contains the c++ implementation
3. a `jupyter notebook` in the directory `eval` to analyze the output of the c++ program
## cpp_code
### prerequisites:
- install [HSTlib](https://github.com/samtools/htslib/tree/develop) and make sure the library path is accessible to the program
- libseqan2-dev

### compile main with
```bash
    g++ -std=gnu++17 main.cpp argparser.cpp -o main -lhts
    g++ -std=gnu++17 main.cpp argparser.cpp ConsecutiveKmers.cpp -o main -lhts -DSEQAN_HAS_ZLIB -lz
```
### compile test with
```bash
    g++ -std=gnu++17 test.cpp argparser.cpp ConsecutiveKmers.cpp -o test -lhts -DSEQAN_HAS_ZLIB -lz
```
### Acknowledgments
- uses the header-only test framework [doctest](https://github.com/doctest/doctest)

## Open questions
- What about `(CAC)_2 CAGCAACAGCAA(CAG)_15 `? Do we also want to recognize the `(CAGCAA)_2` in the middle?
- Use different a score. Perhaps `CAG`-score together with `acc`. Maybe also a score for purity?
- Do we want to consider interruptions that are not multiplicative of `k`?


## TODO
make input -d option for dir to scan for fasta, fastq, and bam files
