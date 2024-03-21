# consecutive_kmers
This is a project that simply detects consecutive k-mer repeats in sequences (sequencing reads).
For example, for sequence `ACTGATGATGGTGATGATGCGAA` and reading frame 0 (although all reading frames are considered) `ACT(GAT)_2 GGT(GAT)_2 GCGAA`
will be the output.

It contains:
1. `brainstorming` directory where new ideas are explored
2. the `cpp_code` directory that contains the c++ implementation
3. a `jupyter notebook` in the directory `eval` to analyze the output of the c++ program
### cpp_code
compile with
```bash
    g++ -std=gnu++17 main.cpp argparser.cpp -o main
```



