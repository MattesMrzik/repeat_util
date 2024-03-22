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

### Open questions
- What about `(CAC)_2 CAGCAACAGCAA(CAG)_15 `? Do we also want to recognize the `(CAGCAA)_2` in the middle?
- Use different a score. Perhaps `CAG`-score together with `acc`.
- Do we want to consider interruptions that are not multiplicative of `k`?


