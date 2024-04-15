# repeat_util

This is some code I developed to detect and classify short tandem repeats in read data. This is part of my contribution to the [Yuen Lab](https://lab.research.sickkids.ca/yuen/team/) at the [Sick Kids Hospital](https://www.sickkids.ca/) in
Toronto as a volunteer under the supervision of [Induja Chandrakumar](https://scholar.google.com/citations?user=9WoRnf8AAAAJ).
![example workflow](https://github.com/github/docs/actions/workflows/c-cpp.yml/badge.svg)
This repo contains:
 1. `repeat_util.py` which is a Python script containing various utility methods. Not related to any particular project.
 2. The project `consecutive_kmers` that converts sequences like `ACTGATGATGGTGATGATGCGAA` and to `ACT(GAT)_2 GGT(GAT)_2
    GCGAA` and groups them.
