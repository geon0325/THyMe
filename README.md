# THyMe+: Temporal Hypergraph Motifs and Fast Algorithms for Exact Counting
Source code for the paper [THyMe+: Temporal Hypergraph Motifs and Fast Algorithms for Exact Counting](https://github.com/geonlee0325/THyMe), Geon Lee and Kijung Shin, [ICDM 2021](https://icdm2021.auckland.ac.nz/).

*Group interactions arise in our daily lives (email communications,  on-demand ride sharing, comment interactions on online communities, to name a few), and they together form hypergraphs that evolve over time. Given such temporal hypergraphs, how can we describe their underlying design principles? If their sizes and time spans are considerably different, how can we compare their structural and temporal characteristics?*
*In this work, we define 96 *temporal hypergraph motifs* (TH-motifs), and propose the relative occurrences of their instances as an answer to the above questions. TH-motifs categorize the relational and temporal dynamics among three connected hyperedges that appear within a short time. For scalable analysis, we develop THyMe+, a fast and exact algorithm for counting the instances of TH-motifs in massive hypergraphs, and show that THyMe+ is at most *2,163X* *faster* while requiring less space than baseline. Using it, we investigate 11 real-world temporal hypergraphs from various domains. We demonstrate that TH-motifs provide important information useful for downstream tasks and reveal interesting patterns, including the striking similarity between temporal hypergraphs from the same domain.*

## Datasets
* The original datasets are available [here](https://www.cs.cornell.edu/~arb/data/).
* The processed datasets are available in [here](https://github.com/geonlee0325/THyMe/data).

## Running Demo
* There are 10 different **runtype**s which can be set in [run.sh](https://github.com/geonlee0325/covid_segmentation/blob/main/code/run.sh):
```setup
1.  Fitting using LLD/NLLD (our segmentation scheme)
2.  Fitting using LLD/NLLD (single segmentation)
3.  Fitting using LLD/NLLD (incremental segmentation)
4.  Forecasting using LLD/NLLD (our segmentation scheme)
5.  Forecasting using LLD/NLLD (single segmentation)
6.  Fitting with SIR (our segmentation scheme)
7.  Fitting with SIR (single segmentation)
8.  Fitting with SIR (incremental segmentation)
9.  Forecasting with SIR (our segmentation scheme)
10. Forecasting with SIR (single segmentation)
```
* For **runtype**s 1, 2, 3, 4, and 5, execute:
```setup
./run.sh [runtype] [country] [output directory] [LLD or NLLD] [latent dimension k] [error rate (only for runtype 3)]
e.g., ./run.sh 1 japan ./ NLLD 3
```
* For **runtype**s 6, 7, 8, 9, and 10, execute:
```setup
./run.sh [runtype] [country] [output directory] [SIR] [error rate (only for runtype 8)]
e.g., ./run.sh 6 japan ./ SIR
```

## Contact Information
If you have any questions, please contact [Geon Lee](https://geonlee0325.github.io/).
