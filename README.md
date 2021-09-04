# THyMe+: Temporal Hypergraph Motifs and Fast Algorithms for Exact Counting
Source code for the paper [THyMe+: Temporal Hypergraph Motifs and Fast Algorithms for Exact Counting](https://github.com/geonlee0325/THyMe), Geon Lee and Kijung Shin, [ICDM 2021](https://icdm2021.auckland.ac.nz/).

*Group interactions arise in our daily lives (email communications,  on-demand ride sharing, comment interactions on online communities, to name a few), and they together form hypergraphs that evolve over time. Given such temporal hypergraphs, how can we describe their underlying design principles? If their sizes and time spans are considerably different, how can we compare their structural and temporal characteristics?*
*In this work, we define 96 *temporal hypergraph motifs* (TH-motifs), and propose the relative occurrences of their instances as an answer to the above questions. TH-motifs categorize the relational and temporal dynamics among three connected hyperedges that appear within a short time. For scalable analysis, we develop THyMe+, a fast and exact algorithm for counting the instances of TH-motifs in massive hypergraphs, and show that THyMe+ is at most *2,163X* *faster* while requiring less space than baseline. Using it, we investigate 11 real-world temporal hypergraphs from various domains. We demonstrate that TH-motifs provide important information useful for downstream tasks and reveal interesting patterns, including the striking similarity between temporal hypergraphs from the same domain.*

## Datasets
* The original datasets are available [here](https://www.cs.cornell.edu/~arb/data/).
* The processed datasets (email-Enron & contact-primary) are available in [here](https://github.com/geonlee0325/THyMe/tree/main/data).

## How to Run
* To run demos, execute following commend in [code](https://github.com/geonlee0325/THyMe/tree/main/code):
```
./run_(dp/thyme/thymeP).sh
```
* To run DP, THyMe, and THyMe+ with specified dataset and time interval (delta), execute following commends in [code](https://github.com/geonlee0325/THyMe/tree/main/code):
```
g++ -O3 -std=c++11 main_(dp/thyme/thymeP).cpp -o run;
./run (dataset) (delta)
```
* Results will be saved in [results](https://github.com/geonlee0325/THyMe/tree/main/results) as:
```
Filename: (dataset)_(delta)_(dp/thyme/thymeP).txt

Runtime (sec.)
1   (# of instances of TH-motif 1)
2   (# of instances of TH-motif 2)
...
96  (# of instances of TH-motif 96)
```

## Contact Information
If you have any questions, please contact [Geon Lee](https://geonlee0325.github.io/).
