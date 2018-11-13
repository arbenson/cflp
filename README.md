# core-fringe link prediction (cflp)

This code accompanies the paper

- Core-fringe link prediction. Austin R. Benson and Jon Kleinberg. 2018.



### Datasets and reproduce summary statistic tables

Public data is included in the repository. Here is code to reproduce the summary statistics of the datasets that appear in Tables 1 and 2 of the paper.

```julia
include("paper_tables.jl")
summary_stats("email-Enron")   # --> 148 & 18444 & 1344 & 41883 (Table 1)
summary_stats("call-Reality")  # --> 91 & 8927 & 127 & 10512  (Table 2)
```



### Running prediction experiments

Here is the code snippet for how to re-run the prediction experiments for the call-Reality dataset.

```julia
include("predict.jl")

# Re-run experiments for call-Reality dataset.
dataset = "call-Reality"
for rand_order in [false, true]  # random order or not
    for trial = 1:10             # 10 random trials
    	make_predictions(dataset, trial, rand_order)
    end
end
```

The output files are in `output/call-Reality`. The file `output/call-Reality/scores-call-Reality-4.txt` corresponds to the fourth trial with an ordering by most connected. The results look something like the following:

```
$ head -5 output/call-Reality/scores-call-Reality-4.txt

0.6634615384615384 0.6057692307692307 690
0.6423076923076924 0.5865384615384616 2096
0.6903846153846154 0.6461538461538462 3428
0.6884615384615385 0.6596153846153846 4550
0.6923076923076923 0.6653846153846154 5200
```

These first five lines correspond to 0, 1, 2, 3, or 4 fringe nodes. The first column is the number of common neighbors prediction score, the second column is the Jaccard similarity prediction score, and the third column is the number of length-2 paths between end points of edges in the test set.

We have pre-computed all of the scores for all of the public datasets (email-Enron, email-W3C, call-Reality, and text-Reality).



### Reproducing figures

Here is the code snippet for reproducing the figures from empirical observations in Section 2.

```julia
include("paper_plots.jl")

# Figures 2A and 2E (output email-Enron-{Common_Neighbors,Jaccard}.pdf)
empirical_plot("email-Enron", "Common Neighbors")
empirical_plot("email-Enron", "Jaccard")

# Figures 3A and 3E (output email-Enron-1-{Common_Neighbors,Jaccard}.pdf)
empirical_plot("email-Enron-1", "Common Neighbors")
empirical_plot("email-Enron-1", "Jaccard")

# Figures 4A and 4C (output call-Reality-{Common_Neighbors,Jaccard}.pdf)
empirical_plot("call-Reality", "Common Neighbors")
empirical_plot("call-Reality", "Jaccard")
```



Here is the code snippet for reproducing the figures on random graph models in Section 3.

```julia
include("paper_plots.jl")

# Figure 7A and 7B (outputs SBM1.pdf and SBM2.pdf)
p, q, r = 0.5, 0.3, 0.2
num_core = 10
sbm_plot(p, q, r, [0.0, 0.2], num_core, 0:100, "SBM1.pdf")
sbm_plot(p, q, r, [0.1],      num_core, 0:300, "SBM2.pdf")

# Figure 8 (output lattice-9-10.pdf)
b, c = 9, 10
lattice_plot(2:6, b, c, 0:12)
```

