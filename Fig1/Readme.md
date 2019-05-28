This folder contains the files needed to reproduce Figure 1. For each of the three isolates, there is an R script that uses the starting data from Supplement 6 of Chen et al., called supp6.[isolate].raw.csv, and [isolate]_recomb_depths.csv. This last file is computed using the python script to calculate the number of reads for each nucleus at each position.

The R script takes the input .raw.csv and removes sites that do not pass the quality filters. The syntax is quite complicated, so do not hesistate to contact regarding any questions.
