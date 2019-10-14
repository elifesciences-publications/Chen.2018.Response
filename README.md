# Auxier.response
These scripts and files were used to re-analyize the Chen et al., 2018 manuscript titled "Single nucleus sequencing reveals evidence of inter-nucleus recombination in arbuscular mycorrhizal fungi", and is associated with the paper from Auxier and Bazzicalupo, "Comment on ’Single nucleus sequencing reveals evidence of inter-nucleus recombination in arbuscular mycorrhizal
fungi’". eLife, 2019. http://doi.org/10.7554/eLife.47301



This analysis is based on single-nucleus data from three isolates; A4, A5, and SL1. 

The files to make each Figure are in their own folder. In this main folder can be found the python script that does the bulk of the work, this calculates the read depth using the pysam module. This script also calculates the blast hits for each position using biopython.
