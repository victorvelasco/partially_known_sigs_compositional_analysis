# partially_known_sigs_compositional_analysis

Code supplement to Chapter 4 of my PhD thesis

## Analysis of simulated colorectal cancer data (Sect. 4.3)

To reproduce the analysis, run the following command
```R
Rscript script/analysis_colorectal_simulated_single.R patientId
```
Replacing patientId with the corresponding patient id (1..10).

The output is an object of class stanfit and will be saved in results/simulated/ in rds format.

To plot results, run the script plot_colorectal_simulated_single.R in RStudio

## Analysis of real oesophageal cancer data (Sect. 4.4)

To reproduce the analysis, run the following command
```R
Rscript script/analysis_oesophageal_real_single.R patientId
```
Replacing patientId with the corresponding patient id (1..10).

The output is an object of class stanfit and will be saved in results/read/ in rds format.

To plot results, run the script plot_oesophageal_real_single.R in RStudio

