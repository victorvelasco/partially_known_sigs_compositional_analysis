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

## R Session information

```
> sessionInfo()
R version 4.3.1 (2023-06-16)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 10 (buster)

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.3.5.so;  LAPACK version 3.8.0

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/London
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.3.1
> packageVersion("rstan")
[1] ‘2.21.8’
```
