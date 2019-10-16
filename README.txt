README

Source code for manuscript "Dynamic inference in general nested case-control designs" by J. Feifel and 
D. Dobler.

For questions, comments, or remarks about the code please contact J. Feifel (jan.feifel@uni-ulm.de).

The code has been written using R version 3.6.1 (Platform: x86_64-suse-linux-gnu, 64-bit) with package 
versions data.table_1.12.2 and survival_survival_2.44-1.1. 
The wild bootstrap via relation (8) and (9) is implemented in the function_breslow.R. 

In addition, a function deriving the confidence bands as well as the linear and log-transformed confidence 
intervals as discussed in Section 4 is provided (functions_conf_regions.R). The conf_regions function requires 
the list produced by the breslow function (see example). The input data.frame needs to be of the structure as 
displayed in Example_Dynamic_inference_general_NCC.R.

Additional, to the programm code, three data csv-files are available. data_SIR3_full.csv gives the full cohort 
as available within the kmi-package. data_SIR3_full_large.csv contains the same information as
data_SIR3_full.csv although it is inflate to mimic an nested-case control data set with stratification variable 
set. data_SIR3_NCC.csv gives a nested-case control data set with stratification variable set. 

To reproduce the functions in Figure 2 and Figure 3 presented in the main document, just run the main analysis 
file Example_Dynamic_inference_general_NCC.R. Please note, the correct directory has to be pre-specified! 

Note: The functions are restricted to the Cox model with 2 covariates, but extentions should be straightforward.

Information on the publicy available data set icu.pneu can be found in the manuscript, the Web Appendix as well 
as on the help page of the R package kmi (https://cran.r-project.org/web/packages/kmi/index.html).
