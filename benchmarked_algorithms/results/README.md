# Results of symmetry-detection methods benchmarked against MemSTATS
---

The following files indicate how each tested symmetry-detection algorithm was benchmarked against the **specific symmetries** 
in the MemSTATS dataset. The name of each file corresponds to 
the particular subset of symmetries and algorithm options considered. For example, 
*alpha* = alpha-helical transmembrane proteins; 
*beta* = beta-barrel transmembrane proteins;
*internal* = internal symmetry; 
*quat* = quaternary symmetry;
*inferred* = indicates that the EncoMPASS method is taken to include symmetries inferred from homologous chains; 
*zscore8* = the results of the SymD method were considered signficant if their z-TM-score was higher than 8.   
The date in the end of the file name indicates when the benchmarking was completed.

•       MemSTATS-20-May-19_alpha_internal-inferred_zscore10_2019-08-27_results.tsv

•       MemSTATS-20-May-19_alpha_internal-inferred_zscore8_2019-08-27_results.tsv

•       MemSTATS-20-May-19_alpha_internal_zscore10_2019-08-27_results.tsv

•       MemSTATS-20-May-19_alpha_internal_zscore8_2019-08-27_results.tsv

•       MemSTATS-20-May-19_alpha_quat_zscore10_2019-08-27_results.tsv

•       MemSTATS-20-May-19_alpha_quat_zscore8_2019-08-27_results.tsv

•       MemSTATS-20-May-19_beta_internal-inferred_zscore10_2019-08-27_results.tsv

•       MemSTATS-20-May-19_beta_internal-inferred_zscore8_2019-08-27_results.tsv

•       MemSTATS-20-May-19_beta_internal_zscore10_2019-08-27_results.tsv

•       MemSTATS-20-May-19_beta_internal_zscore8_2019-08-27_results.tsv

•       MemSTATS-20-May-19_beta_quat_zscore10_2019-08-27_results.tsv

•       MemSTATS-20-May-19_beta_quat_zscore8_2019-08-27_results.tsv


**Overall statistics** for the performance of the algorithms can be found in the following files:

•       MemSTATS-20-May-19_stats_zscore10_2019-08-27.txt

•       MemSTATS-20-May-19_stats_zscore8_2019-08-27.txt


The results are also stored as a dictionary in a **Python pickle** data stream:

•       benchmark_results_20-May-19_zscore10_2019-08-27.pkl

•       benchmark_results_20-May-19_zscore8_2019-08-27.pkl
