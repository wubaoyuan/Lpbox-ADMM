The codes in this folder will reproduce the experimental results of lpbox ADMM reported in Table 2 of the main manuscript.

There are four code files:
1) script_clustering.m : this script calls 'params_grid_demo_clustering.m' to implement lpbox ADMM, with different p values,
        different UCI datasets, and different parameters.
2) params_grid_LP_box_ADMM_table1.m : called by 'script_clustering.m'
3) result_summarization_table2.m : summarize all results of lpbox ADMM produced by 'script_clustering.m'.
        For each p and UCI data, there are many results with different parameters. Among these results, that with the minimal objective value is selected as the reported result.
4) data_UCI : including four UCI datasets, iris, wine, glass, letter

-- Baoyuan Wu (wubaoyuan1987@gmail.com), 2018/06/08
