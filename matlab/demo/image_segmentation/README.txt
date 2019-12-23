The codes in this folder will reproduce the experimental results reported in Table 1 of the main manuscript. 

There are four code files:
1) compare_all_methods_table1.m : this demo implements all compared methods in Table 1, except for Bethe-ADMM, which is implemented by C language by its authors. 
2) script_table1.m : this script calls 'params_grid_LP_box_ADMM_table1.m' to implement lpbox ADMM, with different p values, 
        different problem size (i.e., n), and different parameters. 
3) params_grid_LP_box_ADMM_table1.m : called by 'script_table1.m'
4) demo_summarize_result_of_lpbox_table1.m : summarize all results of lpbox ADMM produced by 'script_table1.m'. 
        For each p and n, there are many results with different parameters. Among these results, that with the minimal energy is selected as the reported result. 

Script:
>> run('compare_all_methods_table1')
>> run('script_table1')
>> run('demo_summarize_result_of_lpbox_table1')

-- Baoyuan Wu (wubaoyuan1987@gmail.com),  2018/06/08
