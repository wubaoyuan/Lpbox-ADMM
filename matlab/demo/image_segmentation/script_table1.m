
poolobj = parpool(6);
parfor task_id = 1:160
    params_grid_LP_box_ADMM_table1(task_id);
end
delete(poolobj);
