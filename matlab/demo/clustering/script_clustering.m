poolobj = parpool(12);
parfor task_id = 1:16
    params_grid_demo_clustering(task_id);
end
delete(poolobj);

