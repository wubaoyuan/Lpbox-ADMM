clear variables;
prSet(3);

%% input
as = rand(1, 100000) * 10;
m = 10;

%% output
ht = tic;
[bs, idx] = sortTopSlow(as, m);
ti = toc(ht);

ht = tic;
[b2s, idx2] = sortTop(as, m);
ti2 = toc(ht);

equal('bs', bs, b2s);
equal('idx', idx, idx2);

fprintf('time : %f\n', ti);
fprintf('time2: %f\n', ti2);