%% Algorithms (functions) to be evaluated
nAlgorithm = 0;

%% Proposed Method
if 1
    nAlgorithm = nAlgorithm + 1;
    Algorithm(nAlgorithm).fhandle = @RWPRM_ECCV2010;
    Algorithm(nAlgorithm).strName = 'RRWM';
    Algorithm(nAlgorithm).postProcess = @postHungarian;
    Algorithm(nAlgorithm).color = 'r';
    Algorithm(nAlgorithm).lineStyle = '-';
    Algorithm(nAlgorithm).Marker = 'p';
end

%% Spectral Matching Method
if 1
    nAlgorithm = nAlgorithm + 1;
    Algorithm(nAlgorithm).fhandle = @SM;
    Algorithm(nAlgorithm).strName = 'SM';
    Algorithm(nAlgorithm).postProcess = @postHungarian;
    Algorithm(nAlgorithm).color = 'k';
    Algorithm(nAlgorithm).lineStyle = '-';
    Algorithm(nAlgorithm).Marker = 'x';
end

disp(['There are ' num2str(nAlgorithm) ' method'])
for k = 1 : nAlgorithm
    disp(Algorithm(k).strName)
end

settings.Algorithm = Algorithm;
clear Algorithm nAlgorithm k strName
