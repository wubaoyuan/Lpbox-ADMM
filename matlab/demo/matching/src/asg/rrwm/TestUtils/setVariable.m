%% List all control vectors here: assignm control variable as X-axis
nVar = 0;

if ~isempty(settings.ratioFill)
    nVar = nVar + 1;
    variable{nVar} = settings.ratioFill;
    varName{nVar} = 'ratioFill';
    varStr{nVar} = '\fontname{times new roman}Edge density \rho';
    saveStr{nVar} = 'Density';
end

if ~isempty(settings.deformation)
    nVar = nVar + 1;
    variable{nVar} = settings.deformation;
    varName{nVar} = 'deformation';
    varStr{nVar} = '\fontname{times new roman}Deformation \sigma';
    saveStr{nVar} = 'Deform';
end

if ~isempty(settings.nOutlier)
    nVar = nVar + 1;
    variable{nVar} = settings.nOutlier;
    varName{nVar} = 'nOutlier';
    varStr{nVar} = '\fontname{times new roman}# of outliers \it{n_{out}}';
    saveStr{nVar} = 'Outlier';
end

%% Figure the variables out
nVector = 0;
for i = 1 : nVar
    if length(variable{i}) ~= 1
        nVector = nVector + 1;
    end
end
if nVector ~= 1
    error('There is not one control variable!! Fixed it!!')
end

nFixVar = 1;
for i = 1 : nVar
    if length(variable{i}) ~= 1
        Con.var = variable{i};
        Con.name = varName{i};
        Con.str = varStr{i};
        Con.saveStr = saveStr{i};
    else
        Fix(nFixVar).var = variable{i};
        Fix(nFixVar).name = varName{i};
        Fix(nFixVar).str = varStr{i};
        Fix(nFixVar).saveStr = saveStr{i};
        nFixVar = nFixVar + 1;
    end
end

settings.Con = Con;
settings.Fix = Fix;

clear i nFixVar nVar nVector Con Fix saveStr varName varStr variable
