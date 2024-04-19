% Formated data in vector
inputSpecular = zeros([40, 2]);     % 1st column is correlation length, 2nd is rms
outputSpecular = zeros([40, 1]);
output = table2array(Results.Stat.meanROI);

for i = 1:length(corrAll)
    idx = (i-1)*10;
    for j = 1:length(rmsAll)
        InitParam = Results.InitParam;
        inputSpecular(idx + j, 1) = InitParam.Corr{j,i}{1};
        inputSpecular(idx + j, 2) = InitParam.Rms{j,i}{1};
        
        outputSpecular(idx + j, 1) = output(j,i);
    end
end

%% Linear Regression with both parameters
mdl = fitlm(inputSpecular,outputSpecular);

%% Correlation coefficient
coeff = corrcoef([outputSpecular inputSpecular]);

% With the ratio of correlation over rms
ratio = inputSpecular(:,2)./ inputSpecular(:,1);
coeffRatio = corrcoef([outputSpecular ratio]);
