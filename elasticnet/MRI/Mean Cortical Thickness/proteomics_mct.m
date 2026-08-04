clear

% MRI DATA
cd /Users/arshiya/Desktop/'Elastic Nets'/'Data'/
mri=readtable('mri_v2_mean_ct_final.csv');
demo=readtable('demographics_v2.csv');
mri=innerjoin(mri,demo);
proteomics=readtable('mean_imputed_olink_i0.csv');
merged_data=innerjoin(proteomics,mri);

% Define variables
age=merged_data.p21003_i2;
sex=merged_data.p31;
site=merged_data.p54_i2;
clear mri proteomics demo;
Y=merged_data.mean_cortical_thickness;

%% Clinical Model
Xlogist=[age, strcmp(sex, 'Male'), strcmp(site, 'Newcastle (imaging)'), strcmp(site, 'Reading (imaging)'), strcmp(site, 'Cheadle (imaging)'), ]; %  , madrs_cut.MDMIS_01, madrs_cut.AIS_01]; %, madrs_cut.MDMIS_01, madrs_cut.EXEC_01, madrs_cut.AIS_01
varnames=[{'age'}, {'sex'},{'newcastle'} ,{'reading'} ,{'cheadle'} ]; %    , {'MemD'}, {'Attn'}]; % {'MemD'}, {'Exec'}, {'Attn'}
Ylogist=Y;

% Removing NaNs
Xlogist(isnan(Ylogist),:)=[];
Ylogist(isnan(Ylogist))=[];
ix=sum(isnan(Xlogist)')';
Xlogist(ix>0,:)=[];
Ylogist(ix>0)=[];

% Running elastic net
[B,FitInfo] = lassoglm(Xlogist,Ylogist,'normal','CV',10, 'PredictorNames', varnames,'alpha', 0.001);
idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
MinModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinDeviance)~=0)
idxLambda1SE = FitInfo.Index1SE;
sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0)
B0 = FitInfo.Intercept(idxLambdaMinDeviance);
coef = [B0; B(:,idxLambdaMinDeviance)];
yhat = [ones(length(Xlogist),1),Xlogist]*coef; %y_pred = [ones(size(X_test ,1),1) X_test]*BETA;
r_yhat=corr(Ylogist,yhat)

% Plotting figure
figure(13); scatter(yhat, Ylogist, '+','k'), lsline

%% Cross Validation
permutation_index = randperm(length(Ylogist));%for randfold=1%:20
coef_all=[];
confusionmatrices=[];
yhat_all=[];
ytest_all=[];
yhatBinom_all=[];
n=length(Ylogist);
train_n=round(length(Ylogist)/10) %
for i=1:10 %8
    i
    if(i==1);
        ix=zeros([n,1]);
        ix(1:train_n)=1;
    elseif (i==10)
        ix=zeros([n,1]);
        ix(train_n*(i-1)+1:n)=1;
    else
        ix=zeros([n,1]);
        ix(train_n*(i-1)+1:train_n*(i-1)+train_n)=1;
    end;
    ix=ix(permutation_index);
    XTest=Xlogist(ix==1,:);
    XTrain=Xlogist(ix==0,:);
    yTest=Ylogist(ix==1);
    yTrain=Ylogist(ix==0);

    [B,FitInfo] = lassoglm(XTrain,yTrain,'normal','CV',10, 'PredictorNames', varnames,'alpha', 0.001);%0.0001 - ideal thresh for best model assessment
    idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
    B0 = FitInfo.Intercept(idxLambdaMinDeviance);
    coef = [B0; B(:,idxLambdaMinDeviance)]; coef_all=[coef_all,coef];

    % Predicted vs observed
    yhat = [ones(length(XTest),1),XTest]*coef;
    yhat_all=[yhat_all; yhat]; ytest_all=[ytest_all; yTest];
    r_x_split(i)=corr(yhat, yTest);
end

r_xval_clin=corr(yhat_all, ytest_all)

%% Clinical and proteomics
Ylogist=Y;
X = merged_data{:, 2:1460}; % double check this one
protein_ids=merged_data.Properties.VariableNames(2:1460);

Xlogist=[X, age, strcmp(sex, 'Male'), strcmp(site, 'Newcastle (imaging)'), strcmp(site, 'Reading (imaging)'), strcmp(site, 'Cheadle (imaging)')];
varnames=[protein_ids, {'age'}, {'sex'},{'newcastle'} ,{'reading'} ,{'cheadle'}];
    
% Removing NaNs
Xlogist(isnan(Ylogist),:)=[];
Ylogist(isnan(Ylogist))=[];
ix=sum(isnan(Xlogist)')';
Xlogist(ix>0,:)=[];
Ylogist(ix>0)=[];

% Elastic net regression
[B,FitInfo] = lassoglm(Xlogist,Ylogist,'normal','CV',10, 'PredictorNames', varnames,'alpha', 0.1, "UseCovariance",true);%, 'MaxIter',50);, ,'RelTol',0.01
idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
MinModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinDeviance)~=0)
idxLambda1SE = FitInfo.Index1SE;
sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0)
B0 = FitInfo.Intercept(idxLambdaMinDeviance);
coef = [B0; B(:,idxLambdaMinDeviance)];
yhat = [ones(length(Xlogist),1),Xlogist]*coef;
r_yhat=corr(Ylogist,yhat)
%figure(15); scatter(yhat, Ylogist, '+','b'), lsline

%% Cross-validation with proteomics + clinical

permutation_index = randperm(length(Ylogist));%for randfold=1%:20
coef_all_prot=[];
confusionmatrices=[];
yhat_all_prot=[];
ytest_all_prot=[];
yhatBinom_all=[];
n=length(Ylogist);
train_n=round(length(Ylogist)/10) %
for i=1:10 %8
    i
    if(i==1);
        ix=zeros([n,1]);
        ix(1:train_n)=1;
    elseif (i==10)
        ix=zeros([n,1]);
        ix(train_n*(i-1)+1:n)=1;
    else
        ix=zeros([n,1]);
        ix(train_n*(i-1)+1:train_n*(i-1)+train_n)=1;
    end;
    ix=ix(permutation_index);
        %permutation_index = randperm(length(Ylogist));ix=zeros([length(Ylogist), 1]);
        %ix(permutation_index(1:50))=1;  %ix(permutation_index(31:304))=0;
    XTest=Xlogist(ix==1,:);
    XTrain=Xlogist(ix==0,:);
    yTest=Ylogist(ix==1);
    yTrain=Ylogist(ix==0);

    [B,FitInfo] = lassoglm(XTrain,yTrain,'normal','CV',10, 'PredictorNames', varnames,'alpha', 0.01, "UseCovariance",true);%0.0001 - ideal thresh for best model assessment  'RelTol',0.01
    idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
    B0 = FitInfo.Intercept(idxLambdaMinDeviance);
    coef = [B0; B(:,idxLambdaMinDeviance)]; coef_all_prot=[coef_all_prot,coef];
    %predicted vs observed
    yhat = [ones(length(XTest(:,1)),1),XTest]*coef;
    yhat_all_prot=[yhat_all_prot; yhat]; 
    ytest_all_prot=[ytest_all_prot; yTest];
    r_x_split_prot(i)=corr(yhat, yTest);
end

r_xval_prot=corr(yhat_all_prot, ytest_all_prot)

%% Saving results
save clinical_prot_elnet.mat
