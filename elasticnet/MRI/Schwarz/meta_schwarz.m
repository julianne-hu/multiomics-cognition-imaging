clear 

% Import data 
cd /Users/arshiya/Desktop/'Elastic Nets'/'Data'/
mri=readtable('mri_v2_schwarz_final.csv');
demo=readtable('demographics_v2.csv');
mri=innerjoin(mri, demo);

%% Just for metabolics 
metabolics=readtable('mean_imputed_metabolomic_i0_v2.csv');

%% extract ids and replace with actual names
protein_ids = metabolics.Properties.VariableNames(2:252);

% extract the field_id number from "pXXXXX_i0"
extracted_ids = regexp(protein_ids, 'p(\d+)_i0', 'tokens');
extracted_ids = cellfun(@(c) str2double(c{1}), extracted_ids);

% extract the ids from data dictionary
nightingale = readtable('Nightingale_biomarker_groups.txt');
field_ids = nightingale.field_id;
titles = string(nightingale.title);

% For each extracted ID, find the matching title
new_names = strings(size(extracted_ids));
for i = 1:length(extracted_ids)
    idx = find(field_ids == extracted_ids(i));
    if ~isempty(idx)
        new_names(i) = titles(idx);
    else
        new_names(i) = protein_ids{i}; % fallback to original name
    end
end

%assign new names to the table
metabolics.Properties.VariableNames(2:252) = cellstr(new_names);

%metabolics(:,find(contains(metabolics.Properties.VariableNames, 'i1')))=[];
%metabolics(isnan(metabolics.p20280_i0),:)=[];
%fields 23474 : 23467

merged_data=innerjoin(metabolics,mri);

%% pausing for sanity check

clear cognitive metabolics demo;
age=merged_data.p21003_i2;
sex=merged_data.p31;
site=merged_data.p54_i2;

%% Defining Y
Y=merged_data.schwarz_ct;

%% Clinical 
Xlogist=[age, strcmp(sex, 'Male'), strcmp(site, 'Newcastle (imaging)'), strcmp(site, 'Reading (imaging)'), strcmp(site, 'Cheadle (imaging)'), ]; %  , madrs_cut.MDMIS_01, madrs_cut.AIS_01]; %, madrs_cut.MDMIS_01, madrs_cut.EXEC_01, madrs_cut.AIS_01
varnames=[{'age'}, {'sex'},{'newcastle'} ,{'reading'} ,{'cheadle'} ]; %    , {'MemD'}, {'Attn'}]; % {'MemD'}, {'Exec'}, {'Attn'}
Ylogist=Y;

% removing NaNs
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

% Plotting figurre
% figure(13); scatter(yhat, Ylogist, '+','k'), lsline

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

    %predicted vs observed
    yhat = [ones(length(XTest),1),XTest]*coef;
    yhat_all=[yhat_all; yhat]; ytest_all=[ytest_all; yTest];
    r_x_split_clin(i)=corr(yhat, yTest);
end

r_xval_clin=corr(yhat_all, ytest_all)

%% rerun with metabolomics + clinical
Ylogist=Y;
X=merged_data{:, 2:252};
protein_ids=merged_data.Properties.VariableNames(2:252);

Xlogist=[X, age, strcmp(sex, 'Male'), strcmp(site, 'Newcastle (imaging)'), strcmp(site, 'Reading (imaging)'), strcmp(site, 'Cheadle (imaging)')];
varnames=[protein_ids, {'age'}, {'sex'},{'newcastle'} ,{'reading'} ,{'cheadle'}];

Xlogist(isnan(Ylogist),:)=[];
Ylogist(isnan(Ylogist))=[];
ix=sum(isnan(Xlogist)')';
Xlogist(ix>0,:)=[];
Ylogist(ix>0)=[];

[B,FitInfo] = lassoglm(Xlogist,Ylogist,'normal','CV',10, 'PredictorNames', varnames,'alpha', 0.1,"UseCovariance",true); %, 'RelTol',0.01
    %lassoPlot(B,FitInfo,'PlotType','CV'); legend('show','Location','best') % show legend
idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
MinModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinDeviance)~=0)
idxLambda1SE = FitInfo.Index1SE;
sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0)
B0 = FitInfo.Intercept(idxLambdaMinDeviance);
coef = [B0; B(:,idxLambdaMinDeviance)];
yhat = [ones(length(Xlogist),1),Xlogist]*coef; %y_pred = [ones(size(X_test ,1),1) X_test]*BETA;
r_yhat=corr(Ylogist,yhat)
%figure(14); scatter(yhat, Ylogist, '+','k'), lsline


%Xlogist=Xlogist(:,B(:,idxLambdaMinDeviance)~=0);varnames=varnames(B(:,idxLambdaMinDeviance)~=0); %parsimonious model test > set alpha to very low 0.001

%% cross-validation with metabolomics + clinical
permutation_index = randperm(length(Ylogist));%for randfold=1%:20
coef_all=[];
confusionmatrices=[];
yhat_all_meta=[];
ytest_all_meta=[];
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

    [B,FitInfo] = lassoglm(XTrain,yTrain,'normal','CV',10, 'PredictorNames', varnames,'alpha', 0.01, "UseCovariance",true);%0.0001 - ideal thresh for best model assessment , 'RelTol',0.01
    idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
    B0 = FitInfo.Intercept(idxLambdaMinDeviance);
    coef = [B0; B(:,idxLambdaMinDeviance)]; coef_all=[coef_all,coef];
    %predicted vs observed
    yhat = [ones(length(XTest),1),XTest]*coef;
    yhat_all_meta=[yhat_all_meta; yhat];
    ytest_all_meta=[ytest_all_meta; yTest];
    r_x_split_meta(i)=corr(yhat, yTest);
end

r_xval_meta=corr(yhat_all_meta, ytest_all_meta)
%figure; scatter(yhat_all, ytest_all, '+','k'), lsline

%% saving 

save meta_elnet_schwarz.mat
%% split

r_yhat.^2

r_xval.^2

r_x_split_clin

r_x_split_meta
figure(12);
for test=1:4;
subplot(2,2,test);
y=mean(r_x_split_clin(test,:));
err=std(r_x_split_clin(test,:));
errorbar(1, y, err,"square"); xlim([0 3])
y=mean(r_x_split_meta(test,:));
err=std(r_x_split_meta(test,:));
hold on; errorbar(1.5, y, err,"square"); xlim([0 3]); ylim([0.25 0.5])
end