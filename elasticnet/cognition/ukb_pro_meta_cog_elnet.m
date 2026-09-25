clear
cd /Users/juliannehu/Documents/UKB_proteomics/data/raw
cognitive=readtable('cognitive_v2.csv');
demo=readtable('demographics_v2.csv');
cognitive=innerjoin(cognitive, demo);

cognitive.p6350_i2(cognitive.p6350_i2>(nanmean(cognitive.p6350_i2)+4*nanstd(cognitive.p6350_i2)) | cognitive.p6350_i2<100)=NaN; %Duration to complete alphanumeric path
%cognitive.p6348_i2(cognitive.p6348_i2<100 | cognitive.p6348_i2> (nanmean(cognitive.p6348_i2)+3*nanstd(cognitive.p6348_i2)) )=NaN; %Duration to complete numeric/easy path
cognitive.tmt_cor=(cognitive.p6350_i2+5*cognitive.p6351_i2);% - (cognitive.x6348_2_0 +5*cognitive.x6349_2_0) ;

%proteomics
cd /Users/juliannehu/Documents/UKB_proteomics/data/cleaned
proteomics=readtable('mean_imputed_olink_i0.csv');
merged_data=innerjoin(proteomics,cognitive);
%metabolomics
cd /Users/juliannehu/Documents/UKB_proteomics/data/cleaned
metabolics=readtable('mean_imputed_metabolomic_i0_v2.csv');

%% extract ids and replace with actual names
protein_ids = metabolics.Properties.VariableNames(2:252);
% extract the field_id number from "pXXXXX_i0"
extracted_ids = regexp(protein_ids, 'p(\d+)_i0', 'tokens');
extracted_ids = cellfun(@(c) str2double(c{1}), extracted_ids);
% extract the ids from data dictionary
cd /Users/juliannehu/Documents/UKB_proteomics/data/dictionaries
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
merged_data=innerjoin(merged_data,metabolics);
clear cognitive metabolics demo;
age=merged_data.p21003_i2;
sex=merged_data.p31;
site=merged_data.p54_i2;

% define matrix to store all minmodel and sparsemodel features for each test
MinModelPredictors_all = cell(4,2);
sparseModelPredictors_all = cell(4,2);

test=1;
for test=1:4
    if test==1;
        Y=merged_data.tmt_cor;
    elseif test==2;
        Y=merged_data.p20016_i2;
    elseif test==3;
        Y=merged_data.p20197_i2;
    elseif test==4;
        Y=merged_data.p23324_i2;
    end
    disp("Running test " + test)

    %% run with clinical
    Xlogist=[age, strcmp(sex, 'Male'), strcmp(site, 'Newcastle (imaging)'), strcmp(site, 'Reading (imaging)'), strcmp(site, 'Cheadle (imaging)'), ]; %  , madrs_cut.MDMIS_01, madrs_cut.AIS_01]; %, madrs_cut.MDMIS_01, madrs_cut.EXEC_01, madrs_cut.AIS_01
    varnames=[{'age'}, {'sex'},{'newcastle'} ,{'reading'} ,{'cheadle'} ]; %    , {'MemD'}, {'Attn'}]; % {'MemD'}, {'Exec'}, {'Attn'}
    Ylogist=Y;
    %violinplot(Ylogist)

    Xlogist(isnan(Ylogist),:)=[];
    Ylogist(isnan(Ylogist))=[];
    ix=sum(isnan(Xlogist)')';
    Xlogist(ix>0,:)=[];
    Ylogist(ix>0)=[];

    [B,FitInfo] = lassoglm(Xlogist,Ylogist,'normal','CV',10, 'PredictorNames', varnames,'alpha', 0.1);
    %lassoPlot(B,FitInfo,'PlotType','CV'); legend('show','Location','best') % show legend
    idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
    MinModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinDeviance)~=0)
    idxLambda1SE = FitInfo.Index1SE;
    sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0)
    B0 = FitInfo.Intercept(idxLambdaMinDeviance);
    coef = [B0; B(:,idxLambdaMinDeviance)];
    yhat = [ones(length(Xlogist),1),Xlogist]*coef; %y_pred = [ones(size(X_test ,1),1) X_test]*BETA;
    r_yhat(test, 1)=corr(Ylogist,yhat)
    figure(11); subplot(4,1,test); scatter(yhat, Ylogist, '+','k'), lsline
    % store minmodel and sparsemodel features for each test
    MinModelPredictors_all{test,1} = MinModelPredictors;
    sparseModelPredictors_all{test,1} = sparseModelPredictors;

    %% cross-validation clinical
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
        %permutation_index = randperm(length(Ylogist));ix=zeros([length(Ylogist), 1]);
        %ix(permutation_index(1:50))=1;  %ix(permutation_index(31:304))=0;
        XTest=Xlogist(ix==1,:);
        XTrain=Xlogist(ix==0,:);
        yTest=Ylogist(ix==1);
        yTrain=Ylogist(ix==0);

        [B,FitInfo] = lassoglm(XTrain,yTrain,'normal','CV',10, 'PredictorNames', varnames,'alpha', 0.1);%0.0001 - ideal thresh for best model assessment
        idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
        B0 = FitInfo.Intercept(idxLambdaMinDeviance);
        coef = [B0; B(:,idxLambdaMinDeviance)]; coef_all=[coef_all,coef];
        %predicted vs observed
        yhat = [ones(length(XTest),1),XTest]*coef;
        yhat_all=[yhat_all; yhat]; ytest_all=[ytest_all; yTest];
        r_x_split_clin(test,i)=corr(yhat, yTest);
    end

    r_xval(test,1)=corr(yhat_all, ytest_all)
    figure(1); subplot(4,1,test); scatter(yhat_all, ytest_all, '+','k'), lsline


    %% rerun with proteomics + metabolomics + clinical
    Ylogist=Y;
    X=merged_data{:,[2:1460, 1565:1815]};
    protein_ids=merged_data.Properties.VariableNames([2:1460, 1565:1815]);
    %m=nanmean(X);
    %for i=1:length(protein_ids)
    %    X(isnan(X(:,i)),i)=m(i);
    %end
    %r=corr(X, Ylogist,'rows','pairwise'); thr=sortrows(abs(r)); thr=thr(2500); ix=abs(r)>thr;
    %X=X(:,ix);
    %protein_ids=protein_ids(ix);

    %Xm=merged_data{:,3025:3275};
    %metabolome_ids=merged_data.Properties.VariableNames(3025:3275);
    %r=corr(Xm, Ylogist,'rows','pairwise'); thr=sortrows(abs(r)); thr=thr(150); ix=abs(r)>thr;
    %Xm=Xm(:,ix);
    %metabolome_ids=metabolome_ids(ix);

    %X=[X,Xm];
    %protein_ids=[protein_ids, metabolome_ids];
    Xlogist=[X, age, strcmp(sex, 'Male'), strcmp(site, 'Newcastle (imaging)'), strcmp(site, 'Reading (imaging)'), strcmp(site, 'Cheadle (imaging)')];
    varnames=[protein_ids, {'age'}, {'sex'},{'newcastle'} ,{'reading'} ,{'cheadle'}];
    %violinplot(Ylogist)

    Xlogist(isnan(Ylogist),:)=[];
    Ylogist(isnan(Ylogist))=[];
    ix=sum(isnan(Xlogist)')';
    Xlogist(ix>0,:)=[];
    Ylogist(ix>0)=[];

    %idx = fscmrmr(Xlogist,Ylogist);
    %Xlogist=Xlogist(:,idx>200); varnames=varnames(idx>200);

    [B,FitInfo] = lassoglm(Xlogist,Ylogist,'normal','CV',10, 'PredictorNames', varnames,'alpha', 0.5,"UseCovariance",true, 'RelTol',0.01);
    %lassoPlot(B,FitInfo,'PlotType','CV'); legend('show','Location','best') % show legend
    idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
    MinModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinDeviance)~=0)
    idxLambda1SE = FitInfo.Index1SE;
    sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0)
    B0 = FitInfo.Intercept(idxLambdaMinDeviance);
    coef = [B0; B(:,idxLambdaMinDeviance)];
    yhat = [ones(length(Xlogist),1),Xlogist]*coef; %y_pred = [ones(size(X_test ,1),1) X_test]*BETA;
    r_yhat(test,2)=corr(Ylogist,yhat)
    figure(12); subplot(4,1,test);scatter(yhat, Ylogist, 'Marker', '+', 'MarkerEdgeColor', [0 0.5 0.5]); lsline

    % store minmodel and sparsemodel features for each test
    MinModelPredictors_all{test,2} = MinModelPredictors;
    sparseModelPredictors_all{test,2} = sparseModelPredictors;

    %Xlogist=Xlogist(:,B(:,idxLambdaMinDeviance)~=0);varnames=varnames(B(:,idxLambdaMinDeviance)~=0); %parsimonious model test > set alpha to very low 0.001

    %% cross-validation with proteomics + metabolomics + clinical
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
        %permutation_index = randperm(length(Ylogist));ix=zeros([length(Ylogist), 1]);
        %ix(permutation_index(1:50))=1;  %ix(permutation_index(31:304))=0;
        XTest=Xlogist(ix==1,:);
        XTrain=Xlogist(ix==0,:);
        yTest=Ylogist(ix==1);
        yTrain=Ylogist(ix==0);

        [B,FitInfo] = lassoglm(XTrain,yTrain,'normal','CV',10, 'PredictorNames', varnames,'alpha', 0.5, "UseCovariance",true, 'RelTol',0.01);        
        idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
        B0 = FitInfo.Intercept(idxLambdaMinDeviance);
        coef = [B0; B(:,idxLambdaMinDeviance)]; coef_all=[coef_all,coef];
        %predicted vs observed
        yhat = [ones(length(XTest(:,1)),1),XTest]*coef;
        yhat_all=[yhat_all; yhat];
        ytest_all=[ytest_all; yTest];
        r_x_split_meta(test,i)=corr(yhat, yTest);
    end

    figure(2); subplot(4,1,test);scatter(yhat_all, ytest_all, 'Marker', '+', 'MarkerEdgeColor', [0 0.5 0.5]); lsline
    r_xval(test,2)=corr(yhat_all, ytest_all)

end

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


figure; scatter(yhat_all, ytest_all, '+','MarkerEdgeColor',[0 0.5 0.5]), lsline
corr(yhat_all, ytest_all)

