clear
% Importing files
cd /Users/juliannehu/Documents/UKB_proteomics/data/raw
cognitive=readtable('cognitive_v2.csv');
demo=readtable('demographics_v2.csv');
cd /Users/juliannehu/Documents/UKB_proteomics/data
demo_bmi=readtable('ukb_demo.csv');
demo_other_covariates=readtable('ukb_demo_ses_edu_cvd.csv');
demo_all_covariates=innerjoin(demo_bmi,demo_other_covariates);
demo=outerjoin(demo,demo_all_covariates);
demo.Properties.VariableNames{'eid_demo'} = 'eid';
cd /Users/juliannehu/Documents/UKB_proteomics/data/cleaned
proteomics=readtable('mean_imputed_olink_i0.csv');
cognitive=innerjoin(cognitive, demo);

% Calculating trail making test scores
cognitive.p6350_i2(cognitive.p6350_i2>(nanmean(cognitive.p6350_i2)+4*nanstd(cognitive.p6350_i2)) | cognitive.p6350_i2<100)=NaN; % remove outliers
cognitive.tmt_cor=(cognitive.p6350_i2+5*cognitive.p6351_i2); % add penalty for errors
%figure;histogram(Y(:,1))
% Joining files
merged_data=innerjoin(proteomics,cognitive);
clear cognitive proteomics demo demo_bmi demo_other_covariates demo_all_covariates;

% Variables
i0age=merged_data.p21003_i0; 
i2age=merged_data.p21003_i2;
sex=merged_data.p31; 
site=merged_data.p54_i2;
bmi=merged_data.x21001_0_0;
bmi(bmi < 0) = NaN;
education=merged_data.x6138_0_0;
education(education < 0) = NaN;
cvd_cols = merged_data(:, 1590:1743);
cvd = double(any(~ismissing(cvd_cols), 2));
ses=merged_data.x738_0_0;
ses(ses < 0) = NaN;

% X = proteomics
X = merged_data{:, 2:1460}; 
protein_ids=merged_data.Properties.VariableNames(2:1460);
% Y: alphanumerical trails, fluid intelligence,  word pairs (PAL), symbol digit substitution
Y=[merged_data.tmt_cor, merged_data.p20016_i2, merged_data.p20197_i2, merged_data.p23324_i2];%   Y=[cognitive.tmt_cor, cognitive.x21004_2_0, cognitive.x20197_2_0]; % cognitive.x399_2_2+cognitive.x399_2_1, 

%%% Removing NaNs
naninx=sum(isnan(Y)')' | sum(isnan(X)')' | isnan(bmi) | isnan(education) | isnan(ses);
Y=Y(naninx==0,:); 
X=X(naninx==0,:); 
i0age=i0age(naninx==0); 
i2age=i2age(naninx==0); 
sex=sex(naninx==0);
site=site(naninx==0); 
bmi=bmi(naninx==0);
education=education(naninx==0);
cvd=cvd(naninx==0);
ses=ses(naninx==0);

% Standardizing data
Y=zscore(Y); 
X=zscore(X); 

% Flip direction of trail making test (higher score = worse performance)
Y(:,1)=Y(:,1)*-1;

% Set # of PLS components
ncomp=3;

% defining other covariates
i0age_squared = i0age.^2;
i2age_squared = i2age.^2;
sex=categorical(sex);
sex_num = double(sex == 'Male');
i0age_by_sex = i0age .* sex_num;
i2age_by_sex = i2age .* sex_num;

% Regress out age, sex, and site (covariates)
mdl = fitlm(table(i2age, sex, site, i2age_squared, i2age_by_sex, bmi, education, cvd, ses, Y(:,1)) ); Y(:,1)=mdl.Residuals.Raw;
mdl = fitlm(table(i2age, sex, site, i2age_squared, i2age_by_sex, bmi, education, cvd, ses, Y(:,2)) ); Y(:,2)=mdl.Residuals.Raw;
mdl = fitlm(table(i2age, sex, site, i2age_squared, i2age_by_sex, bmi, education, cvd, ses, Y(:,3)) ); Y(:,3)=mdl.Residuals.Raw;
mdl = fitlm(table(i2age, sex, site, i2age_squared, i2age_by_sex, bmi, education, cvd, ses, Y(:,4)) ); Y(:,4)=mdl.Residuals.Raw;
clear x; 
for i=1:length(X(1,:)) 
    mdl = fitlm(table(i0age, sex, site, i0age_squared, i0age_by_sex, bmi, education, cvd, ses, X(:,i))); 
    x(:,i)=mdl.Residuals.Raw;
end

% Run PLS
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);
PCTVAR

% Plotting each of the components' correlation with cognitive tests
figure;
imagesc(corr(XS(:,1:3),Y)); 
colormap bone; 
colorbar; 
corr(XS, Y) %1.TMT  2.GF  3.PAL  4.DSST
combs=allcomb([1:3], [1:4]); 
figure; clf;
for i=1:length(combs)
    hold off; 
    ix1=combs(i,1);
    ix2=combs(i,2); 
    subplot(3,4,i); 
    scatter(XS(:,ix1), Y(:,ix2), 5, 'k', 'filled');
end; 
clear mdl ypred yci mdl ix1 ix2 combs

% Calculate correlations
corr(i0age, XS)
corr(XS, Y)
max(corr(X, Y))
min(corr(X, Y))

index=strcmp(protein_ids, 'il6');
corr(X(:, index), Y)
figure;dscatter(X(:, index), Y(:, 4));lsline

index=strcmp(protein_ids, 'vegfd');
corr(X(:, index), Y)
figure;dscatter(X(:, index), Y(:, 1));lsline

%% permutation testing
permutations=5000;   
allobservations=Y; 
for ncomp=1:3
    parfor n = 1:permutations;
    % selecting either next combination, or random permutation
    permutation_index = randperm(length(allobservations));
    % creating random sample based on permutation index
    randomSample = allobservations(permutation_index,:);
    % running the PLS for this permutation
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,randomSample,ncomp);
    Rsq(n) = sum(PCTVAR(2,:));
    Rsq1(n) = sum(PCTVAR(1,:));
    end
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);
    p(ncomp)=sum(sum(PCTVAR(2,:))<Rsq')/permutations
    p_1(ncomp)=sum(sum(PCTVAR(1,:))<Rsq1')/permutations
end
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);PCTVAR

% figure to show the permutation distribution vs actual distribution
figure; 
histogram(Rsq, 'FaceColor', [0.5 0.5 0.5]); 
xlim([0.01 0.035]); 
xline(0.0277, 'r', 'LineWidth', 2); 
% Only show x and y axes
ax = gca;
ax.Box = 'off'; % removes top and right lines
ax.TickDir = 'out'; % ticks point outward
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
sum(PCTVAR')
%% bootstrapping to get the func connectivity weights for PLS1, 2 and 3
dim=3
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,dim);PCTVAR
PLS1w=stats.W(:,1);
PLS2w=stats.W(:,2);
PLS3w=stats.W(:,3);

bootnum=5000;
PLS1weights=[];
PLS2weights=[];
PLS3weights=[];

parfor i=1:bootnum
    i;
    myresample = randsample(size(x,1),size(x,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=x(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data 
    newW=stats.W(:,1);%extract PLS1 weights
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run
    newW=stats.W(:,2);%extract PLS2 weights
    if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run    
    
    newW=stats.W(:,3);%extract PLS2 weights
    if corr(PLS3w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS3weights=[PLS3weights,newW]; %store (ordered) weights from this bootstrap run    
end

%calculate bootstrap standard error
PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');
PLS3sw=std(PLS3weights');

%calculate bootstrap ratio (a z-score)
plsweights1=PLS1w./PLS1sw';
plsweights2=PLS2w./PLS2sw'; 
plsweights3=PLS3w./PLS3sw';

%flip sign of PLS 2 and 3 weights to match upregulated vs downregulated
plsweights2 = plsweights2*-1;
stats.W(:,2) = -stats.W(:,2);
plsweights3 = plsweights3*-1;
stats.W(:,3) = -stats.W(:,3);

% filtering bootstrap results to identify variables that reliably contribute to each component
sum(plsweights1 > 3)
sum(plsweights1 < -3)

sum(plsweights2 > 3)
sum(plsweights2 < -3)

sum(plsweights3 > 3)
sum(plsweights3 < -3)

protein_ids(plsweights1 > 3)'
protein_ids(plsweights1 < -3)'

protein_ids(plsweights2 > 3)'
protein_ids(plsweights2 < -3)'

protein_ids(plsweights3 > 3)'
protein_ids(plsweights3 < -3)'



%% visuals
% %FIGURE WITH OLINK EXPLORE
 cd /Users/juliannehu/Documents/UKB_proteomics/data
 cardiometabolic = readtable('Olink_Cardiometabolic.csv');  
 cardiometabolic_names = string(cardiometabolic.Gene);
 cardiometabolic_2 = readtable('Olink_Cardiometabolic_II.csv');  
 cardiometabolic_2_names = string(cardiometabolic_2.Gene);
 neurology = readtable('Olink_Neurology.csv');  
 neurology_names = string(neurology.Gene);
 neurology_2 = readtable('Olink_Neurology_II.csv');  
 neurology_2_names = string(neurology_2.Gene);
 oncology = readtable('Olink_Oncology.csv');  
 oncology_names = string(oncology.Gene);
 oncology_2 = readtable('Olink_Oncology_II.csv');  %
 oncology_2_names = string(oncology_2.Gene);
 inflammation = readtable('Olink_Inflammation.csv');  
 inflammation_names = string(inflammation.Gene);
 inflammation_2 = readtable('Olink_Inflammation_II.csv');  
 inflammation_2_names = string(inflammation_2.Gene);

 protein_ids_str = string(protein_ids);

% identify cardiometabolic proteins
 is_cardiometabolic = ismember(lower(protein_ids_str), lower(cardiometabolic_names));
% identify cardiometabolic II proteins
 is_cardiometabolic_2 = ismember(lower(protein_ids_str), lower(cardiometabolic_2_names));
% identify neurology proteins
 is_neurology = ismember(lower(protein_ids_str), lower(neurology_names));
% identify neurology II proteins
 is_neurology_2 = ismember(lower(protein_ids_str), lower(neurology_2_names));
% identify oncology proteins
 is_oncology = ismember(lower(protein_ids_str), lower(oncology_names));
% identify oncology II proteins
 is_oncology_2 = ismember(lower(protein_ids_str), lower(oncology_2_names));
% identify inflammation proteins
 is_inflammation = ismember(lower(protein_ids_str), lower(inflammation_names));
% identify inflammation II proteins
 is_inflammation_2 = ismember(lower(protein_ids_str), lower(inflammation_2_names));

 
 %Nightingale groups
protein_group = repmat("Uncategorized", length(protein_ids), 1);
protein_group(is_inflammation | is_inflammation_2) = "Inflammation";
protein_group(is_oncology | is_oncology_2) = "Oncology";
protein_group(is_cardiometabolic | is_cardiometabolic_2) = "Cardiometabolic";
protein_group(is_neurology | is_neurology_2) = "Neurology";
protein_ids = protein_ids(:);
PLS_table = table(protein_ids, protein_group, plsweights1, plsweights2, plsweights3, ...
    'VariableNames', {'Protein', ...
                     'Category', ......
                     'PLS1_weight', ...
                      'PLS2_weight', 'PLS3_weight'});
writetable(PLS_table, 'PLS_protein_weights_i0proteomics_i2cognitive.csv');

% Make sure all category vectors are columns
is_oncology = is_oncology(:);
is_oncology_2 = is_oncology_2(:);
is_cardiometabolic = is_cardiometabolic(:);
is_cardiometabolic_2 = is_cardiometabolic_2(:);
is_neurology = is_neurology(:);
is_neurology_2 = is_neurology_2(:);
is_inflammation = is_inflammation(:);
is_inflammation_2 = is_inflammation_2(:);

%VOLCANO PLOTS FOR PLS 2 WITH TRANSPARENCY FOR NON SIG
 p_pls2 = 2 * (1 - normcdf(abs(plsweights2)));
% Define significance threshold
sig_threshold = 2 * (1 - normcdf(3));

% Protein categories and colors
protein_categories = {
    is_oncology,           [0.80 0.00 0.70];
    is_oncology_2,         [0.80 0.00 0.70];
    is_cardiometabolic,    [0.10 0.30 0.85];
    is_cardiometabolic_2,  [0.10 0.30 0.85];
    is_neurology,          [0.00 0.60 0.20];
    is_neurology_2,        [0.00 0.60 0.20];
    is_inflammation,       [0.85 0.10 0.10];
    is_inflammation_2,     [0.85 0.10 0.10];
};

figure; hold on;

for i = 1:size(protein_categories,1)
    idx = protein_categories{i,1};
    color = protein_categories{i,2};
    
    % Points above significance
    sig_idx = idx & (p_pls2 <= sig_threshold);
    scatter(stats.W(sig_idx,2), -log10(p_pls2(sig_idx)), 25, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);

    % Points below significance
    nonsig_idx = idx & (p_pls2 > sig_threshold);
    scatter(stats.W(nonsig_idx,2), -log10(p_pls2(nonsig_idx)), 25, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
end

yline(-log10((1 - normcdf(3))), 'k--', 'LineWidth', 1.2); 
%xlabel('PLS Weights (I0-I2 Component 2 Age Regressed)');
ylabel('-log_{10}(p)');
 

%VOLCANO PLOTS FOR PLS 3 WITH TRANSPARENCY FOR NON SIG
 p_pls3 = 2 * (1 - normcdf(abs(plsweights3)));
% Define significance threshold
sig_threshold = 2*(1 - normcdf(3));

% Protein categories and colors
protein_categories = {
    is_oncology,           [0.80 0.00 0.70];
    is_oncology_2,         [0.80 0.00 0.70];
    is_cardiometabolic,    [0.10 0.30 0.85];
    is_cardiometabolic_2,  [0.10 0.30 0.85];
    is_neurology,          [0.00 0.60 0.20];
    is_neurology_2,        [0.00 0.60 0.20];
    is_inflammation,       [0.85 0.10 0.10];
    is_inflammation_2,     [0.85 0.10 0.10];
};

figure; hold on;

for i = 1:size(protein_categories,1)
    idx = protein_categories{i,1};
    color = protein_categories{i,2};
    
    % Points above significance
    sig_idx = idx & (p_pls3 <= sig_threshold);
    scatter(stats.W(sig_idx,3), -log10(p_pls3(sig_idx)), 25, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);

    % Points below significance
    nonsig_idx = idx & (p_pls3 > sig_threshold);
    scatter(stats.W(nonsig_idx,3), -log10(p_pls3(nonsig_idx)), 25, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
end

yline(-log10(2*(1 - normcdf(3))), 'k--', 'LineWidth', 1.2); 
%xlabel('PLS Weights (I0-I2 Component 3 Age Regressed)');
ylabel('-log_{10}(p)');






%FIND KEY PROTEINS FOR POST PROCESSING ON POWERPOINT
sig_threshold = (1 - normcdf(3));
is_sig = p_pls3 <= sig_threshold;

up_regulated   = stats.W(:,3) > 0;
down_regulated = stats.W(:,3) < 0;

inflam_all = is_inflammation | is_inflammation_2;
oncol_all = is_oncology | is_oncology_2;
cardio_all = is_cardiometabolic | is_cardiometabolic_2;
neuro_all = is_neurology | is_neurology_2;

extract_proteins = @(mask, category_name) table(...
    protein_ids(mask & is_sig), ...
    repmat({category_name}, sum(mask & is_sig),1), ...
    stats.W(mask & is_sig,3), ...
    p_pls3(mask & is_sig), ...
    'VariableNames', {'Protein','Category','PLS_Weight','pValue'});

%% extract tables for each category
tbl_oncol  = extract_proteins(oncol_all, 'Oncology');
tbl_cardio = extract_proteins(cardio_all, 'Cardiometabolic');
tbl_neuro  = extract_proteins(neuro_all, 'Neurology');
tbl_inflam = extract_proteins(inflam_all, 'Inflammation');
all_tbl = [tbl_inflam; tbl_oncol; tbl_cardio; tbl_neuro];

%% separate upregulated and downregulated
up_tbl   = all_tbl(all_tbl.PLS_Weight > 0, :);
down_tbl = all_tbl(all_tbl.PLS_Weight < 0, :);

%% sort by |PLS weight| (strongest effect at top)
up_tbl   = sortrows(up_tbl, 'PLS_Weight', 'descend');  
down_tbl = sortrows(down_tbl, 'PLS_Weight');           

%% Display results
disp('Upregulated Proteins:');
disp(up_tbl);

disp('Downregulated Proteins:');
disp(down_tbl);

writetable(up_tbl, 'i0proteomics_i2cognitive_pls3_significant_upregulated_proteins.csv');
writetable(down_tbl, 'i0proteomics_i2cognitive_pls3_significant_downregulated_proteins.csv');

%SCATTER PLOT FOR PLS 3 XS and YS 
figure;
ix1=3; ix2=3
scatter(XS(:,ix1), YS(:,ix2), 5,'k', 'filled'); 
mdl= fitlm(XS(:, ix1), YS(:,ix2)); [ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
plot(XS(:, ix1), ypred, 'k', 'LineWidth', 2);
plot(XS(:, ix1), yci(:,1), 'b', 'LineWidth', 0.5);
plot(XS(:, ix1), yci(:,2), 'b', 'LineWidth', 0.5);
xlabel(sprintf('X Scores (LV%d)', ix1), 'FontSize', 12);
ylabel(sprintf('Y Scores (LV%d)', ix2), 'FontSize', 12);
corr(XS(:,ix1), YS(:,ix2))


%SCATTER PLOT FOR PLS 2 XS and YS 
figure;
ix1=2; ix2=2
scatter(XS(:,ix1), YS(:,ix2), 5,'k', 'filled'); 
mdl= fitlm(XS(:, ix1), YS(:,ix2)); [ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
plot(XS(:, ix1), ypred, 'k', 'LineWidth', 2);
plot(XS(:, ix1), yci(:,1), 'b', 'LineWidth', 0.5);
plot(XS(:, ix1), yci(:,2), 'b', 'LineWidth', 0.5);
%xlabel(sprintf('X Scores (LV%d)', ix1), 'FontSize', 12);
%ylabel(sprintf('Y Scores (LV%d)', ix2), 'FontSize', 12);
corr(XS(:,ix1), YS(:,ix2))

%SCATTER PLOT FOR PLS 1 XS and YS 
figure;
ix1=1; ix2=1
scatter(XS(:,ix1), YS(:,ix2), 5,'k', 'filled'); 
mdl= fitlm(XS(:, ix1), YS(:,ix2)); [ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
plot(XS(:, ix1), ypred, 'k', 'LineWidth', 2);
plot(XS(:, ix1), yci(:,1), 'b', 'LineWidth', 0.5);
plot(XS(:, ix1), yci(:,2), 'b', 'LineWidth', 0.5);
%xlabel(sprintf('X Scores (LV%d)', ix1), 'FontSize', 12);
%ylabel(sprintf('Y Scores (LV%d)', ix2), 'FontSize', 12);
corr(XS(:,ix1), YS(:,ix2))

%BARS FOR UPREGULATED PROTEINS
inflam_all = is_inflammation | is_inflammation_2;
oncol_all  = is_oncology | is_oncology_2;
cardio_all = is_cardiometabolic | is_cardiometabolic_2;
neuro_all  = is_neurology | is_neurology_2;

is_sig = p_pls3 <= sig_threshold;
is_pos = stats.W(:,3) > 0;

num_pos_inflammation_sig = sum(inflam_all & is_sig & is_pos)
num_pos_oncology_sig = sum(oncol_all & is_sig & is_pos)
num_pos_cardiometabolic_sig = sum(cardio_all & is_sig & is_pos)
num_pos_neurology_sig = sum(neuro_all & is_sig & is_pos)
values = [num_pos_inflammation_sig num_pos_oncology_sig num_pos_cardiometabolic_sig num_pos_neurology_sig];

figure;
b = bar(1, values, 'stacked');   % one bar at x = 1

% Assign 4 distinct colors
b(1).FaceColor = [0.85 0.10 0.10];  % inflammation
b(2).FaceColor = [0.80 0.00 0.70];  % oncology
b(3).FaceColor = [0.10 0.30 0.85];  % cardiometabolic
b(4).FaceColor = [0.00 0.60 0.20];  % neurology

ax = gca;
ax.XColor = 'none';
ax.YColor = 'none';
ax.Box = 'off';

%BARS FOR DOWNREGULATED PROTEINS
inflam_all = is_inflammation | is_inflammation_2;
oncol_all  = is_oncology | is_oncology_2;
cardio_all = is_cardiometabolic | is_cardiometabolic_2;
neuro_all  = is_neurology | is_neurology_2;

is_sig = p_pls3 <= sig_threshold;
is_neg = stats.W(:,3) < 0;

num_neg_inflammation_sig = sum(inflam_all & is_sig & is_neg)
num_neg_oncology_sig = sum(oncol_all & is_sig & is_neg)
num_neg_cardiometabolic_sig = sum(cardio_all & is_sig & is_neg)
num_neg_neurology_sig = sum(neuro_all & is_sig & is_neg)
values = [num_neg_inflammation_sig num_neg_oncology_sig num_neg_cardiometabolic_sig num_neg_neurology_sig];

figure;
b = bar(1, values, 'stacked');   % one bar at x = 1

% Assign 4 distinct colors
b(1).FaceColor = [0.85 0.10 0.10];  % inflammation
b(2).FaceColor = [0.80 0.00 0.70];  % oncology
b(3).FaceColor = [0.10 0.30 0.85];  % cardiometabolic
b(4).FaceColor = [0.00 0.60 0.20];  % neurology

ax = gca;
ax.XColor = 'none';
ax.YColor = 'none';
ax.Box = 'off';

%SCATTER PLOT FOR TFF1 AND TMT
i = find(strcmp(protein_ids, 'tff1'))
x_tff1 = x(:,i)
y = Y(:,1)
figure;
dscatter(x_tff1, y);
lsline;
corr(x_tff1, y)

%SCATTER PLOT FOR TFF1 AND GF
i = find(strcmp(protein_ids, 'tff1'))
x_tff1 = x(:,i)
y = Y(:,2)
figure;
dscatter(x_tff1, y);
lsline;
corr(x_tff1, y)

%SCATTER PLOT FOR TFF1 AND PAL
i = find(strcmp(protein_ids, 'tff1'))
x_tff1 = x(:,i)
y = Y(:,3)
figure;
dscatter(x_tff1, y);
lsline;
corr(x_tff1, y)

%SCATTER PLOT FOR TFF1 AND SDS
i = find(strcmp(protein_ids, 'tff1'))
x_tff1 = x(:,i)
y = Y(:,4)
figure;
dscatter(x_tff1, y);
lsline;
corr(x_tff1, y)

%SCATTER PLOT FOR FASLG AND TMT
i = find(strcmp(protein_ids, 'faslg'))
x_faslg = x(:,i)
y = Y(:,1)
figure;
dscatter(x_faslg, y);
lsline;
corr(x_faslg, y)

%SCATTER PLOT FOR FASLG AND GF
i = find(strcmp(protein_ids, 'faslg'))
x_faslg = x(:,i)
y = Y(:,2)
figure;
dscatter(x_faslg, y);
lsline;
corr(x_faslg, y)

%SCATTER PLOT FOR FASLG AND PAL
i = find(strcmp(protein_ids, 'faslg'))
x_faslg = x(:,i)
y = Y(:,3)
figure;
dscatter(x_faslg, y);
lsline;
corr(x_faslg, y)

%SCATTER PLOT FOR FASLG AND SDS
i = find(strcmp(protein_ids, 'faslg'))
x_faslg = x(:,i)
y = Y(:,4)
figure;
dscatter(x_faslg, y);
lsline;
corr(x_faslg, y)

%SCATTER PLOT FOR ESAM AND TMT
i = find(strcmp(protein_ids, 'esam'))
x_notch3 = x(:,i)
y = Y(:,1)
figure;
dscatter(x_notch3, y);
lsline;
corr(x_notch3, y)

%SCATTER PLOT FOR ESAM AND SDS
i = find(strcmp(protein_ids, 'esam'))
x_notch3 = x(:,i)
y = Y(:,4)
figure;
dscatter(x_notch3, y);
lsline;
corr(x_notch3, y)

%SCATTER PLOT FOR MMP12 AND TMT
i = find(strcmp(protein_ids, 'mmp12'))
x_mmp12 = x(:,i)
y = Y(:,1)
figure;
dscatter(x_mmp12, y);
lsline;
corr(x_mmp12, y)

%SCATTER PLOT FOR MMP12 AND SDS
i = find(strcmp(protein_ids, 'mmp12'))
x_mmp12 = x(:,i)
y = Y(:,4)
figure;
dscatter(x_mmp12, y);
lsline;
corr(x_mmp12, y)


%FIND KEY PROTEINS FOR POST PROCESSING ON POWERPOINT
up_regulated = stats.W(:,2) > 3;
down_regulated = stats.W(:,2) < -3;

inflam_all = is_inflammation | is_inflammation_2;
oncol_all = is_oncology | is_oncology_2;
cardio_all = is_cardiometabolic | is_cardiometabolic_2;
neuro_all = is_neurology | is_neurology_2;

extract_proteins = @(mask, category_name) table(...
    protein_ids(mask & is_sig), ...
    repmat({category_name}, sum(mask & is_sig),1), ...
    stats.W(mask & is_sig,2), ...
    'VariableNames', {'Protein','Category','PLS_Weight'});

%% extract tables for each category
tbl_oncol  = extract_proteins(oncol_all, 'Oncology');
tbl_cardio = extract_proteins(cardio_all, 'Cardiometabolic');
tbl_neuro  = extract_proteins(neuro_all, 'Neurology');
tbl_inflam = extract_proteins(inflam_all, 'Inflammation');
all_tbl = [tbl_inflam; tbl_oncol; tbl_cardio; tbl_neuro];

%% separate upregulated and downregulated
up_tbl   = all_tbl(all_tbl.PLS_Weight > 0, :);
down_tbl = all_tbl(all_tbl.PLS_Weight < 0, :);

%% sort by |PLS weight| (strongest effect at top)
up_tbl   = sortrows(up_tbl, 'PLS_Weight', 'descend');  
down_tbl = sortrows(down_tbl, 'PLS_Weight');           

%% Display results
disp('Upregulated Proteins:');
disp(up_tbl);

disp('Downregulated Proteins:');
disp(down_tbl);

writetable(up_tbl, 'i0proteomics_i2cognitive_pls2_significant_upregulated_proteins_no_age.csv');
writetable(down_tbl, 'i0proteomics_i2cognitive_pls2_significant_downregulated_proteins_no_age.csv');

%SCATTER PLOT FOR PLS 2 XS and YS 
figure;
ix1=2; ix2=2
scatter(XS(:,ix1), YS(:,ix2), 5,'k', 'filled'); 
mdl= fitlm(XS(:, ix1), YS(:,ix2)); [ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
plot(XS(:, ix1), ypred, 'k', 'LineWidth', 2);
plot(XS(:, ix1), yci(:,1), 'b', 'LineWidth', 0.5);
plot(XS(:, ix1), yci(:,2), 'b', 'LineWidth', 0.5);
xlabel(sprintf('X Scores (LV%d)', ix1), 'FontSize', 12);
ylabel(sprintf('Y Scores (LV%d)', ix2), 'FontSize', 12);
corr(XS(:,ix1), YS(:,ix2))

%BARS FOR UPREGULATED PROTEINS
inflam_all = is_inflammation | is_inflammation_2;
oncol_all  = is_oncology | is_oncology_2;
cardio_all = is_cardiometabolic | is_cardiometabolic_2;
neuro_all  = is_neurology | is_neurology_2;

is_sig = p_pls2 <= sig_threshold;
is_pos = stats.W(:,2) > 0;

num_pos_inflammation_sig = sum(inflam_all & is_sig & is_pos)
num_pos_oncology_sig = sum(oncol_all & is_sig & is_pos)
num_pos_cardiometabolic_sig = sum(cardio_all & is_sig & is_pos)
num_pos_neurology_sig = sum(neuro_all & is_sig & is_pos)
values = [num_pos_inflammation_sig num_pos_oncology_sig num_pos_cardiometabolic_sig num_pos_neurology_sig];

figure;
b = bar(1, values, 'stacked');   % one bar at x = 1

% Assign 4 distinct colors
b(1).FaceColor = [0.85 0.10 0.10];  % inflammation
b(2).FaceColor = [0.80 0.00 0.70];  % oncology
b(3).FaceColor = [0.10 0.30 0.85];  % cardiometabolic
b(4).FaceColor = [0.00 0.60 0.20];  % neurology

ax = gca;
ax.XColor = 'none';
ax.YColor = 'none';
ax.Box = 'off';

%BARS FOR DOWNREGULATED PROTEINS
inflam_all = is_inflammation | is_inflammation_2;
oncol_all  = is_oncology | is_oncology_2;
cardio_all = is_cardiometabolic | is_cardiometabolic_2;
neuro_all  = is_neurology | is_neurology_2;

is_sig = p_pls2 <= sig_threshold;
is_neg = stats.W(:,2) < 0;

num_neg_inflammation_sig = sum(inflam_all & is_sig & is_neg)
num_neg_oncology_sig = sum(oncol_all & is_sig & is_neg)
num_neg_cardiometabolic_sig = sum(cardio_all & is_sig & is_neg)
num_neg_neurology_sig = sum(neuro_all & is_sig & is_neg)
values = [num_neg_inflammation_sig num_neg_oncology_sig num_neg_cardiometabolic_sig num_neg_neurology_sig];

figure;
b = bar(1, values, 'stacked');   % one bar at x = 1

% Assign 4 distinct colors
b(1).FaceColor = [0.85 0.10 0.10];  % inflammation
b(2).FaceColor = [0.80 0.00 0.70];  % oncology
b(3).FaceColor = [0.10 0.30 0.85];  % cardiometabolic
b(4).FaceColor = [0.00 0.60 0.20];  % neurology

ax = gca;
ax.XColor = 'none';
ax.YColor = 'none';
ax.Box = 'off';

