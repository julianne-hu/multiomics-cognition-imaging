clear
cd /Users/juliannehu/Documents/UKB_proteomics/data/raw
mri=readtable('mri_v2.csv');
demo=readtable('demographics_v2.csv');
cd /Users/juliannehu/Documents/UKB_proteomics/data
demo_bmi=readtable('ukb_demo.csv');
demo_other_covariates=readtable('ukb_demo_ses_edu_cvd.csv');
demo_all_covariates=innerjoin(demo_bmi,demo_other_covariates);
demo=outerjoin(demo,demo_all_covariates);
demo.Properties.VariableNames{'eid_demo'} = 'eid';
mri=innerjoin(mri, demo);

cd /Users/juliannehu/Documents/UKB_proteomics/data/dictionaries
mri_dictionary=readtable('mri_dictionary.csv');

cd /Users/juliannehu/Documents/UKB_proteomics/data/cleaned
proteomics=readtable('mean_imputed_olink_i0.csv');
merged_data=innerjoin(proteomics,mri);
clear mri proteomics demo;

%remove outlier participant eid 3814881
outlier = 29675;  % the row where the participant to remove is found
% Remove the row
merged_data(outlier, :) = [];


i0age=merged_data.p21003_i0; 
i2age=merged_data.p21003_i2; 
sex=merged_data.p31; 
site=merged_data.p54_i2; 
motion=merged_data.p24419_i2; 
icv=merged_data.p26521_i2;
merged_data.delta_time=i2age-i0age; 
delta_time=merged_data.delta_time;
eid=merged_data.eid;
bmi=merged_data.x21001_0_0;
bmi(bmi < 0) = NaN;
education=merged_data.x6138_0_0;
education(education < 0) = NaN;
cvd_cols = merged_data(:, 1672:1825);
cvd = double(any(~ismissing(cvd_cols), 2));
ses=merged_data.x738_0_0;
ses(ses < 0) = NaN;

%%% requires tabular_ukb to have ran
%ix=(cellfun('isempty', clinical)' | sum(isnan(ms_ordered)')'>0 ) ; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
%ix=(cellfun('isempty', clinical)'| isnan(ica_partial(:,1))) ; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
%cognitive=cognitive; cognitive(ix,:)=[]; minimal=minimal_ordered; minimal(ix,:)=[]; MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[];


%X=ms_ordered; X(ix,:)=[]; 
%X=merged_data{:,2:2924}; 
protein_ids=merged_data.Properties.VariableNames(2:1460);
mri_ids=merged_data.Properties.VariableNames(1463:1608);
cols = 1463:1608;
[isFound, loc] = ismember(mri_ids, mri_dictionary.name);
mri_names = cell(size(mri_ids));
mri_names(isFound) = mri_dictionary.title(loc(isFound));
mri_names = matlab.lang.makeValidName(mri_names);
new_headers(cols) = mri_names;
merged_data.Properties.VariableNames(1463:1608) = new_headers(cols);


%protein_ids=merged_data.Properties.VariableNames(2:2924);
%m=nanmean(X);

%for i=1:length(protein_ids)
%    X(isnan(X(:,i)),i)=m(i);
%end
%X = merged_data{:, 2:1455};
X = merged_data{:, 2:1460};

Y=merged_data{:,1463:2:1608};
cols=1463:2:1608;
cols_exclude=[1483 1484];
cols(ismember(cols, cols_exclude))= [];
Y= merged_data{:, cols};
mri_variables=merged_data.Properties.VariableNames(cols)';
%%% remove nans
naninx=sum(isnan(Y)')' | isnan(bmi) | isnan(education) | isnan(ses); 
Y=Y(naninx==0,:);  
X=X(naninx==0,:); 
sex=sex(naninx==0);
site=site(naninx==0);
motion=motion(naninx==0);
icv=icv(naninx==0); %minimal=minimal(naninx==0,:);c=c(naninx==0);MDD_prs=MDD_prs(naninx==0);
delta_time=delta_time(naninx==0);
i0age=i0age(naninx==0)
i2age=i2age(naninx==0)
bmi=bmi(naninx==0);
education=education(naninx==0);
cvd=cvd(naninx==0);
ses=ses(naninx==0);

Y=zscore(Y); X=zscore(X); 
ncomp=3

% defining other covariates
i0age_squared = i0age.^2;
i2age_squared = i2age.^2;
sex=categorical(sex);
sex_num = double(sex == 'Male');
i0age_by_sex = i0age .* sex_num;
i2age_by_sex = i2age .* sex_num;

%%% regress out sex and site and potentially other covariates
clear x; 
for i=1:length(X(1,:)); 
    mdl = fitlm(table(delta_time, sex, site, bmi, education, cvd, ses, X(:,i))); 
    x(:,i)=mdl.Residuals.Raw;
end
clear y; 
for i=1:length(Y(1,:)); 
    mdl = fitlm(table(delta_time, sex, site, motion, icv, bmi, education, cvd, ses, Y(:,i))); 
    y(:,i)=mdl.Residuals.Raw;
end
Y=y;

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);PCTVAR
%%% explore the PLS components and their correlation with Y
%XS VS YS FIGURE:%               figure(31); scatter(XS(ismember(c','dep'),1), YS(ismember(c','dep'),1), 5,'b', 'filled'); hold on;scatter(XS(ismember(c','depanx'),1), YS(ismember(c','depanx'),1), 5, 'r', 'filled'); hold on;scatter(XS(ismember(c','anx'),1), YS(ismember(c','anx'),1), 5, 'k', 'filled'); mdl= fitlm(XS(:, 1), YS(:,1)); [ypred,yci] = predict(mdl,XS(:, 1), 'Alpha',0.001); hold on; plot(XS(:, 1), ypred, 'k', 'LineWidth', 2);plot(XS(:, 1), yci(:,1), 'k', 'LineWidth', 0.5);plot(XS(:, 1), yci(:,2), 'k', 'LineWidth', 0.5);
% and their distributions: %        figure; plot_histogram_shaded(XS,'Alpha',0.3,'color',[0.5 0.5 0.5], 'Normalization', 'pdf');figure; plot_histogram_shaded(YS,'Alpha',0.3,'color',[0.5 0.5 0.5], 'Normalization', 'pdf');
figure(1);imagesc(corr(XS(:,1:3),YS)); colormap bone ; colorbar; corr(XS, Y) 
combs=allcomb([1:3], [1:3]);figure(13); 
for i=1:length(combs)
hold off; ix1=combs(i,1);ix2=combs(i,2); subplot(3,4,i); 
scatter(XS(:,ix1), YS(:,ix2), 5,'b', 'filled'); 
hold on;scatter(XS(:,ix1), YS(:,ix2), 5, 'r', 'filled')
hold on;scatter(XS(:,ix1), YS(:,ix2), 5, 'k', 'filled');
%mdl= fitlm(XS(:, ix1), Y(:,ix2)); [ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
%plot(XS(:, ix1), ypred, 'k', 'LineWidth', 2);
%plot(XS(:, ix1), yci(:,1), 'k', 'LineWidth', 0.5);
%plot(XS(:, ix1), yci(:,2), 'k', 'LineWidth', 0.5);
end; clear mdl ypred yci mdl ix1 ix2 combs

% checking how many standard deviations away the outlier is
max(YS(:,1)) / std(YS(:,1))
mean(YS(:,1))


%corr(age, XS)
max(corr(X, Y))

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
xlim([0.005 0.052]); 
xline(0.0511, 'r', 'LineWidth', 2); 
% Only show x and y axes
ax = gca;
ax.Box = 'off';          % removes top and right lines
ax.TickDir = 'out';      % ticks point outward
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
sum(PCTVAR')

%% bootstrapping to get the func connectivity weights for PLS1, 2 and 3
dim=3
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,dim);PCTVAR
PLS1w=stats.W(:,1);
PLS2w=stats.W(:,2);
PLS3w=stats.W(:,2);

bootnum=5000;
PLS1weights=[];
PLS2weights=[];
PLS3weights=[];

parfor i=1:bootnum
    i
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

PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');
PLS3sw=std(PLS3weights');

plsweights1=PLS1w./PLS1sw';
plsweights2=PLS2w./PLS2sw'; 
plsweights3=PLS3w./PLS3sw';

sum(plsweights1>3)
sum(plsweights1<-3)

sum(plsweights3>3)
sum(plsweights3<-3)

protein_ids(plsweights1>3)'
protein_ids(plsweights1<-3)'

protein_ids(plsweights2>3)'
protein_ids(plsweights2<-3)'

protein_ids(plsweights3>3)'
protein_ids(plsweights3<-3)' 


%Nightingale groups
protein_group = repmat("Uncategorized", length(protein_ids), 1);

protein_group(is_inflammation | is_inflammation_2)       = "Inflammation";
protein_group(is_oncology | is_oncology_2)               = "Oncology";
protein_group(is_cardiometabolic | is_cardiometabolic_2) = "Cardiometabolic";
protein_group(is_neurology | is_neurology_2)             = "Neurology";

protein_ids = protein_ids(:);
PLS_table = table(protein_ids, protein_group, plsweights1, plsweights2, plsweights3, ...
    'VariableNames', {'Protein', ...
                     'Category', ......
                     'PLS1_weight', ...
                      'PLS2_weight', 'PLS3_weight'});
writetable(PLS_table, 'PLS_protein_weights_i0proteomics_i2cognitive.csv');







protein_ids = protein_ids(:);
PLS_table = table(protein_ids, plsweights1, plsweights2, plsweights3, ...
    'VariableNames', {'Protein', ...
                      'PLS1_weight', ...
                      'PLS2_weight', 'PLS3_weight'});
writetable(PLS_table, 'PLS_protein_weights_i0proteomics_i2mri_with_age.csv');


%summary_stats=table;
%summary_stats.protein_name=protein_ids;
%summary_stats.plsweights1=plsweights1;
%summary_stats.beta_weights=stats.W(:,1);
%writetable(summary_stats, '')

%%
%figure; scatter(merged_data.gfap, merged_data.MeanThicknessOfInsula_rightHemisphere__Instance2, 'filled','k')
%corr(merged_data.gfap, merged_data.MeanThicknessOfInsula_rightHemisphere__Instance2, 'rows','pairwise')
%corr(merged_data.nefl, merged_data.MeanThicknessOfInsula_rightHemisphere__Instance2, 'rows','pairwise')
%corr(merged_data.tnfrsf10a, merged_data.MeanThicknessOfInsula_rightHemisphere__Instance2, 'rows','pairwise')
%corr(merged_data.tnfrsf11b, merged_data.MeanThicknessOfInsula_rightHemisphere__Instance2, 'rows','pairwise')

%FIGURE WITH OLINK EXPLORE
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
% 
protein_ids_str = string(protein_ids);
% 
% % identify cardiometabolic proteins
 is_cardiometabolic = ismember(lower(protein_ids_str), lower(cardiometabolic_names));
% % identify cardiometabolic II proteins
 is_cardiometabolic_2 = ismember(lower(protein_ids_str), lower(cardiometabolic_2_names));
% % identify neurology proteins
 is_neurology = ismember(lower(protein_ids_str), lower(neurology_names));
% % identify neurology II proteins
 is_neurology_2 = ismember(lower(protein_ids_str), lower(neurology_2_names));
% % identify oncology proteins
 is_oncology = ismember(lower(protein_ids_str), lower(oncology_names));
% % identify oncology II proteins
 is_oncology_2 = ismember(lower(protein_ids_str), lower(oncology_2_names));
 % % identify inflammation proteins
 is_inflammation = ismember(lower(protein_ids_str), lower(inflammation_names));
% % identify inflammation II proteins
 is_inflammation_2 = ismember(lower(protein_ids_str), lower(inflammation_2_names));
% 


%VOLCANO PLOTS FOR PLS 1 WITH TRANSPARENCY FOR NON SIG
p_pls1 = 2 * (1 - normcdf(abs(plsweights1)));
 % Define significance threshold
sig_threshold = (1 - normcdf(3));
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
   
    % Force column vector
    idx = idx(:);
    p_pls1 = p_pls1(:);

    % Skip if idx has no true values
    if ~any(idx)
        continue;
    end

    % Points above significance
    sig_idx = idx & (p_pls1 <= sig_threshold);
    scatter(stats.W(sig_idx,1), -log10(p_pls1(sig_idx)), 25, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);

    % Points below significance
    nonsig_idx = idx & (p_pls1 > sig_threshold);
    scatter(stats.W(nonsig_idx,1), -log10(p_pls1(nonsig_idx)), 25, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
end

yline(-log10(1 - normcdf(3)), 'k--', 'LineWidth', 1.2); 
%xlabel('PLS Weights (I0-I2 Component 3 Age Regressed)');
ylabel('-log_{10}(p)');

% Highlight protein nefl in yellow
target = "il6";
idx = strcmp(protein_ids, target);
scatter(stats.W(idx,1), -log10(p_pls1(idx)), ...
        120, 'y', 'filled', 'MarkerEdgeColor','k');   % yellow + black edge
% Add datatips to all *existing* scatter points
dtt = dataTipTextRow('Protein', protein_ids);
allScatters = findobj(gca,'Type','Scatter');
for s = 1:length(allScatters)
    allScatters(s).DataTipTemplate.DataTipRows(end+1) = dtt;
end



%VOLCANO PLOTS FOR PLS 2 WITH TRANSPARENCY FOR NON SIG
 p_pls2 = 2 * (1 - normcdf(abs(plsweights2)));
 % Define significance threshold
sig_threshold = (1 - normcdf(3));
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
    % Force column vector
    idx = idx(:);
    p_pls1 = p_pls1(:);

    % Skip if idx has no true values
    if ~any(idx)
        continue;
    end
    % Points above significance
    sig_idx = idx & (p_pls2 <= sig_threshold);
    scatter(stats.W(sig_idx,2), -log10(p_pls2(sig_idx)), 25, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);

    % Points below significance
    nonsig_idx = idx & (p_pls2 > sig_threshold);
    scatter(stats.W(nonsig_idx,2), -log10(p_pls2(nonsig_idx)), 25, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
end

yline(-log10(1 - normcdf(3)), 'k--', 'LineWidth', 1.2); 
%xlabel('PLS Weights (I0-I2 Component 3 Age Regressed)');
ylabel('-log_{10}(p)');



%VOLCANO PLOTS FOR PLS 3 WITH TRANSPARENCY FOR NON SIG
 p_pls3 = 2 * (1 - normcdf(abs(plsweights3)));
 % Define significance threshold
sig_threshold = (1 - normcdf(3));
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
    % Force column vector
    idx = idx(:);
    p_pls1 = p_pls1(:);

    % Skip if idx has no true values
    if ~any(idx)
        continue;
    end
    % Points above significance
    sig_idx = idx & (p_pls3 <= sig_threshold);
    scatter(stats.W(sig_idx,3), -log10(p_pls3(sig_idx)), 25, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);

    % Points below significance
    nonsig_idx = idx & (p_pls3 > sig_threshold);
    scatter(stats.W(nonsig_idx,3), -log10(p_pls3(nonsig_idx)), 25, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
end

yline(-log10(1 - normcdf(3)), 'k--', 'LineWidth', 1.2); 
%xlabel('PLS Weights (I0-I2 Component 3 Age Regressed)');
ylabel('-log_{10}(p)');

%FIND KEY PROTEINS FOR POST PROCESSING ON POWERPOINT
sig_threshold = 2*(1 - normcdf(3));
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

%SCATTER PLOT FOR PLS 1 XS AND YS WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(age_i2-min(age_i2))./max(age_i2-min(age_i2));
colormap_custom(:,3)=1-age_i2./max(age_i2);
figure;
ix1=1; ix2=1
scatter(XS(:,ix1), YS(:,ix2), 5, colormap_custom, 'filled'); 
mdl= fitlm(XS(:, ix1), YS(:,ix2)); [ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
plot(XS(:, ix1), ypred, 'k', 'LineWidth', 2);
plot(XS(:, ix1), yci(:,1), 'b', 'LineWidth', 0.5);
plot(XS(:, ix1), yci(:,2), 'b', 'LineWidth', 0.5);
corr(XS(:,ix1), YS(:,ix2))
% Add colorbar 
% Build a 256-level version of your SAME color gradient
age_norm = linspace(0,1,256)';
cmap_bar = [ ...
    (age_norm - min(age_norm)) ./ max(age_norm - min(age_norm)), ... % red
    zeros(256,1), ...
    1 - age_norm ./ max(age_norm) ... % blue 
];
% Create colorbar
hold on;
h = imagesc([max(XS(:,ix1))+1 max(XS(:,ix1))+2], ...
            [min(YS(:,ix2)) max(YS(:,ix2))], age_norm);
set(h,'Visible','off');
colormap(cmap_bar);
c = colorbar;
%c.Label.String = 'Age (years)';
c.Ticks = linspace(0,1,5);
c.TickLabels = round(linspace(min(age_i2), max(age_i2), 5));




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

%%% IDENTIFYING OUTLIER BASED ON THE ABOVE SCATTER PLOT
suspicious_x = 0.00516082;
suspicious_y = 2801.32;
diffs = abs(XS(:,ix1) - suspicious_x) + abs(YS(:,ix2) - suspicious_y);
[~, idx_outlier] = min(diffs) %row 5061 is the outlier

%Replot scatter without participant 5061 to confirm
idx_outlier = 5061;
keep_idx = setdiff(1:size(XS,1), idx_outlier);
figure;
scatter(XS(keep_idx, ix1), YS(keep_idx, ix2), 8, 'k', 'filled');
hold on;
% Refit model without the outlier
mdl2 = fitlm(XS(keep_idx, ix1), YS(keep_idx, ix2));
[ypred2, yci2] = predict(mdl2, XS(keep_idx, ix1), 'Alpha', 0.001);
plot(XS(keep_idx, ix1), ypred2, 'k', 'LineWidth', 2);
plot(XS(keep_idx, ix1), yci2(:,1), 'b', 'LineWidth', 0.5);
plot(XS(keep_idx, ix1), yci2(:,2), 'b', 'LineWidth', 0.5);
title('Scatter without participant 5061');

%figure out which participant row 5061 corresponds to
eid=eid(naninx==0); 
eid(5061) %participant eid 3814881







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



%SCATTER PLOT FOR PLS 3 XS and YS 
figure;
ix1=3; ix2=3
scatter(XS(:,ix1), YS(:,ix2), 5,'k', 'filled'); 
mdl= fitlm(XS(:, ix1), YS(:,ix2)); [ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
plot(XS(:, ix1), ypred, 'k', 'LineWidth', 2);
plot(XS(:, ix1), yci(:,1), 'b', 'LineWidth', 0.5);
plot(XS(:, ix1), yci(:,2), 'b', 'LineWidth', 0.5);
%xlabel(sprintf('X Scores (LV%d)', ix1), 'FontSize', 12);
%ylabel(sprintf('Y Scores (LV%d)', ix2), 'FontSize', 12);
corr(XS(:,ix1), YS(:,ix2))



%SCATTER PLOT FOR NEFL AND YS1 WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(age_i2-min(age_i2))./max(age_i2-min(age_i2));
colormap_custom(:,3)=1-age_i2./max(age_i2);
i = find(strcmp(protein_ids, 'nefl'))
x_nefl = x(:,i)
y = YS(:,1);
figure;
scatter(x_nefl, y, 5, colormap_custom, 'filled'); 
lsline;
corr(x_nefl, y)

%SCATTER PLOT FOR TNFRSF1A AND YS1 WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(age_i2-min(age_i2))./max(age_i2-min(age_i2));
colormap_custom(:,3)=1-age_i2./max(age_i2);
i = find(strcmp(protein_ids, 'tnfrsf1a'))
x_tnfrsf1a = x(:,i)
y = YS(:,1);
figure;
scatter(x_tnfrsf1a, y, 5, colormap_custom, 'filled'); 
lsline;
corr(x_tnfrsf1a, y)


%SCATTER PLOT FOR NEFL AND thickness left superior frontal WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(age_i2-min(age_i2))./max(age_i2-min(age_i2));
colormap_custom(:,3)=1-age_i2./max(age_i2);
i = find(strcmp(protein_ids, 'nefl'))
x_nefl = x(:,i)
j = find(strcmp(mri, 'MeanThicknessOfSuperiorfrontal_leftHemisphere__Instance2'))
y = Y(:,j)
figure;
scatter(x_nefl, y, 5, colormap_custom, 'filled'); 
lsline;
corr(x_nefl, y)


%SCATTER PLOT FOR TNFRSF1A AND thickness left superior frontal WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(age_i2-min(age_i2))./max(age_i2-min(age_i2));
colormap_custom(:,3)=1-age_i2./max(age_i2);
i = find(strcmp(protein_ids, 'tnfrsf1a'))
x_tnfrsf1a = x(:,i)
j = find(strcmp(mri, 'MeanThicknessOfSuperiorfrontal_leftHemisphere__Instance2'))
y = Y(:,j)
figure;
scatter(x_tnfrsf1a, y, 5, colormap_custom, 'filled'); 
lsline;
corr(x_tnfrsf1a, y)



%BARS FOR UPREGULATED PROTEINS
inflam_all = is_inflammation | is_inflammation_2;
oncol_all  = is_oncology | is_oncology_2;
cardio_all = is_cardiometabolic | is_cardiometabolic_2;
neuro_all  = is_neurology | is_neurology_2;
is_sig = p_pls1 <= sig_threshold;
is_pos = stats.W(:,1) > 0;
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
is_sig = p_pls1 <= sig_threshold;
is_neg = stats.W(:,1) < 0;
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




% PREPARING DATA FOR PLOTTING RESULTS WITH GGSEG
[r,p]=corr(XS, Y)
r = r.';
p = p.';
mri_variables = mri_variables(:); 
for i = 1:ncomp
    % Extract the correlations and p-values for component c
    correlation = r(:, i);
    pvalue = p(:, i);
    mri = mri_variables;
    %mri = repmat(i, length(mri_variables), 1);  % k x 1
    % Build the table
    mri_results = table(mri, correlation, pvalue, ...
        'VariableNames', {'BrainRegion','Correlation','PValue'});
    % Write to CSV
    writetable(mri_results, sprintf('mri_results_i0i2_pls%d.csv', i));
end




%% MEDIATION MODEL %%
X_mediation = X; 
%X_mediation = XS(:1);
Y_mediation = YS(:,1);
age_i2=merged_data.p21003_i2;
age_i2=age_i2(naninx==0);
M_mediation = age_i2;

[paths, stats] = mediation(X_mediation, Y_mediation, M_mediation, 'plots', 'verbose', 'boot', 'bootreps', 5000);

%SCATTER PLOT FOR NEFL AND thickness left superior frontal WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(age_i2-min(age_i2))./max(age_i2-min(age_i2));
colormap_custom(:,3)=1-age_i2./max(age_i2);
i = find(strcmp(protein_ids, 'paep'))
x_nefl = x(:,i)
j = find(strcmp(mri, 'MeanThicknessOfSuperiorfrontal_leftHemisphere__Instance2'))
y = Y(:,j)
figure;
scatter(x_nefl, y, 5, colormap_custom, 'filled'); 
lsline;
corr(x_nefl, y)



%SCATTER PLOT FOR NEFL AND thickness left superior frontal WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(age_i2-min(age_i2))./max(age_i2-min(age_i2));
colormap_custom(:,3)=1-age_i2./max(age_i2);
i = find(strcmp(protein_ids, 'snap29'))
x_nefl = x(:,i)
j = find(strcmp(mri, 'MeanThicknessOfPrecuneus_leftHemisphere__Instance2'))
y = Y(:,j)
figure;
scatter(x_nefl, y, 5, colormap_custom, 'filled'); 
lsline;
corr(x_nefl, y)

