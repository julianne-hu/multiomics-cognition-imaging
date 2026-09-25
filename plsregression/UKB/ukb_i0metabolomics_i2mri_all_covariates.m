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
metabolomics=readtable('mean_imputed_metabolomic_i0_v2.csv');
merged_data=innerjoin(metabolomics,mri);
clear mri metabolomics demo;

i0age=merged_data.p21003_i0; 
i2age=merged_data.p21003_i2; 
sex=merged_data.p31; 
site=merged_data.p54_i2; 
motion=merged_data.p24419_i2; 
icv=merged_data.p26521_i2; 
merged_data.delta_time=i2age-i0age; 
delta_time=merged_data.delta_time;
bmi=merged_data.x21001_0_0;
bmi(bmi < 0) = NaN;
education=merged_data.x6138_0_0;
education(education < 0) = NaN;
cvd_cols = merged_data(:, 464:617);
cvd = double(any(~ismissing(cvd_cols), 2));
ses=merged_data.x738_0_0;
ses(ses < 0) = NaN;

%%% requires tabular_ukb to have ran
%ix=(cellfun('isempty', clinical)' | sum(isnan(ms_ordered)')'>0 ) ; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
%ix=(cellfun('isempty', clinical)'| isnan(ica_partial(:,1))) ; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
%cognitive=cognitive; cognitive(ix,:)=[]; minimal=minimal_ordered; minimal(ix,:)=[]; MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[];


%X=ms_ordered; X(ix,:)=[]; 
%X=merged_data{:,2:2924}; 
protein_ids=merged_data.Properties.VariableNames(2:252);
mri_ids=merged_data.Properties.VariableNames(255:400);
cols = 255:400;
[isFound, loc] = ismember(mri_ids, mri_dictionary.name);
mri_names = cell(size(mri_ids));
mri_names(isFound) = mri_dictionary.title(loc(isFound));
mri_names = matlab.lang.makeValidName(mri_names);
new_headers(cols) = mri_names;
merged_data.Properties.VariableNames(255:400) = new_headers(cols);


%protein_ids=merged_data.Properties.VariableNames(2:2924);
%m=nanmean(X);

%for i=1:length(protein_ids)
%    X(isnan(X(:,i)),i)=m(i);
%end
%X = merged_data{:, 2:1455};
X = merged_data{:, 2:252};

Y=merged_data{:,255:2:400};
cols=255:2:400;
cols_exclude=[275 276];
cols(ismember(cols, cols_exclude))= [];
Y=merged_data{:, cols};
mri_variables=merged_data.Properties.VariableNames(cols)';

%sex=categorical(sex);
%site=categorical(site);


naninx=any(isnan(Y),2) | any(isnan(X), 2) | any(ismissing(site),2) | any(isnan(motion), 2) | any(isnan(bmi),2) | any(isnan(education),2) | any(isnan(ses),2); %|sum(ismissing(sex)')' | sum(ismissing(icv)')' | sum(ismissing(i0age)')' | sum(ismissing(i2age)')' | sum(ismissing(delta_time)')'
Y=Y(~naninx,:);
X=X(~naninx,:); 
i0age=i0age(~naninx); 
i2age=i2age(~naninx); 
sex=sex(~naninx);
site=site(~naninx); 
motion=motion(~naninx); 
icv=icv(~naninx); 
delta_time=delta_time(~naninx); 
bmi = bmi(~naninx);
education = education(~naninx);
cvd = cvd(~naninx);
ses = ses(~naninx);



%%% remove rows with missing data
%valid_X = any(isnan(X),2)
%valid_Y = any(isnan(Y),2);
%valid_sex = isundefined(sex);
%valid_site = isundefined(site);
%valid_motion= isundefined(motion);
%valid_icv= isundefined(icv);
%valid_i0age= isundefined(i0age);
%valid_i2age= isundefined(i2age);
%valid_delta_time= isundefined(delta_time);
%rowsToRemove = valid_X | valid_Y | valid_sex | valid_site | valid_motion | valid_icv | valid_i0age | valid_i2age | valid_delta_time
%naninx = ~rowsToRemove; 

% Identify rows with any NaN
%rowsWithAnyNaN = any(isnan(Y),2);
% Identify rows that are entirely NaN
%rowsWithAllNaN = all(isnan(Y),2);
% Combine - TRUE if row should be removed
%rowsToRemove = rowsWithAnyNaN | rowsWithAllNaN;
% Keep only rows that are not in rowsToRemove
%naninx = ~rowsToRemove;
%Y=Y(naninx,:);
%X=X(naninx,:);
%sex=sex(naninx);
%site=site(naninx);
%motion=motion(naninx);
%icv=icv(naninx);
%i0age=i0age(naninx);
%i2age=i2age(naninx);
%delta_time=delta_time(naninx);

%example corr for precuneus and LDL without regressing things out 
x_val = merged_data.p23412_i0;
y = merged_data.MeanThicknessOfPrecentral_leftHemisphere__Instance2;
figure;
scatter(x_val, y);
lsline;
corr(x_val, y, 'Rows','pairwise')


%naninx=sum(isnan(Y)')'; 
%Y=Y(naninx==0,:);  
%X=X(naninx==0,:); 
%sex=sex(naninx==0);
%site=site(naninx==0);
%motion=motion(naninx==0);
%icv=icv(naninx==0); 
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
    mdl = fitlm(table(delta_time, i0age, sex, site, i0age_squared, i0age_by_sex, bmi, education, cvd, ses, X(:,i))); 
    x(:,i)=mdl.Residuals.Raw;
end
clear y; 
for i=1:length(Y(1,:)); 
    mdl = fitlm(table(delta_time, i2age, sex, site, i2age_squared, i2age_by_sex, motion, icv, bmi, education, cvd, ses, Y(:,i))); 
    y(:,i)=mdl.Residuals.Raw;
end


Y=y;




% removing rows with NaNs after regression?????
% Identify rows with any NaN
%rowsWithAnyNaN_Y = any(isnan(Y),2);
%rowsWithAnyNaN_x = any(isnan(x),2);
% Identify rows that are entirely NaN
%rowsWithAllNaN_Y = all(isnan(Y),2);
%rowsWithAllNaN_x = all(isnan(x),2);
% Combine - TRUE if row should be removed
%rowsToRemove = rowsWithAnyNaN_Y | rowsWithAnyNaN_x | rowsWithAllNaN_Y | rowsWithAllNaN_x;
% Keep only rows that are not in rowsToRemove
%naninx = ~rowsToRemove;
%naninx = ~any(isnan(Y),2) | ~any(isnan(x),2);
%Y=Y(naninx,:);
%x=x(naninx,:);
%sex=sex(naninx);
%site=site(naninx);
%motion=motion(naninx);
%icv=icv(naninx);
%i0age=i0age(naninx);
%i2age=i2age(naninx);
%delta_time=delta_time(naninx);



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
xlim([0 0.0019]); 
xline(0.0019, 'r', 'LineWidth', 2); 
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

PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');
PLS3sw=std(PLS3weights');

plsweights1=PLS1w./PLS1sw';
plsweights2=PLS2w./PLS2sw'; 
plsweights3=PLS3w./PLS3sw';

sum(plsweights1>3)
sum(plsweights1<-3)

sum(plsweights2>3)
sum(plsweights2<-3)

sum(plsweights3>3)
sum(plsweights3<-3)

protein_ids(plsweights1>3)'
protein_ids(plsweights1<-3)'

protein_ids(plsweights2>3)'
protein_ids(plsweights2<-3)'

protein_ids(plsweights3>3)'
protein_ids(plsweights3<-3)' 

protein_ids = protein_ids(:);
 cd /Users/juliannehu/Documents/UKB_proteomics/data/dictionaries
 nightingale = readtable('Nightingale_biomarker_groups.txt');
 proteins = regexprep(protein_ids, {'^p', '_i\d+$'}, '');
 [found, idx] = ismember(proteins, string(nightingale.field_id));
 met_names = strings(size(proteins));
 met_names(found) = string(nightingale.title(idx(found)));


PLS_table = table(met_names, plsweights1, plsweights2, plsweights3, ...
    'VariableNames', {'Metabolite', ...
                      'PLS1_weight', ...
                      'PLS2_weight', 'PLS3_weight'});
writetable(PLS_table, 'PLS_metabolite_weights_i0metabolomics_i2mri_with_age.csv');

%metabolomic_ids=merged_data.Properties.VariableNames(2:252);

%protein_ids = metabolomic_ids.title(:)
%PLS_table = table(protein_ids, plsweights1, plsweights2, plsweights3, ...
 %   'VariableNames', {'Protein', ...
  %                    'PLS1_weight', ...
   %                   'PLS2_weight', 'PLS3_weight'});
%writetable(PLS_table, 'PLS_protein_weights.csv');



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

% %FIGURE WITH NIGHTINGALE GROUPINGS
 cd /Users/juliannehu/Documents/UKB_proteomics/data/dictionaries
 nightingale = readtable('Nightingale_biomarker_groups.txt');
  
 
  %proteins = string(proteins)
 %met_title = string(nightingale.title)
 %met_group = string(nightingale.Subgroup)
%protein_ids = protein_ids(:);

% Define significance threshold
sig_threshold = (1 - normcdf(3));

% Find out the group for each metabolite
  proteins = regexprep(protein_ids, {'^p', '_i\d+$'}, '');
  [found, idx] = ismember(proteins, string(nightingale.field_id));
 met_group = strings(size(proteins));
 met_group(found) = string(nightingale.Subgroup(idx(found)));
 
% List of groups to merge together
lipids_groups = ["Other lipids","Total lipids"];
tri_phospho_groups = ["Triglycerides","Phospholipids"];
lipo_groups = ["Lipoprotein particle concentrations","Lipoprotein particle sizes","Apolipoproteins"];
chol_groups = ["Cholesterol", "Cholesteryl esters", "Free cholesterol"];
hdl_chol_groups = ["Very large HDL (average diameter 14.3 nm)", "Large HDL (average diameter 12.1 nm)", "Medium HDL (average diameter 10.9 nm)", "Small HDL (average diameter 8.7 nm)", "Very large HDL ratios", "Large HDL ratios", "Medium HDL ratios", "Small HDL ratios", "HDL"];
ldl_chol_groups = ["Chylomicrons and extremely large VLDL (particle diameters from 75 nm upwards)", "Very large VLDL (average diameter 64 nm)", "Large VLDL (average diameter 53.6 nm)", "Medium VLDL (average diameter 44.5 nm)", "Small VLDL (average diameter 36.8 nm)", "Very small VLDL (average diameter 31.3 nm)", "Large LDL (average diameter 25.5 nm)", "Medium LDL (average diameter 23 nm)", "Small LDL (average diameter 18.7 nm)", "Chylomicrons and extremely large VLDL ratios", "Very large VLDL ratios", "Large VLDL ratios", "Medium VLDL ratios", "Small VLDL ratios", "Very small VLDL ratios", "IDL ratios", "Large LDL ratios", "Medium LDL ratios", "Small LDL ratios", "IDL (average diameter 28.6 nm)", "VLDL", "LDL"];
fatty_groups = ["Fatty acids", "Fatty acid ratios"];
amino_groups = ["Amino acids", "Branched-chain amino acids", "Aromatic amino acids"];
other_groups = ["Glycolysis related metabolites","Ketone bodies","Fluid balance","Inflammation"];

is_lipids = ismember(met_group, lipids_groups);
is_tri_phospho = ismember(met_group, tri_phospho_groups);
is_lipo = ismember(met_group, lipo_groups);
is_chol = ismember(met_group, chol_groups);
is_hdl = ismember(met_group, hdl_chol_groups);
is_ldl = ismember(met_group, ldl_chol_groups);
is_fatty = ismember(met_group, fatty_groups);
is_amino = ismember(met_group, amino_groups);
is_other = ismember(met_group, other_groups);

% Metabolite categories and colors
metabolite_categories = {
    is_lipids,        [0.30 0.70 0.30]; % green shade
    is_tri_phospho,   [0.70 0.90 0.70]; % light green
    is_lipo,          [0.00 0.35 0.00]; % dark green
    is_chol,          [1.00 0.55 0.00]; % orange
    is_hdl,           [0.25 0.50 0.85]; % blue
    is_ldl,           [1.00 0.00 0.00]; % red
    is_fatty,         [0.95 0.85 0.25]; % yellow
    is_amino,         [0.85 0.10 0.85]; % magenta
    is_other,         [0.85 0.85 0.85]; % gray
};

% VOLCANO PLOT FOR PLS 1
p_pls1 = 2 * (1 - normcdf(abs(plsweights1)));
figure; hold on;
for i = 1:size(metabolite_categories,1)
    idx = metabolite_categories{i,1};
    color = metabolite_categories{i,2};
    
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
%xlabel('PLS Weights (I0-I2 Component 1 Age Regressed)');
ylabel('-log_{10}(p)');


% VOLCANO PLOT FOR PLS 2
p_pls2 = 2 * (1 - normcdf(abs(plsweights2)));
figure; hold on;
for i = 1:size(metabolite_categories,1)
    idx = metabolite_categories{i,1};
    color = metabolite_categories{i,2};
    
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
%xlabel('PLS Weights (I0-I2 Component 1 Age Regressed)');
ylabel('-log_{10}(p)');


% VOLCANO PLOT FOR PLS 3
p_pls3 = 2 * (1 - normcdf(abs(plsweights3)));
figure; hold on;
for i = 1:size(metabolite_categories,1)
    idx = metabolite_categories{i,1};
    color = metabolite_categories{i,2};
    
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
%xlabel('PLS Weights (I0-I2 Component 1 Age Regressed)');
ylabel('-log_{10}(p)');



h = scatter(stats.W(:,2), -log10(p_pls2), 10);
h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Protein', protein_ids);


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


%SCATTER PLOT WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(age_i2-min(age_i2))./max(age_i2-min(age_i2));
colormap_custom(:,3)=1-age_i2./max(age_i2);
figure;scatter(XS(:,ix1), YS(:,ix2), 5, colormap_custom, 'filled'); 



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



%SCATTER PLOT FOR PLS 1 XS AND YS WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(i2age-min(i2age))./max(i2age-min(i2age));
colormap_custom(:,3)=1-i2age./max(i2age);
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
c.Label.String = 'Age (years)';
c.Ticks = linspace(0,1,5);
c.TickLabels = round(linspace(min(i2age), max(i2age), 5));



%SCATTER PLOT FOR LDL cholesterol AND YS1 WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(i2age-min(i2age))./max(i2age-min(i2age));
colormap_custom(:,3)=1-i2age./max(i2age);
i = find(strcmp(met_names, 'LDL cholesterol'))
x_val = x(:,i)
y = YS(:,1)
figure;
scatter(x_val, y, 5, colormap_custom, 'filled'); 
lsline;
corr(x_val, y)


%SCATTER PLOT FOR HDL cholesterol AND YS1 WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(i2age-min(i2age))./max(i2age-min(i2age));
colormap_custom(:,3)=1-i2age./max(i2age);
i = find(strcmp(met_names, 'HDL cholesterol'))
x_val = x(:,i)
y = YS(:,1)
figure;
scatter(x_val, y, 5, colormap_custom, 'filled'); 
lsline;
corr(x_val, y)


%SCATTER PLOT FOR LDL cholesterol AND thickness of insula WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(i2age-min(i2age))./max(i2age-min(i2age));
colormap_custom(:,3)=1-i2age./max(i2age);
i = find(strcmp(met_names, 'HDL cholesterol'))
x_val = x(:,i)
j = find(strcmp(mri, 'MeanThicknessOfInsula_rightHemisphere__Instance2'))
y = Y(:,j)
figure;
scatter(x_val, y, 5, colormap_custom, 'filled'); 
lsline;
corr(x_val, y)



%SCATTER PLOT FOR Glycoprotein acetyls AND YS2 WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(i2age-min(i2age))./max(i2age-min(i2age));
colormap_custom(:,3)=1-i2age./max(i2age);
i = find(strcmp(met_names, 'Glycoprotein acetyls'))
x_val = x(:,i)
y = YS(:,2)
figure;
scatter(x_val, y, 5, colormap_custom, 'filled'); 
lsline;
corr(x_val, y)


%SCATTER PLOT FOR Glucose AND YS2 WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(i2age-min(i2age))./max(i2age-min(i2age));
colormap_custom(:,3)=1-i2age./max(i2age);
i = find(strcmp(met_names, 'Glucose'))
x_val = x(:,i)
y = YS(:,2)
figure;
scatter(x_val, y, 5, colormap_custom, 'filled'); 
lsline;
corr(x_val, y)


%SCATTER PLOT FOR Ratio of linoleic acid to total fatty acids AND thickness of INSULA WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(i2age-min(i2age))./max(i2age-min(i2age));
colormap_custom(:,3)=1-i2age./max(i2age);
i = find(strcmp(met_names, 'Ratio of linoleic acid to total fatty acids'))
x_val = x(:,i)
j = find(strcmp(mri, 'MeanThicknessOfInsula_rightHemisphere__Instance2'))
y = Y(:,j)
figure;
scatter(x_val, y, 5, colormap_custom, 'filled'); 
lsline;
corr(x_val, y)


%correlation matrix
r=corr(X,Y);
figure; imagesc(r'); colorbar
cd /Users/juliannehu/Documents/UKB_proteomics/results/met-mri
correlation_matrix = array2table(r', ...
    'VariableNames', met_names, ...
    'RowNames', mri_variables);
writetable(correlation_matrix, 'ukb_i0metabolomics_i2mri_correlation_matrix.csv', 'WriteRowNames', true);


% correlation matrix with clustering
% 1. Compute the correlation matrix
corr_map = corr(X,Y);
% 2. Calculate distances for x (1 - correlation) and get optimal leaf order for visualization
distX = pdist(corr_map,'correlation');
treeX = linkage(distX,'average');
leafOrderX = optimalleaforder(treeX,distX);
% 2. Calculate distances for y (1 - correlation) and get optimal leaf order for visualization
distY = pdist(corr_map','correlation');
treeY = linkage(distY,'average');
leafOrderY = optimalleaforder(treeY,distY);
% 5. Reorder the correlation matrix using the leaf order
sorted_corr_map = corr_map(leafOrderX, leafOrderY);
% 6. Reorder variables names using the leaf order
sorted_met_names = met_names(leafOrderX);
sorted_mri_variables = mri_variables(leafOrderY);
% 7. Visualize
figure; imagesc(sorted_corr_map')
colorbar
% saving correlations as csv
cd /Users/juliannehu/Documents/UKB_proteomics/results/met-mri
correlation_matrix = array2table(sorted_corr_map', ...
    'VariableNames', sorted_met_names, ...
    'RowNames', sorted_mri_variables);
writetable(correlation_matrix, 'clustered_ukb_i0metabolomics_i2mri_correlation_matrix.csv', 'WriteRowNames', true);




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
    writetable(mri_results, sprintf('mri_results_with_age_i0i2_pls%d.csv', i));
end

%%
i = find(strcmp(met_names, 'Glucose'))
x_val = x(:,i);
i = find(strcmp(mri_variables, 'MeanThicknessOfPrecuneus_leftHemisphere__Instance2'))
y = Y(:,i);
figure;
dscatter(x_val, y);
lsline;
corr(x_val, y)

figure;
dscatter(x_val, YS(:,1));
lsline;
corr(x_val, YS(:,1))




%% MEDIATION MODEL %%
X_mediation = X; 
%X_mediation = XS(:1);
Y_mediation = YS(:,1);
age_i2=merged_data.p21003_i2;
age_i2=age_i2(naninx==0);
M_mediation = age_i2;

[paths, stats] = mediation(X_mediation, Y_mediation, M_mediation, 'plots', 'verbose', 'boot', 'bootreps', 5000);

[r,p]=corr(X,Y)
%74x250 after covariates 