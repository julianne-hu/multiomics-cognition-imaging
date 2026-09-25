clear
cd /Users/juliannehu/Documents/UKB_proteomics/data/raw
cognitive=readtable('cognitive_v2.csv');
demo=readtable('demographics_v2.csv');
cognitive=innerjoin(cognitive, demo);

%trail making
cognitive.p6350_i2(cognitive.p6350_i2> (nanmean(cognitive.p6350_i2)+4*nanstd(cognitive.p6350_i2)) | cognitive.p6350_i2<100)=NaN; %Duration to complete alphanumeric path 
%cognitive.p6348_i2(cognitive.p6348_i2<100 | cognitive.p6348_i2> (nanmean(cognitive.p6348_i2)+3*nanstd(cognitive.p6348_i2)) )=NaN; %Duration to complete numeric/easy path 
cognitive.tmt_cor=(cognitive.p6350_i2 + 50*cognitive.p6351_i2);% - (cognitive.x6348_2_0 +5*cognitive.x6349_2_0) ;


cd /Users/juliannehu/Documents/UKB_proteomics/data/cleaned
metabolomics=readtable('mean_imputed_metabolomic_i0_v2.csv');


merged_data=innerjoin(metabolomics,cognitive);
clear cognitive metabolomics demo;
i0age=merged_data.p21003_i0; i2age=merged_data.p21003_i2; sex=categorical(merged_data.p31); site=categorical(merged_data.p54_i2);

%%% requires tabular_ukb to have ran
%ix=(cellfun('isempty', clinical)' | sum(isnan(ms_ordered)')'>0 ) ; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
%ix=(cellfun('isempty', clinical)'| isnan(ica_partial(:,1))) ; a=age; a(ix)=[];s=sex;s(ix)=[];clin=clinical; clin(ix)=[];
%cognitive=cognitive; cognitive(ix,:)=[]; minimal=minimal_ordered; minimal(ix,:)=[]; MDD_prs=MDD_prs_ordered; MDD_prs(ix)=[];

%%%%%%IMPUTATION --- using mean imputation instead%%%%%%%
%X=ms_ordered; X(ix,:)=[]; 
%X=merged_data{:,2:2924}; 
metabolomic_ids=merged_data.Properties.VariableNames(2:252);
%m=nanmean(X);
%for i=1:length(protein_ids)
%    X(isnan(X(:,i)),i)=m(i);
%end


%%% x variable: proteomics
X = merged_data{:, 2:252};
%%% y variables: alphanumtrails, fluid intelligence,  word pairs (PAL), symbol digit, %%%% pairs match tower  cognitive.x21004_2_0,
Y=[merged_data.tmt_cor, merged_data.p20016_i2, merged_data.p20197_i2, merged_data.p23324_i2];%   Y=[cognitive.tmt_cor, cognitive.x21004_2_0, cognitive.x20197_2_0]; % cognitive.x399_2_2+cognitive.x399_2_1, 
%Y=age;
%Y=[merged_data.tmt_cor, merged_data.p20197_i2, merged_data.p23324_i2];%   Y=[cognitive.tmt_cor, cognitive.x21004_2_0, cognitive.x20197_2_0]; % cognitive.x399_2_2+cognitive.x399_2_1, 
%%% no controls? remove them here
% index=~(  ismember(clin','ahc')   ); % %
% index=~ismember(clin','ahc') & ~ismember(clin', 'depanx') & ~ismember(clin', 'ahc') & ~ismember(clin', 'str');
%Y=Y(index,:); X=X(index,:); minimal=minimal(index,:); c=clin(index); MDD_prs=MDD_prs(index);
%%% remove nans

valid_Y = all(~isnan(Y),2);
valid_sex = ~isundefined(sex);
valid_site = ~isundefined(site);
naninx= valid_Y & valid_sex & valid_site; 
Y=Y(naninx,:); 
X=X(naninx,:); 
i0age=i0age(naninx); 
i2age=i2age(naninx); 
sex=sex(naninx);
site=site(naninx); 

%Y=Y(naninx,:); X=X(naninx,:); sex=sex(naninx); site=site(naninx);
Y=zscore(Y); X=zscore(X); Y(:,1)=Y(:,1)*-1;
ncomp=3

%%% regress out sex and age site and potentially other covariates
mdl = fitlm(table(i2age,sex, site,Y(:,1)) ); Y(:,1)=mdl.Residuals.Raw;
mdl = fitlm(table(i2age,sex, site,Y(:,2)) ); Y(:,2)=mdl.Residuals.Raw;
mdl = fitlm(table(i2age,sex, site,Y(:,3)) ); Y(:,3)=mdl.Residuals.Raw;
mdl = fitlm(table(i2age,sex, site,Y(:,4)) ); Y(:,4)=mdl.Residuals.Raw;
clear x; 
for i=1:length(X(1,:)); 
    mdl = fitlm(table(i0age,sex, site,X(:,i)) ); 
    x(:,i)=mdl.Residuals.Raw;
end
%x=X;

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);PCTVAR
%%% explore the PLS components and their correlation with Y
%XS VS YS FIGURE:%               figure(31); scatter(XS(ismember(c','dep'),1), YS(ismember(c','dep'),1), 5,'b', 'filled'); hold on;scatter(XS(ismember(c','depanx'),1), YS(ismember(c','depanx'),1), 5, 'r', 'filled'); hold on;scatter(XS(ismember(c','anx'),1), YS(ismember(c','anx'),1), 5, 'k', 'filled'); mdl= fitlm(XS(:, 1), YS(:,1)); [ypred,yci] = predict(mdl,XS(:, 1), 'Alpha',0.001); hold on; plot(XS(:, 1), ypred, 'k', 'LineWidth', 2);plot(XS(:, 1), yci(:,1), 'k', 'LineWidth', 0.5);plot(XS(:, 1), yci(:,2), 'k', 'LineWidth', 0.5);
% and their distributions: %        figure; plot_histogram_shaded(XS,'Alpha',0.3,'color',[0.5 0.5 0.5], 'Normalization', 'pdf');figure; plot_histogram_shaded(YS,'Alpha',0.3,'color',[0.5 0.5 0.5], 'Normalization', 'pdf');
%% plotting each of the components' correlation with cognitive tests
figure(1);imagesc(corr(XS(:,1:3),Y)); colormap bone ; colorbar; corr(XS, Y) %1.TMT  2.Tower  3.PAL  4.DSST
combs=allcomb([1:3], [1:4]); figure(13); clf; % Clear the figure to ensure it is ready for new plots
for i=1:length(combs)
hold off; ix1=combs(i,1);ix2=combs(i,2); subplot(3,4,i); 
scatter(XS(:,ix1), Y(:,ix2), 5,'b', 'filled'); 
hold on;scatter(XS(:,ix1), Y(:,ix2), 5, 'r', 'filled')
hold on;scatter(XS(:,ix1), Y(:,ix2), 5, 'k', 'filled');
%mdl= fitlm(XS(:, ix1), Y(:,ix2)); [ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
%plot(XS(:, ix1), ypred, 'k', 'LineWidth', 2);
%plot(XS(:, ix1), yci(:,1), 'k', 'LineWidth', 0.5);
%plot(XS(:, ix1), yci(:,2), 'k', 'LineWidth', 0.5);
end; clear mdl ypred yci mdl ix1 ix2 combs

corr(XS, Y)
max(corr(X, Y))
min(corr(X, Y))

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
figure(2); histogram(Rsq); xlim([0 0.004]); xline(0.0038, 'r'); sum(PCTVAR')
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

PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');
PLS3sw=std(PLS3weights');

plsweights1=PLS1w./PLS1sw';
plsweights2=PLS2w./PLS2sw'; 
plsweights3=PLS3w./PLS3sw';

% filtering bootstrap results to identify variables that reliably
% contribute to each component
sum(plsweights1 > 3)
sum(plsweights1 < -3)

sum(plsweights2 > 3)
sum(plsweights2 < -3)

sum(plsweights3 > 3)
sum(plsweights3 < -3)

metabolomic_ids(plsweights1 > 3)'
metabolomic_ids(plsweights1 < -3)'

metabolomic_ids(plsweights2 > 3)'
metabolomic_ids(plsweights2 < -3)'

metabolomic_ids(plsweights3 > 3)'
metabolomic_ids(plsweights3 < -3)'


protein_ids = metabolic_ids.title(:)
PLS_table = table(protein_ids, plsweights1, plsweights2, plsweights3, ...
    'VariableNames', {'Protein', ...
                      'PLS1_weight', ...
                      'PLS2_weight', 'PLS3_weight'});
writetable(PLS_table, 'PLS_protein_weights.csv');


% %FIGURE WITH NIGHTINGALE GROUPINGS
 cd /Users/juliannehu/Documents/UKB_proteomics/data/dictionaries
 nightingale = readtable('Nightingale_biomarker_groups.txt');
 met_title = string(nightingale.title)
 met_group = string(nightingale.Subgroup)

% Define significance threshold
sig_threshold = (1 - normcdf(3));

% Find out the group for each metabolite
protein_titles = (protein_ids);  
[found, idx] = ismember(protein_titles, met_title);  
protein_group = strings(size(protein_titles));
protein_group(found) = met_group(idx(found));

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

is_lipids = ismember(protein_group, lipids_groups);
is_tri_phospho = ismember(protein_group, tri_phospho_groups);
is_lipo = ismember(protein_group, lipo_groups);
is_chol = ismember(protein_group, chol_groups);
is_hdl = ismember(protein_group, hdl_chol_groups);
is_ldl = ismember(protein_group, ldl_chol_groups);
is_fatty = ismember(protein_group, fatty_groups);
is_amino = ismember(protein_group, amino_groups);
is_other = ismember(protein_group, other_groups);

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

h = scatter(stats.W(:,2), -log10(p_pls2), 10);
h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Protein', protein_ids);



%p_pls2 = 2 * (1 - normcdf(abs(plsweights2)));

%figure;scatter(stats.W(:,2),-log10(p_pls2),'filled','k','SizeData', 10);
%hold on; ix=abs(plsweights2)>3;
%scatter(stats.W(ix,2),-log10(p_pls2(ix)),'b','SizeData', 10);

%yline(-log10(2*(1 - normcdf(3))), 'r--', 'LineWidth', 1.2); 

%xlabel('PLS Weights (I0-I2 Component 2 Age Regressed)');
%ylabel('-log_{10}(p)');






%VOLCANO PLOTS FOR PLS 3 WITH TRANSPARENCY FOR NON SIG
 p_pls1 = 2 * (1 - normcdf(abs(plsweights1)));
% Define significance threshold
sig_threshold = 2*(1 - normcdf(3));
% Set colors
color_sig = [0.2 0.6 0.8];  
color_nonsig = [0.2 0.6 0.8];
sig_idx = p_pls1 <= sig_threshold;
nonsig_idx = p_pls1 > sig_threshold;

figure; hold on;

% Significant points
scatter(stats.W(sig_idx,3), -log10(p_pls1(sig_idx)), 25, ...
    'MarkerFaceColor', color_sig, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);

% Non-significant points
scatter(stats.W(nonsig_idx,3), -log10(p_pls1(nonsig_idx)), 25, ...
    'MarkerFaceColor', color_nonsig, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);

yline(-log10(sig_threshold), 'k--', 'LineWidth', 1.2); 
%xlabel('PLS Weights (I0-I2 PLS 1 Age Regressed)');
ylabel('-log_{10}(p)');







figure;scatter(stats.W(:,3),-log10(p_pls1),'filled','k');
hold on; ix=abs(plsweights3)>3;
scatter(stats.W(ix,3),-log10(p_pls1(ix)),'b');

yline(-log10(2*(1 - normcdf(3))), 'r--', 'LineWidth', 1.2); 

xlabel('PLS Weights (I0-I2 Component 3 Age Regressed)');
ylabel('-log_{10}(p)');








%SCATTER PLOT FOR PLS 1 XS and YS 
figure;
ix1=1; ix2=1
dscatter(XS(:,ix1), YS(:,ix2)); 
mdl= fitlm(XS(:, ix1), YS(:,ix2)); [ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
plot(XS(:, ix1), ypred, 'k', 'LineWidth', 2);
plot(XS(:, ix1), yci(:,1), 'b', 'LineWidth', 0.5);
plot(XS(:, ix1), yci(:,2), 'b', 'LineWidth', 0.5);
xlabel(sprintf('X Scores (LV%d)', ix1), 'FontSize', 12);
ylabel(sprintf('Y Scores (LV%d)', ix2), 'FontSize', 12);
corr(XS(:,ix1), YS(:,ix2))


%%% y variables: alphanumtrails, fluid intelligence,  word pairs (PAL), symbol digit
%SCATTER PLOT FOR glycoprotein acetyls AND TMT
i = find(strcmp(protein_names, 'Glycoprotein acetyls'))
x_val = x(:,i)
y = Y(:,1)
figure;
dscatter(x_val, y);
lsline;
corr(x_val, y)

%SCATTER PLOT FOR glycoprotein acetyls AND TMT
i = find(strcmp(protein_names, 'Cholesterol in large HDL'))
x_val = x(:,i)
y = Y(:,1)
figure;
dscatter(x_val, y);
lsline;
corr(x_val, y)

%SCATTER PLOT FOR Docosahexaenoic acid AND TMT
i = find(strcmp(protein_names, 'Docosahexaenoic acid'))
x_val = x(:,i)
y = Y(:,1)
figure;
dscatter(x_val, y);
lsline;
corr(x_val, y)

%SCATTER PLOT FOR valine AND TMT
i = find(strcmp(protein_names, 'Valine'))
x_val = x(:,i)
y = Y(:,1)
figure;
dscatter(x_val, y);
lsline;
corr(x_val, y)

%SCATTER PLOT FOR valine AND SDS
i = find(strcmp(protein_names, 'Glucose'))
x_val = x(:,i)
y = Y(:,4)
figure;
dscatter(x_val, y);
lsline;
corr(x_val, y)

%SCATTER PLOT FOR Polyunsaturated fatty acid AND SDS
i = find(strcmp(protein_names, 'Polyunsaturated fatty acids'))
x_val = x(:,i)
y = Y(:,4)
figure;
dscatter(x_val, y);
lsline;
corr(x_val, y)
