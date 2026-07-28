clear

% Uploading files
cd /Users/arshiya/Desktop/CLSA/2601007_CAMH_PZhukovsky_BL/
cognitive_baseline=readtable('2601007_CAMH_PZhukovsky_Baseline_CoPv7-2.csv');

cd /Users/arshiya/Desktop/CLSA/2601007_CAMH_PZhukovsky_FUP2
cognitive=readtable('2601007_CAMH_PZhukovsky_FUP2_CoPv4-2.csv');

cognitive=innerjoin(cognitive_baseline, cognitive);
cognitive=cognitive_baseline;

cd /Users/arshiya/Desktop/CLSA/'2.CLSA METABOLOMICS_v2_July2024'/
metabolomics=readtable('CLSA_QC_NORM_COMMON_META_NO_XENOBIOTICS_TRANSPOSED_RUN2.csv');

metabolomics = renamevars(metabolomics, 'Var1', 'ADM_METABOLON2_COM');

% Displaying ChemIDs as the headers
opts = detectImportOptions('CLSA_QC_NORM_COMMON_META_NO_XENOBIOTICS_TRANSPOSED_RUN2.csv');
opts.VariableNamesLine = 1;
opts.DataLines = [2, Inf];
metabolomics = readtable('CLSA_QC_NORM_COMMON_META_NO_XENOBIOTICS_TRANSPOSED_RUN2.csv', opts);

% Final merged dataset 
merged_data=innerjoin(metabolomics,cognitive, 'Keys','ADM_METABOLON2_COM');

%% Covariates
i0age=merged_data.AGE_NMBR_COM; 
%i2age=merged_data.AGE_NMBR_COF2; 
sex=categorical(merged_data.SEX_ASK_COM); 
site=categorical(merged_data.GEOSTRATA_COM);

% defining metabolomics
metabolomic_ids=merged_data.Properties.VariableNames(2:558);

%%% x variable: metabolites
X = merged_data{:, 2:558};
Y=[merged_data.COG_REYII_SCORE_COM, merged_data.COG_MAT_SCORE_COM];

% Removing NaNs
valid_Y = all(~isnan(Y),2);valid_sex = ~isundefined(sex);valid_site = ~isundefined(site);
naninx = valid_Y & valid_sex & valid_site;
Y = Y(naninx,:); X = X(naninx,:);
i0age = i0age(naninx); %i2age = i2age(naninx);
sex = sex(naninx); site = site(naninx);

Y=zscore(Y); X=zscore(X);%  Y(:,1)=Y(:,1)*-1;
ncomp=3

mdl = fitlm(table(i0age, sex, site,Y(:,1)) ); Y(:,1)=mdl.Residuals.Raw;
mdl = fitlm(table(i0age, sex, site, Y(:,2))); Y(:,2) = mdl.Residuals.Raw;
clear x; for i=1:length(X(1,:)); mdl = fitlm(table(i0age, sex, site,X(:,i)) ); x(:,i)=mdl.Residuals.Raw;end
%x=X;

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);PCTVAR

%% Plotting each of the components' correlation with cognitive tests
figure(1);imagesc(corr(XS(:,1:3),Y)); colormap bone ; colorbar; corr(XS, Y) 
combs=allcomb([1:3], [1:4]); figure(13); clf; 
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

figure;scatter(XS(:,1), YS(:,1));lsline
%% Permutation Testing
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
figure(2); histogram(Rsq); xlim([0 0.02]); xline(0.0305, 'r'); sum(PCTVAR')
%% Bootstrapping to get the func connectivity weights for PLS1, 2 and 3

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

% Filtering bootstrap results to identify variables that reliably
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

% Flip PLS1 and PLS2 weights 
plsweights1= plsweights1*-1;
stats.W(:,1) = -stats.W(:,1);
plsweights2 = plsweights2*-1;
stats.W(:,2) = -stats.W(:,2);

protein_ids = metabolomic_ids(:)
PLS_table = table(protein_ids, plsweights1, plsweights2, plsweights3, ...
    'VariableNames', {'CHEM_ID', ...
                      'PLS1_weight', ...
                      'PLS2_weight', 'PLS3_weight'});
writetable(PLS_table, 'CLSA_PLS_weights_5000.csv');


%% Volcano plots
anno=readtable('ANNOTATION TABLE_V2.CSV');
anno_ids=string(anno.CHEM_ID);
anno_pathway = string(anno.SUPER_PATHWAY);
anno_lipid_pathway = string(anno.SUB_PATHWAY);

% Define significance threshold
sig_threshold = (1 - normcdf(3));

% Find out the group for each metabolite
met_titles = string(metabolomic_ids);
met_titles = erase(string(metabolomic_ids(:)), 'x');
[found, idx] = ismember(met_titles,anno_ids);
protein_group = strings(size(met_titles));
protein_group(found) = anno_pathway(idx(found));

lipid_subgroup = strings(size(met_titles));
lipid_subgroup(found) = anno_lipid_pathway(idx(found));
lipid_subgroup(~ismember(protein_group,"Lipid"))="";

% Lipid groups to merge together
other_lipids = {'Androgenic Steroids', 'Carnitine Metabolism', 'Ceramides', ...
    'Corticosteroids', 'Diacylglycerol', 'Dihydroceramides', 'Endocannabinoid', ...
    'Glycerolipid Metabolism', 'Hexosylceramides (HCER)', 'Inositol Metabolism', ...
    'Ketone Bodies', 'Lactosylceramides (LCER)', 'Lysophospholipid', 'Lysoplasmalogen', ...
    'Mevalonate Metabolism', 'Monoacylglycerol', 'Phosphatidylcholine (PC)', ...
    'Phosphatidylethanolamine (PE)', 'Phosphatidylinositol (PI)', 'Phospholipid Metabolism', ...
    'Plasmalogen', 'Pregnenolone Steroids', 'Primary Bile Acid Metabolism', ...
    'Progestin Steroids', 'Secondary Bile Acid Metabolism', 'Sphingolipid Synthesis', ...
    'Sphingosines', 'Sterol'};
fatty_acids = {'Fatty Acid Metabolism (Acyl Carnitine, Dicarboxylate)', 'Fatty Acid Metabolism (Acyl Carnitine, Hydroxy)', ...
    'Fatty Acid Metabolism (Acyl Carnitine, Long Chain Saturated)', 'Fatty Acid Metabolism (Acyl Carnitine, Medium Chain)', ...
    'Fatty Acid Metabolism (Acyl Carnitine, Monounsaturated)', 'Fatty Acid Metabolism (Acyl Carnitine, Polyunsaturated)', ...
    'Fatty Acid Metabolism (Acyl Carnitine, Short Chain)', 'Fatty Acid Metabolism (Acyl Glutamine)', 'Fatty Acid Metabolism (Acyl Glycine)', ...
    'Fatty Acid Metabolism (also BCAA Metabolism)', 'Fatty Acid, Amino', 'Fatty Acid, Branched', 'Fatty Acid, Dicarboxylate', ...
    'Fatty Acid, Dihydroxy', 'Fatty Acid, Monohydroxy', 'Long Chain Monounsaturated Fatty Acid', ...
    'Long Chain Polyunsaturated Fatty Acid (n3 and n6)', 'Long Chain Saturated Fatty Acid', ...
    'Medium Chain Fatty Acid', 'Short Chain Fatty Acid'};
sphingo = {'Sphingomyelins', 'Dihydrosphingomyelins'};

is_lipid = ismember(lipid_subgroup, other_lipids);
is_fatty_acids = ismember(lipid_subgroup, fatty_acids);
is_sphingo = ismember(lipid_subgroup, sphingo);
is_peptide = ismember(protein_group, "Peptide");
is_PCM = ismember(protein_group, "Partially Characterized Molecules");
is_nucleotide = ismember(protein_group, "Nucleotide");
is_energy = ismember(protein_group, "Energy");
is_cofacandvit = ismember(protein_group, "Cofactors and Vitamins");
is_carb = ismember(protein_group, "Carbohydrate");
is_amino = ismember(protein_group, "Amino Acid");

% Metabolite categories and colors
metabolite_categories = {
    is_lipid,        [0.30 0.70 0.30]; % green shade
    is_fatty_acids,  [0.70 0.90 0.70]; % light green
    is_sphingo,      [0.00 0.35 0.00]; % dark green
    is_nucleotide,   [1.00 0.55 0.00]; % orange
    is_energy,       [0.25 0.50 0.85]; % blue
    is_cofacandvit,  [1.00 0.00 0.00]; % red
    is_carb,         [0.95 0.85 0.25]; % yellow
    is_amino,        [0.85 0.10 0.85]; % magenta
    is_peptide,      [0.50 0.00 0.50]; % purple
    is_PCM,          [1.00 0.70 0.80]; % light pink
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

%% VOLCANO PLOTS FOR PLS 3 WITH TRANSPARENCY FOR NON SIG
 p_pls1 = 2 * (1 - normcdf(abs(plsweights3)));
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

%% SCATTER PLOT FOR PLS 1 XS and YS 
ix1=1; ix2=1
mdl= fitlm(XS(:, ix1), YS(:,ix2)); 
[ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
[sorted_x, sortIdx] = sort(XS(:, ix1)); ypred_sorted = ypred(sortIdx); yci_sorted = yci(sortIdx, :);
figure;
dscatter(XS(:,ix1), YS(:,ix2));
hold on
plot(sorted_x, ypred_sorted, 'k', 'LineWidth', 2);
plot(sorted_x, yci_sorted(:,1), 'b', 'LineWidth', 0.5);
plot(sorted_x, yci_sorted(:,2), 'b', 'LineWidth', 0.5);
%xlabel(sprintf('X Scores (LV%d)', ix1), 'FontSize', 12);
%ylabel(sprintf('Y Scores (LV%d)', ix2), 'FontSize', 12);
corr(XS(:,ix1), YS(:,ix2))


%% Example Scatter Plots

% SCATTER PLOT FOR DHA AND REYII
i = find(strcmp(metabolomic_ids, 'x266'))
x_val = x(:,i)
y = Y(:,1)
figure;
dscatter(x_val, y);
lsline;
corr(x_val, y)

% SCATTER PLOT FOR DHA and MAT
i = find(strcmp(metabolomic_ids, 'x100000665'))
x_val = x(:,i)
y = Y(:,2)
figure;
scatter(x_val, y, 10, 'filled');
lsline;
corr(x_val, y)

% SCATTER PLOT FOR cholesterol AND REYII
i = find(strcmp(metabolomic_ids, 'x565'))
x_val = x(:,i)
y = Y(:,1)
figure;
dscatter(x_val, y);
lsline;
corr(x_val, y)

% SCATTER PLOT FOR cholesterol AND MAT
i = find(strcmp(metabolomic_ids, 'x565'))
x_val = x(:,i)
y = Y(:,2)
figure;
dscatter(x_val, y);
lsline;
corr(x_val, y)
