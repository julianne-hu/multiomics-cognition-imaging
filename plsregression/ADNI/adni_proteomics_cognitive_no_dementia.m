clear
% Importing files
cognitive=readtable('NEUROBAT.csv');
diagnosis=readtable('DXSUM.csv');
demo=readtable('ADNI_Participant_Demographics.csv');
proteomics=readtable('mean_imputed_somascan_proteomics');

% Extracting diagnosis with dates
diagnosis = diagnosis(:, {'RID', 'DIAGNOSIS', 'EXAMDATE'});

% Calculating trail making test scores
cognitive.TRABSCOR(cognitive.TRABSCOR>(nanmean(cognitive.TRABSCOR)+4*nanstd(cognitive.TRABSCOR)) | cognitive.TRABSCOR<10)=NaN; % remove outliers
cognitive.tmt_cor=(cognitive.TRABSCOR+5*(cognitive.TRABERRCOM+cognitive.TRABERROM)); % add penalty for errors

%Age
% Convert VISDATE to datetime
demo.VISDATE = datetime(demo.VISDATE, 'InputFormat', 'yyyy-MM-dd');
% Convert PTDOB (MM/yyyy) to datetime
demo.PTDOB = datetime(demo.PTDOB, 'InputFormat', 'MM/yyyy');
% Calculate age at visit
demo.AgeAtVisit = years(demo.VISDATE - demo.PTDOB);
%demo.AgeAtVisit = floor(demo.AgeAtVisit);

% Joining files
%pro_demo = innerjoin(proteomics, demo, 'Keys', {'RID', 'VISCODE2'});
pro_demo = innerjoin(proteomics, demo, ...
                     'Keys', {'RID'});
%merged_data = innerjoin(pro_demo, cognitive, 'Keys', {'RID', 'VISCODE'});
merged_data = innerjoin(pro_demo, cognitive, ...
                         'Keys', {'RID'});
clear cognitive proteomics demo;

% Specifying date format
merged_data.EXAMDATE = datetime(merged_data.EXAMDATE, 'InputFormat', 'MM-dd-yyyy');
merged_data.VISDATE_cognitive = datetime(merged_data.VISDATE_cognitive, 'InputFormat', 'yyyy-MM-dd');
% Calculate change in time
merged_data.time_delta_pro_cog= abs(merged_data.EXAMDATE-merged_data.VISDATE_cognitive); 

% Only keeping correct timepoint rows
uniqueIDs = unique(merged_data.RID);
rows_to_keep = false(height(merged_data),1);
for i = 1:length(uniqueIDs)
    rid = uniqueIDs(i);
    % Find rows belonging to single participant
    idx = merged_data.RID == rid;
    % Find the minimum delta time in those
    [~, minIdx] = min(merged_data.time_delta_pro_cog(idx));
    rid_rows = find(idx);
    rows_to_keep(rid_rows(minIdx)) = true;
end
merged_data = merged_data(rows_to_keep, :);

% Convert diagnosis EXAMDATE to datetime
diagnosis.EXAMDATE = datetime(diagnosis.EXAMDATE, 'InputFormat', 'yyyy-MM-dd');
% Sort by RID and examdate
diagnosis = sortrows(diagnosis, {'RID', 'EXAMDATE'});

% Add in diagnosis column to merge file based on dates that match closest together 
diagnosis_cog = NaN(height(merged_data),1);

for i = 1:height(merged_data)
    rid = merged_data.RID(i);

    % Search all diagnosis for this RID
    rows = diagnosis.RID == rid;
    if any(rows)
        % Pulling diagnosis and exam date from diagnosis data
        diagnosis_dates = diagnosis.EXAMDATE(rows);
        diagnosis_labels = diagnosis.DIAGNOSIS(rows);

        % Pulling visit date from cognitive data
        cog_date = merged_data.VISDATE_cognitive(i);

        % Find diagnosis closest to cognitive visit
        [~, idx_x] = min(diagnosis_dates - cog_date);

        % Assigning closest label to participant
        diagnosis_cog(i) = diagnosis_labels(idx_x);
    end
end

% Add to merged data
merged_data.DIAGNOSIS_COG = diagnosis_cog;

% Checking before removal
fprintf('N before removal: %d\n', height(merged_data));

% Exclude participants with dementia 
cog_dx_idx = merged_data.DIAGNOSIS_COG ~=3;
merged_data = merged_data(cog_dx_idx,:);

% Check for removal
fprintf('N after dementia removal: %d\n', height(merged_data));

%Remove PTIDs with the format 381_S_##### as advised by ADNI
merged_data(contains(merged_data.PTID_cognitive, '381_S_'), :) = [];

% Variables
age=merged_data.AgeAtVisit; 
sex=merged_data.PTGENDER; 
site=merged_data.SITEID_pro_demo;
%Replace missing data for sex and site
sex(sex==-4) = NaN;
site(site==-4) = NaN;
% X = proteomics
X = merged_data{:, 8:7008}; 
protein_ids=merged_data.Properties.VariableNames(8:7008);
% Y: alphanumerical trails, RAVLT
Y=[merged_data.tmt_cor, merged_data.AVDEL30MIN];%   Y=[cognitive.tmt_cor, cognitive.x21004_2_0, cognitive.x20197_2_0]; % cognitive.x399_2_2+cognitive.x399_2_1, 

%%% Removing NaNs
%naninx=sum(isnan(Y)')'; 
%Y=Y(naninx==0,:); 
%X=X(naninx==0,:); 
%age=age(naninx==0); 
%sex=sex(naninx==0);
%site=site(naninx==0); 
rows_keep = ~any(isnan([X Y]), 2);
X=X(rows_keep, :);
Y=Y(rows_keep, :);
age=age(rows_keep);
sex=sex(rows_keep);
site=site(rows_keep);

% Convert to categorical
sex=categorical(sex);
site=categorical(site);

% Standardizing data
Y=zscore(Y); 
X=zscore(X); 

% Flip direction of trail making test (originally higher score = worse performance, now higher score = better performance)
Y(:,1)=Y(:,1)*-1;

% Set # of PLS components
ncomp=1;

% Regress out age, sex, and site (covariates)
mdl = fitlm(table(age, sex, site,Y(:,1)), 'CategoricalVars',[2 3]); Y(:,1)=mdl.Residuals.Raw;
mdl = fitlm(table(age, sex, site,Y(:,2)), 'CategoricalVars',[2 3]); Y(:,2)=mdl.Residuals.Raw;
%mdl = fitlm(table(age, sex, site,Y(:,3)) ); Y(:,3)=mdl.Residuals.Raw;
%mdl = fitlm(table(i2age, sex, site,Y(:,4)) ); Y(:,4)=mdl.Residuals.Raw;
clear x; 
for i=1:length(X(1,:)) 
    mdl = fitlm(table(age, sex, site,X(:,i)), 'CategoricalVars', [2 3]); 
    x(:,i)=mdl.Residuals.Raw;
end

%% Double checking for dementia patients
assert(~any(merged_data.DIAGNOSIS_COG == 3), 'Dementia patients still present — check diagnosis matching');

%% Run PLS
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);
PCTVAR

%% Plotting each of the components' correlation with cognitive tests
figure;
imagesc(corr(XS(:,1),Y)); 
colormap bone; 
colorbar; 
corr(XS, Y) %1.TMT  2.GF  3.PAL  4.DSST
combs=allcomb([1], [1:2]); 
figure; clf;
for i=1:length(combs)
    hold off; 
    ix1=combs(i,1);
    ix2=combs(i,2); 
    subplot(3,4,i); 
    scatter(XS(:,ix1), Y(:,ix2), 5, 'k', 'filled');
    %mdl= fitlm(XS(:, ix1), Y(:,ix2)); [ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
    %plot(XS(:, ix1), ypred, 'k', 'LineWidth', 2);
    %plot(XS(:, ix1), yci(:,1), 'k', 'LineWidth', 0.5);
    %plot(XS(:, ix1), yci(:,2), 'k', 'LineWidth', 0.5);
end; 
clear mdl ypred yci mdl ix1 ix2 combs

%% Calculate correlations
corr(age, XS)
corr(XS, Y)
max(corr(X, Y))
min(corr(X, Y))

%% permutation testing
permutations=5000;   
allobservations=Y; 
for ncomp=1
    
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


%% figure to show the permutation distribution vs actual distribution
figure; 
histogram(Rsq, 'FaceColor', [0.5 0.5 0.5]); 
xline(0.1484, 'r', 'LineWidth', 2); 
% Only show x and y axes
ax = gca;
ax.Box = 'off';          % removes top and right lines
ax.TickDir = 'out';      % ticks point outward
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
sum(PCTVAR')

%% bootstrapping to get the func connectivity weights for PLS1
dim=1
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,dim);PCTVAR
PLS1w=stats.W(:,1);

bootnum=5000;
PLS1weights=[];

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
end

PLS1sw=std(PLS1weights');

plsweights1=PLS1w./PLS1sw';

%flip sign of PLS 1, 2 and 3 weights to match upregulated vs downregulated
%plsweights1 = plsweights1*-1;
%stats.W(:,1) = -stats.W(:,1);

% filtering bootstrap results to identify variables that reliably contribute to each component
sum(plsweights1 > 3)
sum(plsweights1 < -3)

protein_ids(plsweights1 > 3)'
protein_ids(plsweights1 < -3)'

%% Adding Gene Names and exracting csv
somascan_dictionary=readtable("ADNI_CSF_SOMAscan7k_Dictionary");

% Format IDs to match dictionary
ids_formatted = strrep(protein_ids, '_', '.'); 
% Match to dictionary
[isFound, loc] = ismember(ids_formatted, somascan_dictionary.Analytes);
protein_names = protein_ids; 
protein_names(isFound) = somascan_dictionary.TargetFullName(loc(isFound));
protein_names=protein_names(:);


[isFound, loc] = ismember(ids_formatted, somascan_dictionary.Analytes);
gene_names = protein_ids; 
gene_names(isFound) = somascan_dictionary.EntrezGeneSymbol(loc(isFound));
gene_names=gene_names(:);

% Create the table
PLS_table = table( ...
    protein_names, ...
    plsweights1, ...
    'VariableNames', {'ProteinName', 'PLS1_weight'} ...
);
PLS_table.GeneName = gene_names;
writetable(PLS_table, 'ADNI_PLS_protein_weights.csv');

%% FIGURE WITH OLINK EXPLORE
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
gene_names_str = lower(string(gene_names));
% 
% % identify cardiometabolic proteins
 is_cardiometabolic = ismember(lower(gene_names_str), lower(cardiometabolic_names));
% % identify cardiometabolic II proteins
 is_cardiometabolic_2 = ismember(lower(gene_names_str), lower(cardiometabolic_2_names));
% % identify neurology proteins
 is_neurology = ismember(lower(gene_names_str), lower(neurology_names));
% % identify neurology II proteins
 is_neurology_2 = ismember(lower(gene_names_str), lower(neurology_2_names));
% % identify oncology proteins
 is_oncology = ismember(lower(gene_names_str), lower(oncology_names));
% % identify oncology II proteins
 is_oncology_2 = ismember(lower(gene_names_str), lower(oncology_2_names));
 % % identify inflammation proteins
 is_inflammation = ismember(lower(gene_names_str), lower(inflammation_names));
% % identify inflammation II proteins
 is_inflammation_2 = ismember(lower(gene_names_str), lower(inflammation_2_names));
%identify non overlapping proteins w/UKB
 no_overlap = ~( ...
    is_cardiometabolic | is_cardiometabolic_2 | ...
    is_neurology | is_neurology_2 | ...
    is_oncology | is_oncology_2 | ...
    is_inflammation | is_inflammation_2 );

%% VOLCANO PLOTS FOR PLS 1 WITH TRANSPARENCY FOR NON SIG
p_pls1 = 2 * (1 - normcdf(abs(plsweights1)));
% Define significance threshold
sig_threshold = (1 - normcdf(3));

% Protein categories and colors
protein_categories = {
    no_overlap,            [0.00 0.00 0.00];
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
    sig_idx = idx & (p_pls1 <= sig_threshold);
    scatter(stats.W(sig_idx,1), -log10(p_pls1(sig_idx)), 10, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);

    % Points below significance
    nonsig_idx = idx & (p_pls1 > sig_threshold);
    scatter(stats.W(nonsig_idx,1), -log10(p_pls1(nonsig_idx)), 10, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
end

yline(-log10((1 - normcdf(3))), 'k--', 'LineWidth', 1.2); 
ylabel('-log_{10}(p)');


%% FIND KEY PROTEINS FOR POST PROCESSING ON POWERPOINT
%PLS 1
sig_threshold = 2*(1 - normcdf(3));
is_sig = p_pls1 <= sig_threshold;

up_mask = stats.W(:,1) > 0 & is_sig;
down_mask = stats.W(:,1) < 0 & is_sig;

protein_category = repmat("None", length(gene_names_str), 1);

protein_category(is_oncology | is_oncology_2) = "Oncology";
protein_category(is_cardiometabolic | is_cardiometabolic_2) = "Cardiometabolic";
protein_category(is_neurology | is_neurology_2) = "Neurology";
protein_category(is_inflammation | is_inflammation_2) = "Inflammation";
protein_category(no_overlap) = "NaN";

%Create tables
up_tbl = table(protein_names(up_mask), ...
               gene_names_str(up_mask), ...
               protein_category(up_mask), ...
               stats.W(up_mask,1), ...
               p_pls1(up_mask), ...
               'VariableNames', {'Protein', 'Gene', 'Category', 'PLS_Weight','pValue'});

down_tbl = table(protein_names(down_mask), ...
                 gene_names_str(down_mask), ...
                 protein_category(down_mask), ...
                 stats.W(down_mask,1), ...
                 p_pls1(down_mask), ...
                 'VariableNames', {'Protein', 'Gene', 'Category', 'PLS_Weight','pValue'});
% sort by |PLS weight| (strongest effect at top)
up_tbl   = sortrows(up_tbl, 'PLS_Weight', 'descend');  
down_tbl = sortrows(down_tbl, 'PLS_Weight');           

% Display results
disp('Upregulated Proteins:');
disp(up_tbl);

disp('Downregulated Proteins:');
disp(down_tbl);

writetable(up_tbl, 'adni_csf_proteomics_cognitive_pls1_significant_upregulated_proteins_no_age.csv');
writetable(down_tbl, 'adni_csf_proteomics_cognitive_pls1_significant_downregulated_proteins_no_age.csv');


%% sort by |PLS weight| (strongest effect at top)
up_tbl   = sortrows(up_tbl, 'PLS_Weight', 'descend');  
down_tbl = sortrows(down_tbl, 'PLS_Weight');           

%% SCATTER PLOT FOR PLS 1 XS and YS 
figure;
ix1=1; ix2=1;
scatter(XS(:,ix1), YS(:,ix2), 5,'k', 'filled'); 
mdl= fitlm(XS(:, ix1), YS(:,ix2)); [ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
plot(XS(:, ix1), ypred, 'k', 'LineWidth', 2);
plot(XS(:, ix1), yci(:,1), 'b', 'LineWidth', 0.5);
plot(XS(:, ix1), yci(:,2), 'b', 'LineWidth', 0.5);
corr(XS(:,ix1), YS(:,ix2))
exportgraphics(gcf, 'scatter1.png','Resolution', 300);

%% Maybe figure out what to do with this block here - confirm no dementia patients remaining

all_rids = merged_data.RID;
rid_vector = all_rids(rows_keep);  % 240 rows
fprintf('The rows are %d\n', height(rid_vector));
subject_dx = merged_data.DIAGNOSIS_COG(rows_keep);  % 1=CN, 2=MCI

%% d
assert(~any(subject_dx == 3), 'Dementia patients still present — check diagnosis matching');
%% Define colors
colors = [0 0 1;      % CN = blue
          0.9 0.6 0.1]; % MCI = orange 

%% Scatter plot for a notch3 and TMT - might change these proteins to the more correlated ones from protein weights
i = find(strcmp(protein_ids, 'X8480_29'));
x_efemp1 = x(:,i);
y_val = Y(:,1);
figure; hold on;
for dx_val = 1:2
    idx_dx = subject_dx == dx_val;  
    scatter(x_efemp1(idx_dx), y_val(idx_dx), 20, colors(dx_val,:), 'filled');
end
%xlabel('Notch3');
%ylabel('TMT score');
%legend({'CN','MCI','Dementia'});
%regression line
coeffs = polyfit(x_efemp1, y_val, 1);  % linear fit
x_fit = linspace(min(x_efemp1), max(x_efemp1), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_efemp1, y_val)
exportgraphics(gcf, 'efemp1_tmt.png','Resolution',300);

%% Scatter plot for a notch3 and SDS
i = find(strcmp(protein_ids, 'X5108_72'));
x_notch3 = x(:,i);
y_val = Y(:,2);
figure; hold on;
for dx_val = 1:3
    idx_dx = subject_dx == dx_val;
    scatter(x_notch3(idx_dx), y_val(idx_dx), 20, colors(dx_val,:), 'filled');
end
%xlabel('Notch3');
%ylabel('SDS score');
%legend({'CN','MCI','Dementia'});
%regression line
coeffs = polyfit(x_notch3, y_val, 1);  % linear fit
x_fit = linspace(min(x_notch3), max(x_notch3), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_notch3, y_val)


%% Scatter plot for a notch3 and TMT
i = find(strcmp(protein_ids, 'X5108_72'));
x_notch3 = x(:,i);
y_val = Y(:,1);
figure; hold on;
for dx_val = 1:3
    idx_dx = subject_dx == dx_val;  
    scatter(x_notch3(idx_dx), y_val(idx_dx), 20, colors(dx_val,:), 'filled');
end
%xlabel('Notch3');
%ylabel('TMT score');
%legend({'CN','MCI','Dementia'});
%regression line
coeffs = polyfit(x_notch3, y_val, 1);  % linear fit
x_fit = linspace(min(x_notch3), max(x_notch3), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_notch3, y_val)

%% Scatter plot for a notch3 and SDS
i = find(strcmp(protein_ids, 'X5108_72'));
x_notch3 = x(:,i);
y_val = Y(:,2);
figure; hold on;
for dx_val = 1:3
    idx_dx = subject_dx == dx_val;
    scatter(x_notch3(idx_dx), y_val(idx_dx), 20, colors(dx_val,:), 'filled');
end
%xlabel('Notch3');
%ylabel('SDS score');
%legend({'CN','MCI','Dementia'});
%regression line
coeffs = polyfit(x_notch3, y_val, 1);  % linear fit
x_fit = linspace(min(x_notch3), max(x_notch3), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_notch3, y_val)

%% Scatter plot for nptxr and TMT
i = find(strcmp(protein_ids, 'X8997_4'));
x_nptxr = x(:,i);
y_val = Y(:,1);
figure; hold on;
for dx_val = 1:2
    idx_dx = subject_dx == dx_val;
    scatter(x_nptxr(idx_dx), y_val(idx_dx), 10, colors(dx_val,:), 'filled');
end
%xlabel('Notch3');
%ylabel('SDS score');
%legend({'CN','MCI','Dementia'});
%regression line
coeffs = polyfit(x_nptxr, y_val, 1);  % linear fit
x_fit = linspace(min(x_nptxr), max(x_nptxr), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_nptxr, y_val)

%% Scatter plot for nptxr and RAVLT
i = find(strcmp(protein_ids, 'X8997_4'));
x_nptxr = x(:,i);
y_val = Y(:,2);
figure; hold on;
for dx_val = 1:2
    idx_dx = subject_dx == dx_val;
    scatter(x_nptxr(idx_dx), y_val(idx_dx), 10, colors(dx_val,:), 'filled');
end
%xlabel('Notch3');
%ylabel('SDS score');
%legend({'CN','MCI','Dementia'});
%regression line
coeffs = polyfit(x_nptxr, y_val, 1);  % linear fit
x_fit = linspace(min(x_nptxr), max(x_nptxr), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_nptxr, y_val)

%% Scatter plot for ywhag and TMT
i = find(strcmp(protein_ids, 'X4179_57'));
x_ywhag = x(:,i);
y_val = Y(:,1);
figure; hold on;
for dx_val = 1:2
    idx_dx = subject_dx == dx_val;
    scatter(x_ywhag(idx_dx), y_val(idx_dx), 10, colors(dx_val,:), 'filled');
end
%xlabel('Notch3');
%ylabel('SDS score');
%legend({'CN','MCI','Dementia'});
%regression line
coeffs = polyfit(x_ywhag, y_val, 1);  % linear fit
x_fit = linspace(min(x_ywhag), max(x_ywhag), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_ywhag, y_val)

%% Scatter plot for ywhag and RAVLT
i = find(strcmp(protein_ids, 'X4179_57'));
x_ywhag = x(:,i);
y_val = Y(:,2);
figure; hold on;
for dx_val = 1:2
    idx_dx = subject_dx == dx_val;
    scatter(x_ywhag(idx_dx), y_val(idx_dx), 10, colors(dx_val,:), 'filled');
end
%xlabel('Notch3');
%ylabel('SDS score');
%legend({'CN','MCI','Dementia'});
%regression line
coeffs = polyfit(x_ywhag, y_val, 1);  % linear fit
x_fit = linspace(min(x_ywhag), max(x_ywhag), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_ywhag, y_val)

%% SCATTER PLOT FOR nptxr AND TMT
i = find(strcmp(protein_ids, 'X8997_4'))
x_nptxr = x(:,i);
y = Y(:,1);
figure;
dscatter(x_nptxr, y);
lsline;
corr(x_nptxr, y)

%% SCATTER PLOT FOR notch3 AND TMT
i = find(strcmp(protein_ids, 'X5108_72'))
x_notch3 = x(:,i);
y = Y(:,1);
figure;
dscatter(x_notch3, y);
lsline;
corr(x_notch3, y)

%% SCATTER PLOT FOR notch3 AND SDS
i = find(strcmp(protein_ids, 'X5108_72'))
x_notch3 = x(:,i);
y = Y(:,2);
figure;
dscatter(x_notch3, y);
lsline;
corr(x_notch3, y)

%% SCATTER PLOT FOR mmp12 AND TMT
i = find(strcmp(protein_ids, 'X4496_60'))
x_mmp12 = x(:,i)
y = Y(:,1)
figure;
dscatter(x_mmp12, y);
lsline;
corr(x_mmp12, y)


%% SCATTER PLOT FOR mmp12 AND SDS
i = find(strcmp(protein_ids, 'X4496_60'))
x_mmp12 = x(:,i)
y = Y(:,2)
figure;
dscatter(x_mmp12, y);
lsline;
corr(x_mmp12, y)

%% SCATTER PLOT FOR IL6 AND YS
%x = merged_data{:, 'il6'};
x = merged_data.il6(naninx == 0);
y = YS(:,3);

valid_idx = ~isnan(y);
x_clean = x(valid_idx);
y_clean = y(valid_idx);

% Create scatter plot
figure;
scatter(x_clean, y_clean, 5, 'k', 'filled');

% Optional: add linear regression line
mdl = fitlm(x_clean, y_clean);
[ypred,yci] = predict(mdl,x_clean, 'Alpha',0.001); hold on
plot(x_clean, ypred, 'k', 'LineWidth', 2);
plot(x_clean, yci(:,1), 'b', 'LineWidth', 0.5);
plot(x_clean, yci(:,2), 'b', 'LineWidth', 0.5);


corr(x_clean, y_clean)
corr(XS,Y)

%% SCATTER PLOT FOR Y SCORES VS RAW Y

yraw = merged_data.p23324_i2(naninx == 0);
y = YS(:,3);
% remove NaNs consistently
valid_idx = ~isnan(y) & ~isnan(yraw);
yraw = yraw(valid_idx);
y = y(valid_idx);

% scatter plot
figure;
scatter(yraw, y, 5, 'k', 'filled');
xlabel('Raw Y Data');
ylabel('Y-scores (Component 3)');
title('Raw Y vs. Y-scores');

corr(yraw, y)

%% Scatter plot for a cdcp1 and tmt
i = find(strcmp(protein_ids, 'X16818_200'));
x_cdcp1 = x(:,i);
y_val = Y(:,1);
figure; hold on;
for dx_val = 1:2
    idx_dx = subject_dx == dx_val;  
    scatter(x_cdcp1(idx_dx), y_val(idx_dx), 20, colors(dx_val,:), 'filled');
end
%xlabel('Notch3');
%ylabel('TMT score');
%legend({'CN','MCI','Dementia'});
%regression line
coeffs = polyfit(x_cdcp1, y_val, 1);  % linear fit
x_fit = linspace(min(x_cdcp1), max(x_cdcp1), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_cdcp1, y_val)
exportgraphics(gcf,'cdcp1.png','Resolution', 300);
% r = 0.0906

%% Scatter plot for a spink6 and tmt
i = find(strcmp(protein_ids, 'X5731_1'));
x_spink6 = x(:,i);
y_val = Y(:,1);
figure; hold on;
for dx_val = 1:2
    idx_dx = subject_dx == dx_val;  
    scatter(x_spink6(idx_dx), y_val(idx_dx), 20, colors(dx_val,:), 'filled');
end
%xlabel('Notch3');
%ylabel('TMT score');
%legend({'CN','MCI','Dementia'});
%regression line
coeffs = polyfit(x_cdcp1, y_val, 1);  % linear fit
x_fit = linspace(min(x_spink6), max(x_spink6), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_spink6, y_val)
exportgraphics(gcf,'spink6.png','Resolution', 300);
% r = 0.0984

%% Scatter plot for a crkl and tmt
i = find(strcmp(protein_ids, 'X16885_49'));
x_MAPRE3 = x(:,i);
y_val = Y(:,1);
figure; hold on;
for dx_val = 1:2
    idx_dx = subject_dx == dx_val;  
    scatter(x_MAPRE3(idx_dx), y_val(idx_dx), 20, colors(dx_val,:), 'filled');
end
coeffs = polyfit(x_MAPRE3, y_val, 1);  % linear fit
x_fit = linspace(min(x_MAPRE3), max(x_MAPRE3), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_MAPRE3, y_val)
exportgraphics(gcf,'MAPRE3.png','Resolution', 300);
% r = -0.1763








