clear
% Importing files
cd /Users/juliannehu/Documents/ADNI/data
mri=readtable('ADNI_Freesurfer_UCSFFSX7_26Feb2026.csv');
vitals = readtable('All_Subjects_VITALS_03Jul2026.csv');
disease = readtable('MEDHIST_11Aug2026.csv');
demo=readtable('ADNI_Participant_Demographics.csv');
cd /Users/juliannehu/Documents/ADNI/data/cleaned
proteomics=readtable('mean_imputed_somascan_proteomics');
cd /Users/juliannehu/Documents/ADNI/data
diagnosis=readtable('DXSUM.csv');

%Extracting diagnosis with dates
diagnosis = diagnosis(:, {'RID', 'DIAGNOSIS', 'EXAMDATE'});

% Joining files
%pro_demo = innerjoin(proteomics, demo, 'Keys', {'RID', 'VISCODE2'});
covariate_data = innerjoin(vitals, disease, 'Keys', 'RID');
covariate_data = innerjoin(covariate_data, demo, 'Keys', 'RID');
pro_demo = innerjoin(proteomics, covariate_data,...
                     'Keys', {'RID'});
%merged_data = innerjoin(pro_demo, cognitive, 'Keys', {'RID', 'VISCODE'});
merged_data = innerjoin(pro_demo, mri, ...
                         'Keys', {'RID'});

%Calculating Age
% Convert PTDOB (MM/yyyy) to datetime
merged_data.PTDOB = datetime(merged_data.PTDOB, 'InputFormat', 'MM/yyyy');
% Convert proteomics examdate to datetime
merged_data.EXAMDATE_pro_demo = datetime(merged_data.EXAMDATE_pro_demo, 'InputFormat', 'yyyy-MM-dd');
% Calculate age at proteomics visit
merged_data.AGE_proteomics = years(merged_data.EXAMDATE_pro_demo - merged_data.PTDOB);
%Convert mri examdate to datetime
merged_data.EXAMDATE_mri = datetime(merged_data.EXAMDATE_mri, 'InputFormat', 'yyyy-MM-dd');
% Calculate age at mri visit
merged_data.AGE_mri = years(merged_data.EXAMDATE_mri - merged_data.PTDOB);
%demo.AgeAtVisit = floor(demo.AgeAtVisit);

clear mri proteomics demo;

% Keep only rows with 3T scans
mri_3t_idx = strcmp(merged_data.FIELD_STRENGTH, '3T');
merged_data = merged_data(mri_3t_idx,:);

% Exclude rows that failed overall QC
mri_qc_idx = ~strcmp(merged_data.OVERALLQC, 'Fail');
merged_data = merged_data(mri_qc_idx,:);

% Exclude 2 MRI variables that have almost all NaN values
merged_data = removevars(merged_data, {'ST8SV','ST68SV'});


% Specifying date format
%merged_data.EXAMDATE = datetime(merged_data.EXAMDATE, 'InputFormat', 'MM-dd-yyyy');
%merged_data.VISDATE_cognitive = datetime(merged_data.VISDATE_cognitive, 'InputFormat', 'yyyy-MM-dd');
% Calculate change in time
merged_data.time_delta_pro_mri=(merged_data.EXAMDATE_mri-merged_data.EXAMDATE_pro_demo); 

%Only keeping correct timepoint rows
uniqueIDs = unique(merged_data.RID);
rows_to_keep = false(height(merged_data),1);
for i = 1:length(uniqueIDs)
    rid = uniqueIDs(i);
    % Find rows belonging to single participant
    idx = merged_data.RID == rid;
    deltas = merged_data.time_delta_pro_mri(idx);
    % Keep only NON-negative deltas
    valid = deltas >= 0;
    if any(valid)
        % Find the minimum delta time in those rows to keep
        [~, minIdx] = min(deltas(valid));
        rid_rows = find(idx);
        rid_valid_rows = rid_rows(valid);
        rows_to_keep(rid_valid_rows(minIdx)) = true;
    else
        %do nothing
    end
end
merged_data = merged_data(rows_to_keep, :);


% Convert diagnosis EXAMDATE to datetime
diagnosis.EXAMDATE = datetime(diagnosis.EXAMDATE, 'InputFormat', 'yyyy-MM-dd');
% Sort by RID and examdate
diagnosis = sortrows(diagnosis, {'RID', 'EXAMDATE'});

%add in diagnosis column to merge file based on dates that match closest together 
diagnosis_pro = NaN(height(merged_data),1);
diagnosis_mri = NaN(height(merged_data),1);

for i = 1:height(merged_data)
    rid = merged_data.RID(i);
    
    % Search all diagnosis for this RID
    rows = diagnosis.RID == rid;

    if any(rows)
        diagnosis_dates  = diagnosis.EXAMDATE(rows);
        diagnosis_labels = diagnosis.DIAGNOSIS(rows);

        % Diagnosis at proteomics timepoint
        pro_date = merged_data.EXAMDATE_pro_demo(i);
        valid_idx = diagnosis_dates >= pro_date;
        if any(valid_idx)
            [~, idx_x] = min(diagnosis_dates(valid_idx) - pro_date);
            temp_labels = diagnosis_labels(valid_idx);
            diagnosis_pro(i) = temp_labels(idx_x);
        else
            diagnosis_pro(i) = NaN;  % no diagnosis on/after visit
        end

        % Diagnosis at mri timepoint
        mri_date = merged_data.EXAMDATE_mri(i);
        valid_idx = diagnosis_dates >= mri_date;  % only dates >= visit
        if any(valid_idx)
            [~, idx_y] = min(diagnosis_dates(valid_idx) - mri_date);
            temp_labels = diagnosis_labels(valid_idx);
            diagnosis_mri(i) = temp_labels(idx_y);
        else
            diagnosis_mri(i) = NaN;
        end
    end
end

% Add to merged_data
merged_data.DIAGNOSIS_PRO = diagnosis_pro;
merged_data.DIAGNOSIS_MRI = diagnosis_mri;

% Exclude participants with dementia (DIAGNOSIS = 3)
mri_dx_idx = merged_data.DIAGNOSIS_PRO ~= 3 & merged_data.DIAGNOSIS_MRI ~= 3;
merged_data = merged_data(mri_dx_idx,:);


% Variables
age_pro=merged_data.AGE_proteomics; 
age_mri=merged_data.AGE_mri; 
sex=merged_data.PTGENDER; 
site=merged_data.SITEID;
icv=merged_data.ST10CV;
%motion=merged_data.
dx_pro=merged_data.DIAGNOSIS_PRO;
dx_mri=merged_data.DIAGNOSIS_MRI;
delta_time=merged_data.time_delta_pro_mri;
%convert delta time to numeric
delta_time_num = days(delta_time); % numeric, in days
education=merged_data.PTEDUCAT;
cvd=merged_data.MH4CARD;
%calculate bmi
%weight calucation
merged_data.weight = NaN(height(merged_data), 1);
%convert lbs to kg
lb_idx = merged_data.VSWTUNIT == 1;
merged_data.weight(lb_idx) = merged_data.VSWEIGHT(lb_idx) * 0.45359237;
%keep kilograms as is
kg_idx = merged_data.VSWTUNIT == 2;
merged_data.weight(kg_idx) = merged_data.VSWEIGHT(kg_idx);
%height calculation
merged_data.height = NaN(height(merged_data), 1);
%convert in to m
in_idx = merged_data.VSHTUNIT == 1;
merged_data.height(in_idx) = merged_data.VSHEIGHT(in_idx) * 0.0254;
% convert cm to m
cm_idx = merged_data.VSHTUNIT == 2;
merged_data.height(cm_idx) = merged_data.VSHEIGHT(cm_idx) / 100;
%bmi calculation
merged_data.BMI = merged_data.weight ./ (merged_data.height .^ 2);
%define bmi variable
bmi=merged_data.BMI;

%Replace missing data for sex and site
sex(sex==-4) = NaN;
site(site==-4) = NaN;
% X = proteomics
X = merged_data{:, 8:7008}; 
protein_ids=merged_data.Properties.VariableNames(8:7008);

% Y: cortical thickness (ends with TA) + subcortical volume (ends with SV)
% get all variable names
var_names = merged_data.Properties.VariableNames;
% Find columns ending with TA and specific subcortical volumes
idx = endsWith(var_names, 'TA') | strcmp(var_names,'ST112SV') | strcmp(var_names,'ST53SV') | strcmp(var_names,'ST11SV') | strcmp(var_names,'ST70SV') | strcmp(var_names,'ST120SV')| strcmp(var_names,'ST61SV') | strcmp(var_names,'ST12SV')| strcmp(var_names,'ST71SV')|...
    strcmp(var_names,'ST71SV') | strcmp(var_names,'ST16SV') | strcmp(var_names,'ST75SV') | strcmp(var_names,'ST29SV') | strcmp(var_names,'ST88SV')|...
    strcmp(var_names,'ST127SV') |strcmp(var_names,'ST30SV')|strcmp(var_names,'ST37SV')|strcmp(var_names,'ST89SV')|strcmp(var_names,'ST8SV')|strcmp(var_names,'ST96SV')|strcmp(var_names,'ST9SV');

% Extract Y (mri data)
Y = merged_data{:, idx};
% Save the names of the Y variables
mri_var = var_names(idx);



%%% Removing NaNs
%naninx=sum(isnan(Y)')'; 
%Y=Y(naninx==0,:); 
%X=X(naninx==0,:); 
%age=age(naninx==0); 
%sex=sex(naninx==0);
%site=site(naninx==0); 
rows_keep = ~any(ismissing([X Y age_pro age_mri dx_pro dx_mri delta_time_num sex site education cvd bmi]), 2);
X=X(rows_keep, :);
Y=Y(rows_keep, :);
age_pro=age_pro(rows_keep);
age_mri=age_mri(rows_keep);
dx_pro=dx_pro(rows_keep);
dx_mri=dx_mri(rows_keep);
delta_time_num=delta_time_num(rows_keep);
sex=sex(rows_keep);
site=site(rows_keep);
icv=icv(rows_keep);
education=education(rows_keep);
cvd=cvd(rows_keep);
bmi=bmi(rows_keep);
%motion=motion(rows_keep);

%Convert to categorical
sex=categorical(sex);
site=categorical(site);
dx_pro=categorical(dx_pro);
dx_mri=categorical(dx_mri);
cvd=categorical(cvd);

% Standardizing data
Y=zscore(Y); 
X=zscore(X); 


% Set # of PLS components
ncomp=1;

%defining other covariates
age_pro_squared = age_pro.^2;
age_mri_squared = age_mri.^2;

%make sex num binary 0/1
sex_num = double(string(sex)) - 1;
age_pro_by_sex = age_pro .* sex_num;
age_mri_by_sex = age_mri .* sex_num;

% Regress out age, sex, and site (covariates)
%mdl = fitlm(table(age_mri, dx_cog, sex, site,Y(:,1)), 'CategoricalVars',[2 3]); Y(:,1)=mdl.Residuals.Raw;
%mdl = fitlm(table(age_mri, dx_cog, sex, site,Y(:,2)), 'CategoricalVars',[2 3]); Y(:,2)=mdl.Residuals.Raw;
%mdl = fitlm(table(age, sex, site,Y(:,3)) ); Y(:,3)=mdl.Residuals.Raw;
%mdl = fitlm(table(i2age, sex, site,Y(:,4)) ); Y(:,4)=mdl.Residuals.Raw;
clear x; 
for i=1:length(X(1,:)) 
    mdl = fitlm(table(delta_time_num, dx_pro, sex, site, bmi, education, cvd, X(:,i))); 
    x(:,i)=mdl.Residuals.Raw;
end

clear y; 
for i=1:length(Y(1,:)); 
    mdl = fitlm(table(delta_time_num, dx_mri, sex, site, icv, bmi, education, cvd, Y(:,i))); 
    y(:,i)=mdl.Residuals.Raw;
end

% Run PLS
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);
PCTVAR

% Plotting each of the components' correlation with cognitive tests
figure;
imagesc(corr(XS(:,1:2),Y)); 
colormap bone; 
colorbar; 
corr(XS, Y) %1.TMT  2.GF  3.PAL  4.DSST
combs=allcomb([1:3], [1:2]); 
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

% Calculate correlations
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


% figure to show the permutation distribution vs actual distribution
figure; 
histogram(Rsq, 'FaceColor', [0.5 0.5 0.5]); 
xlim([0 0.09]); 
xline(0.0501, 'r', 'LineWidth', 2); 
% Only show x and y axes
ax = gca;
ax.Box = 'off';          % removes top and right lines
ax.TickDir = 'out';      % ticks point outward
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
sum(PCTVAR')

%% bootstrapping to get the func connectivity weights for PLS1, 2 and 3
dim=1
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,dim);PCTVAR
PLS1w=stats.W(:,1);
%PLS2w=stats.W(:,2);
%PLS3w=stats.W(:,3);

bootnum=5000;
PLS1weights=[];
%PLS2weights=[];
%PLS3weights=[];

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
    
   
   % newW=stats.W(:,2);%extract PLS2 weights
    %if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
     %   newW=-1*newW;
    %end
    %PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run    
    
    %newW=stats.W(:,3);%extract PLS2 weights
    %if corr(PLS3w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
     %   newW=-1*newW;
    %end
    %PLS3weights=[PLS3weights,newW]; %store (ordered) weights from this bootstrap run    
end

PLS1sw=std(PLS1weights');
%PLS2sw=std(PLS2weights');
%PLS3sw=std(PLS3weights');

plsweights1=PLS1w./PLS1sw';
%plsweights2=PLS2w./PLS2sw'; 
%plsweights3=PLS3w./PLS3sw';

%flip sign of PLS 1, 2 and 3 weights to match upregulated vs downregulated
%plsweights1 = plsweights1*-1;
%stats.W(:,1) = -stats.W(:,1);
%plsweights2 = plsweights2*-1;
%stats.W(:,2) = -stats.W(:,2);
%plsweights3 = plsweights3*-1;
%stats.W(:,3) = -stats.W(:,3);
%flipping other PLS stats for consistency
%XS(:,1) = -XS(:,1);
%XS(:,2) = -XS(:,2);
%XS(:,3) = -XS(:,3);
%XL(:,1) = -XL(:,1);
%XL(:,2) = -XL(:,2);
%XL(:,3) = -XL(:,3);
%YS(:,1) = -YS(:,1);
%YS(:,2) = -YS(:,2);
%YS(:,3) = -YS(:,3);
%YL(:,1) = -YL(:,1);
%YL(:,2) = -YL(:,2); 
%YL(:,3) = -YL(:,3); 

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


cd /Users/juliannehu/Documents/ADNI/data/dictionaries
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
cd /Users/juliannehu/Documents/ADNI/results/pro-mri/
writetable(PLS_table, 'ADNI_PLS_protein_weights_with_age_ventricles_no_dementia.csv');






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



%VOLCANO PLOTS FOR PLS 1 WITH TRANSPARENCY FOR NON SIG
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
%xlabel('PLS Weights (I0-I2 Component 3 Age Regressed)');
ylabel('-log_{10}(p)');

% Highlight protein X3520.58 in yellow
target = "X3520_58";
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
sig_threshold = 2*(1 - normcdf(3));

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
    sig_idx = idx & (p_pls2 <= sig_threshold);
    scatter(stats.W(sig_idx,2), -log10(p_pls2(sig_idx)), 10, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);

    % Points below significance
    nonsig_idx = idx & (p_pls2 > sig_threshold);
    scatter(stats.W(nonsig_idx,2), -log10(p_pls2(nonsig_idx)), 10, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
end

yline(-log10(2*(1 - normcdf(3))), 'k--', 'LineWidth', 1.2); 
%xlabel('PLS Weights (I0-I2 Component 3 Age Regressed)');
ylabel('-log_{10}(p)');


%VOLCANO PLOTS FOR PLS 3 WITH TRANSPARENCY FOR NON SIG
p_pls3 = (1 - normcdf(abs(plsweights3)));
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
    sig_idx = idx & (p_pls3 <= sig_threshold);
    scatter(stats.W(sig_idx,3), -log10(p_pls3(sig_idx)), 10, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);

    % Points below significance
    nonsig_idx = idx & (p_pls3 > sig_threshold);
    scatter(stats.W(nonsig_idx,3), -log10(p_pls3(nonsig_idx)), 10, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
end

yline(-log10(1 - normcdf(3)), 'k--', 'LineWidth', 1.2); 
%xlabel('PLS Weights (I0-I2 Component 3 Age Regressed)');
ylabel('-log_{10}(p)');





%%%% VISUALIZATIONS ----- NO COLOURS %%%%
%% Volcano plot for PLS 1
p_pls1 = 2 * (1 - normcdf(abs(plsweights1)));
% Define significance threshold
sig_threshold = 2*(1 - normcdf(3));
% Set colors
color_sig = [0 0 0];  
color_nonsig = [0 0 0];
%Significance
sig_idx = p_pls1 <= sig_threshold;
nonsig_idx = p_pls1 > sig_threshold;
figure; hold on;
%Plot significant metabolites
scatter(stats.W(sig_idx,1),-log10(p_pls1(sig_idx)), 10, ...
    'MarkerFaceColor', color_sig, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
%Plot nonsignificant metabolites
scatter(stats.W(nonsig_idx,1),-log10(p_pls1(nonsig_idx)), 10, ...
    'MarkerFaceColor', color_nonsig, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
yline(-log10(sig_threshold), 'k--', 'LineWidth', 1.2); 
%xlabel('PLS Weights (I0-I2 PLS 1 Age Regressed)');
ylabel('-log_{10}(p)');


%% Volcano plot for PLS 2
p_pls2 = 2 * (1 - normcdf(abs(plsweights2)));
% Define significance threshold
sig_threshold = 2*(1 - normcdf(3));
% Set colors
color_sig = [0.2 0.6 0.8];  
color_nonsig = [0.2 0.6 0.8];
%Significance
sig_idx = p_pls2 <= sig_threshold;
nonsig_idx = p_pls2 > sig_threshold;
figure; hold on;
%Plot significant metabolites
scatter(stats.W(sig_idx,2),-log10(p_pls2(sig_idx)), 10, ...
    'MarkerFaceColor', color_sig, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
%Plot nonsignificant metabolites
scatter(stats.W(nonsig_idx,2),-log10(p_pls2(nonsig_idx)), 10, ...
    'MarkerFaceColor', color_nonsig, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
yline(-log10(sig_threshold), 'k--', 'LineWidth', 1.2); 
%xlabel('PLS Weights (I0-I2 PLS 1 Age Regressed)');
ylabel('-log_{10}(p)');


%% Volcano plot for PLS 3
p_pls3 = 2 * (1 - normcdf(abs(plsweights3)));
% Define significance threshold
sig_threshold = 2*(1 - normcdf(3));
% Set colors
color_sig = [0.2 0.6 0.8];  
color_nonsig = [0.2 0.6 0.8];
%Significance
sig_idx = p_pls3 <= sig_threshold;
nonsig_idx = p_pls3 > sig_threshold;
figure; hold on;
%Plot significant proteins
scatter(stats.W(sig_idx,3),-log10(p_pls3(sig_idx)), 10, ...
    'MarkerFaceColor', color_sig, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
%Plot nonsignificant proteins
scatter(stats.W(nonsig_idx,3),-log10(p_pls3(nonsig_idx)), 10, ...
    'MarkerFaceColor', color_nonsig, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
yline(-log10(sig_threshold), 'k--', 'LineWidth', 1.2); 
%xlabel('PLS Weights (I0-I2 PLS 1 Age Regressed)');
ylabel('-log_{10}(p)');




%FIND KEY PROTEINS FOR POST PROCESSING ON POWERPOINT
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
%% sort by |PLS weight| (strongest effect at top)
up_tbl   = sortrows(up_tbl, 'PLS_Weight', 'descend');  
down_tbl = sortrows(down_tbl, 'PLS_Weight');           

%% Display results
disp('Upregulated Proteins:');
disp(up_tbl);

disp('Downregulated Proteins:');
disp(down_tbl);

writetable(up_tbl, 'adni_csf_proteomics_cognitive_pls1_significant_upregulated_proteins_no_age.csv');
writetable(down_tbl, 'adni_csf_proteomics_cognitive_pls1_significant_downregulated_proteins_no_age.csv');



%PLS 2
sig_threshold = 2*(1 - normcdf(3));
is_sig = p_pls2 <= sig_threshold;

up_mask = stats.W(:,2) > 0 & is_sig;
down_mask = stats.W(:,2) < 0 & is_sig;

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
               stats.W(up_mask,2), ...
               p_pls2(up_mask), ...
               'VariableNames', {'Protein', 'Gene', 'Category', 'PLS_Weight','pValue'});

down_tbl = table(protein_names(down_mask), ...
                 gene_names_str(down_mask), ...
                 protein_category(down_mask), ...
                 stats.W(down_mask,2), ...
                 p_pls2(down_mask), ...
                 'VariableNames', {'Protein', 'Gene', 'Category', 'PLS_Weight','pValue'});

%% sort by |PLS weight| (strongest effect at top)
up_tbl   = sortrows(up_tbl, 'PLS_Weight', 'descend');  
down_tbl = sortrows(down_tbl, 'PLS_Weight');           

%% Display results
disp('Upregulated Proteins:');
disp(up_tbl);

disp('Downregulated Proteins:');
disp(down_tbl);

writetable(up_tbl, 'adni_csf_proteomics_cognitive_pls2_significant_upregulated_proteins_no_age.csv');
writetable(down_tbl, 'adni_csf_proteomics_cognitive_pls2_significant_downregulated_proteins_no_age.csv');



%PLS 3
sig_threshold = 2*(1 - normcdf(3));
is_sig = p_pls3 <= sig_threshold;

up_mask = stats.W(:,3) > 0 & is_sig;
down_mask = stats.W(:,3) < 0 & is_sig;

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
               stats.W(up_mask,3), ...
               p_pls3(up_mask), ...
               'VariableNames', {'Protein', 'Gene', 'Category', 'PLS_Weight','pValue'});

down_tbl = table(protein_names(down_mask), ...
                 gene_names_str(down_mask), ...
                 protein_category(down_mask), ...
                 stats.W(down_mask,3), ...
                 p_pls3(down_mask), ...
                 'VariableNames', {'Protein', 'Gene', 'Category', 'PLS_Weight','pValue'});

%% sort by |PLS weight| (strongest effect at top)
up_tbl   = sortrows(up_tbl, 'PLS_Weight', 'descend');  
down_tbl = sortrows(down_tbl, 'PLS_Weight');           

%% Display results
disp('Upregulated Proteins:');
disp(up_tbl);

disp('Downregulated Proteins:');
disp(down_tbl);

writetable(up_tbl, 'adni_csf_proteomics_cognitive_pls3_significant_upregulated_proteins_no_age.csv');
writetable(down_tbl, 'adni_csf_proteomics_cognitive_pls3_significant_downregulated_proteins_no_age.csv');



%SCATTER PLOT FOR PLS 1 XS AND YS WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(age_mri-min(age_mri))./max(age_mri-min(age_mri));
colormap_custom(:,3)=1-age_mri./max(age_mri);
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
c.TickLabels = round(linspace(min(age_mri), max(age_mri), 5));



%SCATTER PLOT FOR PLS 2 XS AND YS WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(age_mri-min(age_mri))./max(age_mri-min(age_mri));
colormap_custom(:,3)=1-age_mri./max(age_mri);
figure;
ix1=2; ix2=2
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
c.TickLabels = round(linspace(min(age_mri), max(age_mri), 5));



%SCATTER PLOT FOR PLS 3 XS AND YS WITH COLOURSCHEME BY AGE
colormap_custom=zeros([length(XS(:,1)), 3]);
colormap_custom(:,1)=(age_mri-min(age_mri))./max(age_mri-min(age_mri));
colormap_custom(:,3)=1-age_mri./max(age_mri);
figure;
ix1=3; ix2=3;
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
c.TickLabels = round(linspace(min(age_mri), max(age_mri), 5));






%SCATTER PLOT FOR PLS 1 XS and YS 
figure;
ix1=1; ix2=1;
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
ix1=2; ix2=2;
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
ix1=3; ix2=3;
scatter(XS(:,ix1), YS(:,ix2), 5,'k', 'filled'); 
mdl= fitlm(XS(:, ix1), YS(:,ix2)); [ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
plot(XS(:, ix1), ypred, 'k', 'LineWidth', 2);
plot(XS(:, ix1), yci(:,1), 'b', 'LineWidth', 0.5);
plot(XS(:, ix1), yci(:,2), 'b', 'LineWidth', 0.5);
%xlabel(sprintf('X Scores (LV%d)', ix1), 'FontSize', 12);
%ylabel(sprintf('Y Scores (LV%d)', ix2), 'FontSize', 12);
corr(XS(:,ix1), YS(:,ix2))








%SCATTER PLOTS WITH COLOURS BASED ON DIAGNOSES
cd /Users/juliannehu/Documents/ADNI/data
diagnosis=readtable('DXSUM.csv');

all_rids = merged_data.RID;
rid_vector = all_rids(rows_keep);  % 240 rows
[isFound, idx] = ismember(rid_vector, diagnosis.RID);
subject_dx = diagnosis.DIAGNOSIS(idx);  % 1=CN, 2=MCI, 3=Dementia
[isFound, idx] = ismember(rid_vector, diagnosis.RID);  
%% Define colors
colors = [0 0 1;      % CN = blue
          0.9 0.6 0.1; % MCI = orange
          0.8 0.1 0.1]; % Dementia = red


%Scatter plot for nefl and YS1
i = find(strcmp(protein_ids, 'X10082_251'));
x_val = x(:,i);
y_val = YS(:,1);
figure; hold on;
for dx_val = 1:3
    idx_dx = subject_dx == dx_val;
    scatter(x_val(idx_dx), y_val(idx_dx), 20, colors(dx_val,:), 'filled');
end
%xlabel('Notch3');
%ylabel('SDS score');
%legend({'CN','MCI','Dementia'});
%regression line
coeffs = polyfit(x_val, y_val, 1);  % linear fit
x_fit = linspace(min(x_val), max(x_val), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_val, y_val)


%Scatter plot for nefl and left superior temporal
i = find(strcmp(protein_ids, 'X10082_251'));
x_val = x(:,i);
j = find(strcmp(mri_region, 'Thickness Average (aparc.stats) of LeftSuperiorTemporal'))
y_val = Y(:,j);
figure; hold on;
for dx_val = 1:3
    idx_dx = subject_dx == dx_val;
    scatter(x_val(idx_dx), y_val(idx_dx), 20, colors(dx_val,:), 'filled');
end
%xlabel('Notch3');
%ylabel('SDS score');
%legend({'CN','MCI','Dementia'});
%regression line
coeffs = polyfit(x_val, y_val, 1);  % linear fit
x_fit = linspace(min(x_val), max(x_val), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_val, y_val)


%Scatter plot for tgfb3 and YS1
i = find(strcmp(protein_ids, 'X3520_58'));
x_val = x(:,i);
y_val = YS(:,1);
figure; hold on;
for dx_val = 1:3
    idx_dx = subject_dx == dx_val;
    scatter(x_val(idx_dx), y_val(idx_dx), 20, colors(dx_val,:), 'filled');
end
%xlabel('Notch3');
%ylabel('SDS score');
%legend({'CN','MCI','Dementia'});
%regression line
coeffs = polyfit(x_val, y_val, 1);  % linear fit
x_fit = linspace(min(x_val), max(x_val), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_val, y_val)


%Scatter plot for tgfb3 and left superior temporal
i = find(strcmp(protein_ids, 'X3520_58'));
x_val = x(:,i);
j = find(strcmp(mri_region, 'Thickness Average (aparc.stats) of LeftSuperiorTemporal'))
y_val = Y(:,j);
figure; hold on;
for dx_val = 1:3
    idx_dx = subject_dx == dx_val;
    scatter(x_val(idx_dx), y_val(idx_dx), 20, colors(dx_val,:), 'filled');
end
%xlabel('Notch3');
%ylabel('SDS score');
%legend({'CN','MCI','Dementia'});
%regression line
coeffs = polyfit(x_val, y_val, 1);  % linear fit
x_fit = linspace(min(x_val), max(x_val), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'k-');  % regression line
corr_coef = corr(x_val, y_val)





%SCATTER PLOT FOR Y SCORES VS RAW Y

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





% PREPARING DATA FOR PLOTTING RESULTS WITH GGSEG
[r,p]=corr(XS, Y)
r = r.';
p = p.';
cd /Users/juliannehu/Documents/ADNI/data/dictionaries
data_dictionary=readtable('ADNI_Freesurfer_UCSFFSX7_Data_Dictionary');
% Preallocate a cell array to hold the descriptions
mri_region = cell(length(mri_var),1);
% Loop through each MRI variable
for i = 1:length(mri_var)
    idx = strcmp(data_dictionary.FLDNAME, mri_var{i});  % find matching FLDNAME
    if any(idx)
        mri_region{i} = data_dictionary.TEXT{idx};    % get the corresponding TEXT
    else
        mri_region{i} = '';  % or NaN if no match found
    end
end
 
%mri_var = mri_var(:); 
for i = 1:ncomp
    % Extract the correlations and p-values for component c
    correlation = r(:, i);
    pvalue = p(:, i);
    mri = mri_region;
    %mri = repmat(i, length(mri_variables), 1);  % k x 1
    % Build the table
    mri_results = table(mri, correlation, pvalue, ...
        'VariableNames', {'BrainRegion','correlation','pvalue'});
    % Write to CSV
    cd /Users/juliannehu/Documents/ADNI/results/pro-mri
    writetable(mri_results, sprintf('1mri_results_i0i2_with_age_ventricles_no_dementia_pls%d.csv', i));
end



%% MEDIATION MODEL %%
X_mediation = X; 
%X_mediation = XS(:1);
Y_mediation = YS(:,1);
M_mediation = age_mri;

%[paths, stats] = mediation(X_mediation, Y_mediation, M_mediation, 'plots', 'verbose', 'boot', 'bootreps', 5000);

ix=plsweights1>3.3; sum(ix)
ix=plsweights1<-3.3; sum(ix)
ix=abs(plsweights1)>3.3; sum(ix)
[paths, stats] = mediation(X_mediation(:,ix), Y_mediation, M_mediation, 'plots', 'verbose', 'boot', 'bootreps', 5000);