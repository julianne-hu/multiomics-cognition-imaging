clear

cd /Users/arshiya/Documents/CLSA/2601007_CAMH_PZhukovsky_BL/
data = readtable('2601007_CAMH_PZhukovsky_Baseline_CoPv7-2.csv');

data.BLD_IL6_COM(data.BLD_IL6_COM < -5) = NaN;
data.BLD_TNF_COM(data.BLD_TNF_COM < -5) = NaN;
data.COG_REYII_SCORE_COM(data.COG_REYII_SCORE_COM < -5) = NaN;
data.COG_MAT_SCORE_COM(data.COG_MAT_SCORE_COM < -5) = NaN;
data.BLD_HDL_COM(data.BLD_HDL_COM < -5) = NaN;
data.BLD_LDL_COM(data.BLD_LDL_COM < -5) = NaN;
data.BLD_nonHDL_COM(data.BLD_nonHDL_COM < -5) = NaN;
data.BLD_CHOL_COM(data.BLD_CHOL_COM < -5) = NaN;
tnf_data.BLD_TNF_COM(tnf_data.BLD_TNF_COM < -5) = NaN;
tnf_data.COG_REYII_SCORE_COM(tnf_data.COG_REYII_SCORE_COM < -5) = NaN;


%% glm models

%IL6 w/ REY
mdl = fitlm(data, 'COG_REYII_SCORE_COM ~ BLD_IL6_COM + AGE_NMBR_COM + SEX_ASK_COM + GEOSTRATA_COM')

%IL6 w/ MAT
mdl = fitlm(data, 'COG_MAT_SCORE_COM ~ BLD_IL6_COM + AGE_NMBR_COM + SEX_ASK_COM + GEOSTRATA_COM')

%TNF w/ REY
mdl = fitlm(tnf_data, 'COG_REYII_SCORE_COM ~ BLD_TNF_COM + AGE_NMBR_COM + SEX_ASK_COM + GEOSTRATA_COM')

%TNF w/ MAT
mdl = fitlm(tnf_data, 'COG_MAT_SCORE_COM ~ BLD_TNF_COM + AGE_NMBR_COM + SEX_ASK_COM + GEOSTRATA_COM')

%HDL w/ REY
mdl = fitlm(data, 'COG_REYII_SCORE_COM ~ BLD_HDL_COM + AGE_NMBR_COM + SEX_ASK_COM + GEOSTRATA_COM')

%HDL w/ MAT
mdl = fitlm(data, 'COG_MAT_SCORE_COM ~ BLD_HDL_COM + AGE_NMBR_COM + SEX_ASK_COM + GEOSTRATA_COM')

%LDL w/ REY
mdl = fitlm(data, 'COG_REYII_SCORE_COM ~ BLD_LDL_COM + AGE_NMBR_COM + SEX_ASK_COM + GEOSTRATA_COM')

%LDL w/ MAT
mdl = fitlm(data, 'COG_MAT_SCORE_COM ~ BLD_LDL_COM + AGE_NMBR_COM + SEX_ASK_COM + GEOSTRATA_COM')

%cholesterol w/ REY
mdl = fitlm(data, 'COG_REYII_SCORE_COM ~ BLD_CHOL_COM + AGE_NMBR_COM + SEX_ASK_COM + GEOSTRATA_COM')

%cholesterol w/ MAT
mdl = fitlm(data, 'COG_MAT_SCORE_COM ~ BLD_CHOL_COM + AGE_NMBR_COM + SEX_ASK_COM + GEOSTRATA_COM')

%nonHDL w/ REY
mdl = fitlm(data, 'COG_REYII_SCORE_COM ~ BLD_nonHDL_COM + AGE_NMBR_COM + SEX_ASK_COM + GEOSTRATA_COM')

%nonHDL w/ MAT
mdl = fitlm(data, 'COG_MAT_SCORE_COM ~ BLD_nonHDL_COM + AGE_NMBR_COM + SEX_ASK_COM + GEOSTRATA_COM')

%% LDL vs. REYII
missingCodes = [-8888, -9999, -1111];

x = data.BLD_LDL_COM;
y = data.COG_REYII_SCORE_COM;
age = data.AGE_NMBR_COM;
sex_num = double(categorical(data.SEX_ASK_COM));
geo_cat = categorical(data.GEOSTRATA_COM);
geo_dummies_full = dummyvar(geo_cat);
geo_dummies = geo_dummies_full(:, 2:end);

x(ismember(x, missingCodes)) = NaN;
y(ismember(y, missingCodes)) = NaN;

z = [age, sex_num, geo_dummies];

valid = isfinite(x) & isfinite(y) & all(isfinite(z), 2);
x = x(valid);
y = y(valid);
z = z(valid, :);
age = age(valid);          

colormap_custom = zeros([length(x), 3]);
colormap_custom(:,1) = (age - min(age)) ./ max(age - min(age));
colormap_custom(:,3) = 1 - age ./ max(age);

figure;
scatter(x, y, 5, colormap_custom, 'filled');
lsline;

[R, P] = corr(x, y)
[rho, p_partial] = partialcorr(x, y, z)

%% HDL vs REYII
missingCodes = [-8888, -9999, -1111];

x = data.BLD_HDL_COM;
y = data.COG_REYII_SCORE_COM;
age = data.AGE_NMBR_COM;
sex_num = double(categorical(data.SEX_ASK_COM));
geo_cat = categorical(data.GEOSTRATA_COM);
geo_dummies_full = dummyvar(geo_cat);
geo_dummies = geo_dummies_full(:, 2:end);

x(ismember(x, missingCodes)) = NaN;
y(ismember(y, missingCodes)) = NaN;

z = [age, sex_num, geo_dummies];

valid = isfinite(x) & isfinite(y) & all(isfinite(z), 2);
x = x(valid);
y = y(valid);
z = z(valid, :);
age = age(valid);          

colormap_custom = zeros([length(x), 3]);
colormap_custom(:,1) = (age - min(age)) ./ max(age - min(age));
colormap_custom(:,3) = 1 - age ./ max(age);

figure;
scatter(x, y, 5, colormap_custom, 'filled');
lsline;


[R, P] = corr(x, y)
[rho, p_partial] = partialcorr(x, y, z)

%% nonHDL vs REYII
missingCodes = [-8888, -9999, -1111];

x = data.BLD_nonHDL_COM;
y = data.COG_REYII_SCORE_COM;
age = data.AGE_NMBR_COM;
sex_num = double(categorical(data.SEX_ASK_COM));
geo_cat = categorical(data.GEOSTRATA_COM);
geo_dummies_full = dummyvar(geo_cat);
geo_dummies = geo_dummies_full(:, 2:end);

x(ismember(x, missingCodes)) = NaN;
y(ismember(y, missingCodes)) = NaN;

z = [age, sex_num, geo_dummies];

valid = isfinite(x) & isfinite(y) & all(isfinite(z), 2);
x = x(valid);
y = y(valid);
z = z(valid, :);
age = age(valid);          % <-- add this line, the missing fix

colormap_custom = zeros([length(x), 3]);
colormap_custom(:,1) = (age - min(age)) ./ max(age - min(age));
colormap_custom(:,3) = 1 - age ./ max(age);

figure;
scatter(x, y, 5, colormap_custom, 'filled');
lsline;


[R, P] = corr(x, y)
[rho, p_partial] = partialcorr(x, y, z)

%% cholesterol vs REYII
missingCodes = [-8888, -9999, -1111];

x = data.BLD_CHOL_COM;
y = data.COG_REYII_SCORE_COM;
age = data.AGE_NMBR_COM;
sex_num = double(categorical(data.SEX_ASK_COM));
geo_cat = categorical(data.GEOSTRATA_COM);
geo_dummies_full = dummyvar(geo_cat);
geo_dummies = geo_dummies_full(:, 2:end);

x(ismember(x, missingCodes)) = NaN;
y(ismember(y, missingCodes)) = NaN;

z = [age, sex_num, geo_dummies];

valid = isfinite(x) & isfinite(y) & all(isfinite(z), 2);
x = x(valid);
y = y(valid);
z = z(valid, :);
age = age(valid);          % <-- add this line, the missing fix

colormap_custom = zeros([length(x), 3]);
colormap_custom(:,1) = (age - min(age)) ./ max(age - min(age));
colormap_custom(:,3) = 1 - age ./ max(age);

figure;
scatter(x, y, 5, colormap_custom, 'filled');
lsline;

[R, P] = corr(x, y)
[rho, p_partial] = partialcorr(x, y, z)

%% HDL vs. MAT
missingCodes = [-8888, -9999, -1111];

x = data.BLD_HDL_COM;
y = data.COG_MAT_SCORE_COM;
age = data.AGE_NMBR_COM;
sex_num = double(categorical(data.SEX_ASK_COM));
geo_cat = categorical(data.GEOSTRATA_COM);
geo_dummies_full = dummyvar(geo_cat);
geo_dummies = geo_dummies_full(:, 2:end);

x(ismember(x, missingCodes)) = NaN;
y(ismember(y, missingCodes)) = NaN;

z = [age, sex_num, geo_dummies];

valid = isfinite(x) & isfinite(y) & all(isfinite(z), 2);
x = x(valid);
y = y(valid);
z = z(valid, :);
age = age(valid);          

colormap_custom = zeros([length(x), 3]);
colormap_custom(:,1) = (age - min(age)) ./ max(age - min(age));
colormap_custom(:,3) = 1 - age ./ max(age);

figure;
scatter(x, y, 5, colormap_custom, 'filled');
lsline;

[R, P] = corr(x, y)
[rho, p_partial] = partialcorr(x, y, z)

%% LDL vs. MAT
missingCodes = [-8888, -9999, -1111];

x = data.BLD_LDL_COM;
y = data.COG_MAT_SCORE_COM;
age = data.AGE_NMBR_COM;
sex_num = double(categorical(data.SEX_ASK_COM));
geo_cat = categorical(data.GEOSTRATA_COM);
geo_dummies_full = dummyvar(geo_cat);
geo_dummies = geo_dummies_full(:, 2:end);

x(ismember(x, missingCodes)) = NaN;
y(ismember(y, missingCodes)) = NaN;

z = [age, sex_num, geo_dummies];

valid = isfinite(x) & isfinite(y) & all(isfinite(z), 2);
x = x(valid);
y = y(valid);
z = z(valid, :);
age = age(valid);          

colormap_custom = zeros([length(x), 3]);
colormap_custom(:,1) = (age - min(age)) ./ max(age - min(age));
colormap_custom(:,3) = 1 - age ./ max(age);

figure;
scatter(x, y, 5, colormap_custom, 'filled');
lsline;


[R, P] = corr(x, y)
[rho, p_partial] = partialcorr(x, y, z)

%% nonHDL vs. MAT
missingCodes = [-8888, -9999, -1111];

x = data.BLD_nonHDL_COM;
y = data.COG_MAT_SCORE_COM;
age = data.AGE_NMBR_COM;
sex_num = double(categorical(data.SEX_ASK_COM));
geo_cat = categorical(data.GEOSTRATA_COM);
geo_dummies_full = dummyvar(geo_cat);
geo_dummies = geo_dummies_full(:, 2:end);

x(ismember(x, missingCodes)) = NaN;
y(ismember(y, missingCodes)) = NaN;

z = [age, sex_num, geo_dummies];

valid = isfinite(x) & isfinite(y) & all(isfinite(z), 2);
x = x(valid);
y = y(valid);
z = z(valid, :);
age = age(valid);          

colormap_custom = zeros([length(x), 3]);
colormap_custom(:,1) = (age - min(age)) ./ max(age - min(age));
colormap_custom(:,3) = 1 - age ./ max(age);

figure;
scatter(x, y, 5, colormap_custom, 'filled');
lsline;

[R, P] = corr(x, y)
[rho, p_partial] = partialcorr(x, y, z)

%% cholesterol vs. MAT
missingCodes = [-8888, -9991, -1111];

x = data.BLD_CHOL_COM;
y = data.COG_MAT_SCORE_COM;
age = data.AGE_NMBR_COM;
sex_num = double(categorical(data.SEX_ASK_COM));
geo_cat = categorical(data.GEOSTRATA_COM);
geo_dummies_full = dummyvar(geo_cat);
geo_dummies = geo_dummies_full(:, 2:end);

x(ismember(x, missingCodes)) = NaN;
y(ismember(y, missingCodes)) = NaN;

z = [age, sex_num, geo_dummies];

valid = isfinite(x) & isfinite(y) & all(isfinite(z), 2);
x = x(valid);
y = y(valid);
z = z(valid, :);
age = age(valid);          

colormap_custom = zeros([length(x), 3]);
colormap_custom(:,1) = (age - min(age)) ./ max(age - min(age));
colormap_custom(:,3) = 1 - age ./ max(age);

figure;
scatter(x, y, 5, colormap_custom, 'filled');
lsline;

[R, P] = corr(x, y)
[rho, p_partial] = partialcorr(x, y, z)

%% IL6 vs REYII
missingCodes = [-8888, -9991, -1111, -2222];

x = data.BLD_IL6_COM;
y = data.COG_REYII_SCORE_COM;
age = data.AGE_NMBR_COM;
sex_num = double(categorical(data.SEX_ASK_COM));
geo_cat = categorical(data.GEOSTRATA_COM);
geo_dummies_full = dummyvar(geo_cat);
geo_dummies = geo_dummies_full(:, 2:end);

x(ismember(x, missingCodes)) = NaN;
y(ismember(y, missingCodes)) = NaN;

z = [age, sex_num, geo_dummies];

valid = isfinite(x) & isfinite(y) & all(isfinite(z), 2);
x = x(valid);
y = y(valid);
z = z(valid, :);
age = age(valid);          

colormap_custom = zeros([length(x), 3]);
colormap_custom(:,1) = (age - min(age)) ./ max(age - min(age));
colormap_custom(:,3) = 1 - age ./ max(age);

subplot(2, 2, 1);
scatter(x, y, 5, colormap_custom, 'filled');
lsline;


[R, P] = corr(x, y)
[rho, p_partial] = partialcorr(x, y, z)

%% TNF vs REYII
tnf_data = readtable("tnf_winsorized_run2.csv");

missingCodes = [-8888, -9991, -1111, -2222];

x = tnf_data.BLD_TNF_COM;
y = data.COG_REYII_SCORE_COM;
age = data.AGE_NMBR_COM;
sex_num = double(categorical(data.SEX_ASK_COM));
geo_cat = categorical(data.GEOSTRATA_COM);
geo_dummies_full = dummyvar(geo_cat);
geo_dummies = geo_dummies_full(:, 2:end);

x(ismember(x, missingCodes)) = NaN;
y(ismember(y, missingCodes)) = NaN;

z = [age, sex_num, geo_dummies];

valid = isfinite(x) & isfinite(y) & all(isfinite(z), 2);
x = x(valid);
y = y(valid);
z = z(valid, :);
age = age(valid);          

colormap_custom = zeros([length(x), 3]);
colormap_custom(:,1) = (age - min(age)) ./ max(age - min(age));
colormap_custom(:,3) = 1 - age ./ max(age);

subplot(2, 2, 2);
scatter(x, y, 5, colormap_custom, 'filled');
lsline;


[R, P] = corr(x, y)
[rho, p_partial] = partialcorr(x, y, z)

%% TNF vs MAT
missingCodes = [-8888, -9991, -1111, -2222];

x = tnf_data.BLD_TNF_COM;
y = data.COG_MAT_SCORE_COM;
age = data.AGE_NMBR_COM;
sex_num = double(categorical(data.SEX_ASK_COM));
geo_cat = categorical(data.GEOSTRATA_COM);
geo_dummies_full = dummyvar(geo_cat);
geo_dummies = geo_dummies_full(:, 2:end);

x(ismember(x, missingCodes)) = NaN;
y(ismember(y, missingCodes)) = NaN;

z = [age, sex_num, geo_dummies];

valid = isfinite(x) & isfinite(y) & all(isfinite(z), 2);
x = x(valid);
y = y(valid);
z = z(valid, :);
age = age(valid);          

colormap_custom = zeros([length(x), 3]);
colormap_custom(:,1) = (age - min(age)) ./ max(age - min(age));
colormap_custom(:,3) = 1 - age ./ max(age);

subplot(2, 2, 4);
scatter(x, y, 5, colormap_custom, 'filled');
lsline;

[R, P] = corr(x, y)
[rho, p_partial] = partialcorr(x, y, z)

%% IL6 v MAT
missingCodes = [-8888, -9991, -1111, -2222];

x = data.BLD_IL6_COM;
y = data.COG_MAT_SCORE_COM;
age = data.AGE_NMBR_COM;
sex_num = double(categorical(data.SEX_ASK_COM));
geo_cat = categorical(data.GEOSTRATA_COM);
geo_dummies_full = dummyvar(geo_cat);
geo_dummies = geo_dummies_full(:, 2:end);

x(ismember(x, missingCodes)) = NaN;
y(ismember(y, missingCodes)) = NaN;

z = [age, sex_num, geo_dummies];

valid = isfinite(x) & isfinite(y) & all(isfinite(z), 2);
x = x(valid);
y = y(valid);
z = z(valid, :);
age = age(valid);          

colormap_custom = zeros([length(x), 3]);
colormap_custom(:,1) = (age - min(age)) ./ max(age - min(age));
colormap_custom(:,3) = 1 - age ./ max(age);

subplot(2, 2, 3);
scatter(x, y, 5, colormap_custom, 'filled');
lsline;


[R, P] = corr(x, y)
[rho, p_partial] = partialcorr(x, y, z)