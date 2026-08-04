clear
pro = load('clinical_prot_elnet_schwarz.mat');
meta = load('meta_elnet_schwarz.mat');
pro_meta = load('pro_met_elnet_schwarz.mat');

figure;

%Proteomics
%clinical
y=mean(pro.r_x_split);
err=std(pro.r_x_split);
errorbar(0.5, y, err,"square",'Color','k');
xlim([0 7.5]);
%clinical + proteomics
y=mean(pro.r_x_split_prot);
err=std(pro.r_x_split_prot);
hold on; errorbar(1.5, y, err,"square", 'Color','b');
xlim([0 7.5]);
ylim([0 0.6])
%clinical + proteomics non CV
y=pro.r_yhat;
hold on; errorbar(2, y, 0,"square", 'Color',[0.4 0.5 1]);
xlim([0 7.5]);
ylim([0 0.6])

%Metabolomics
%clinical
y=mean(meta.r_x_split_clin);
err=std(meta.r_x_split_clin);
errorbar(3, y, err,"square",'Color','k');
xlim([0 7.5]);
%clinical + metabolomics
y=mean(meta.r_x_split_meta);
err=std(meta.r_x_split_meta);
hold on; errorbar(4, y, err, "square", 'Color', [0 0.65 0]);
xlim([0 7.5]);
ylim([0 0.6])
%clinical + metabolomics non CV
y=meta.r_yhat;
hold on; errorbar(4.5, y, 0,"square", 'Color', [0.25 0.80 0.25]);
xlim([0 7.5]);
ylim([0 0.6])

%Proteomics + Metabolomics
%clinical
y=mean(pro_meta.r_x_split_clin);
err=std(pro_meta.r_x_split_clin);
errorbar(5.5, y, err,"square",'Color','k');
xlim([0 7.5]);
%clinical + proteomics + metabolomics
y=mean(pro_meta.r_x_split_pro_met);
err=std(pro_meta.r_x_split_pro_met);
hold on; errorbar(6.5, y, err, "square", 'Color', [0 0.5 0.5]);
xlim([0 7.5]);
ylim([0 0.6])
%clinical + proteomics + metabolomics non CV
y=pro_meta.r_yhat;
hold on; errorbar(7, y, 0,"square", 'Color', [0.4 0.7 0.7]);
xlim([0 7.5]);
ylim([0 0.6])

%pretty aesthetics
set(gca, 'XTick', []);
set(gca, 'YAxisLocation', 'left'); % ensure only left axis exists
set(gca, 'box', 'off'); % removes the right axis line + ticks
% Set y-axis ticks and labels
yticks(0:0.1:0.6);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6'});
% Add light dotted horizontal reference lines
hold on
yvals = 0.1:0.1:0.5;
for yline_val = yvals
    yline(yline_val, 'k:', 'LineWidth', 0.5);
end

%% he
savefig('schwarz_figure.fig');

%% for subplot
save('schwarz_figures.mat', figure);