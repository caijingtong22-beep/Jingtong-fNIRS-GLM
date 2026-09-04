function out = run_q9_hiload_group_analysis()
%RUN_Q9_HILOAD_GROUP_ANALYSIS Group tests for direct HILoad - ExperimentalBaseline.
% Reads only the new 29-subject direct-HILoad first-level contrast outputs.
% Q9 is intentionally not split into Procedure, DOM, and AOE.

project = 'C:\aliceproject\Jingtong-fNIRS-GLM';
cohort = fullfile(project,'results','q9_hiload_vs_experimentalbaseline','cohort_specification');
outdir = fullfile(cohort,'group_analysis_n29');
if ~isfolder(outdir), mkdir(outdir); end
hbs = ["HbO","HbR"];

allrows = table();
for hb = hbs
    ff = dir(fullfile(cohort,'subject*',sprintf('*_%s_q9_hiload_minus_experimentalbaseline_channel_contrasts.csv',char(hb))));
    assert(numel(ff)==29,'Expected 29 %s contrast files, found %d.',hb,numel(ff));
    for i = 1:numel(ff)
        T = readtable(fullfile(ff(i).folder,ff(i).name),'TextType','string');
        assert(height(T)==61,'Expected 61 channels in %s.',ff(i).name);
        allrows = [allrows; T]; %#ok<AGROW>
    end
end
allrows.subject_id = string(allrows.subject_id);
allrows.hemoglobin = string(allrows.hemoglobin);
allrows.channel_type = string(allrows.channel_type);
allrows.roi = string(allrows.roi);
allrows.laterality = string(allrows.laterality);
subjects = unique(allrows.subject_id,'stable');
assert(numel(subjects)==29,'Expected 29 unique eligible subjects.');

% Channel-level tests: 61-channel BH family separately for HbO and HbR.
channel_summary = table();
for hb = hbs
    H = allrows(allrows.hemoglobin==hb,:);
    key = H(1:61,{'channel_index','channel_id','channel_type','roi','laterality'});
    key = sortrows(key,'channel_index');
    pvals = nan(height(key),1); rows = table();
    for c = 1:height(key)
        x = H.contrast_beta(H.channel_index==key.channel_index(c));
        [t,p] = one_sample_t(x); pvals(c)=p;
        rows = [rows; table(hb,key.channel_index(c),key.channel_id(c),key.channel_type(c), ...
            key.roi(c),key.laterality(c),numel(x),mean(x),std(x),std(x)/sqrt(numel(x)),t,p, ...
            sum(x>0),sum(x<0),'VariableNames',{'hemoglobin','channel_index','channel_id', ...
            'channel_type','roi','laterality','n_subjects','mean_contrast_beta','sd_contrast_beta', ...
            'se_contrast_beta','t_statistic','p_value','n_positive_beta','n_negative_beta'})]; %#ok<AGROW>
    end
    rows.q_value_bh_fdr = bh_fdr(pvals);
    channel_summary = [channel_summary; rows]; %#ok<AGROW>
end
writetable(channel_summary,fullfile(outdir,'q9_hiload_minus_experimentalbaseline_group_channel_summary_61ch_n29.csv'));

% ROI means: ordinary within-subject mean across allocation-defined long channels.
long = allrows(allrows.channel_type=="Long",:);
roi_keys = unique(long(long.hemoglobin=="HbO",{'roi','laterality'}),'rows','stable');
individual_roi = table();
for s = 1:numel(subjects)
    for hb = hbs
        H = long(long.subject_id==subjects(s) & long.hemoglobin==hb,:);
        assert(height(H)==45,'Subject %s / %s is missing long channels.',subjects(s),hb);
        for r = 1:height(roi_keys)
            x = H.contrast_beta(H.roi==roi_keys.roi(r) & H.laterality==roi_keys.laterality(r));
            individual_roi = [individual_roi; table(subjects(s),hb,roi_keys.roi(r),roi_keys.laterality(r), ...
                numel(x),mean(x),median(x),'VariableNames',{'subject_id','hemoglobin','roi','laterality', ...
                'n_long_channels','mean_contrast_beta','median_contrast_beta'})]; %#ok<AGROW>
        end
    end
end
writetable(individual_roi,fullfile(outdir,'q9_hiload_minus_experimentalbaseline_individual_ROI_long45_n29.csv'));

roi_summary = table();
for hb = hbs
    H = individual_roi(individual_roi.hemoglobin==hb,:); pvals = nan(height(roi_keys),1); rows=table();
    for r = 1:height(roi_keys)
        x = H.mean_contrast_beta(H.roi==roi_keys.roi(r) & H.laterality==roi_keys.laterality(r));
        [t,p]=one_sample_t(x); pvals(r)=p;
        rows=[rows;table(hb,roi_keys.roi(r),roi_keys.laterality(r),numel(x),mean(x),std(x), ...
            std(x)/sqrt(numel(x)),t,p,sum(x>0),sum(x<0),'VariableNames',{'hemoglobin','roi','laterality', ...
            'n_subjects','mean_contrast_beta','sd_contrast_beta','se_contrast_beta','t_statistic','p_value', ...
            'n_positive_beta','n_negative_beta'})]; %#ok<AGROW>
    end
    rows.q_value_bh_fdr=bh_fdr(pvals); roi_summary=[roi_summary;rows]; %#ok<AGROW>
end
writetable(roi_summary,fullfile(outdir,'q9_hiload_minus_experimentalbaseline_group_ROI_long45_n29.csv'));

save(fullfile(outdir,'q9_hiload_minus_experimentalbaseline_group_analysis_n29.mat'), ...
    'channel_summary','individual_roi','roi_summary','subjects');
plot_q9_hiload_topography(channel_summary,project,outdir);

out=struct('output_dir',outdir,'n_subjects',numel(subjects), ...
    'channel_summary',channel_summary,'roi_summary',roi_summary);
fprintf('Q9 direct-HILoad group analysis complete: n=%d; %s\n',numel(subjects),outdir);
end

function [t,p] = one_sample_t(x)
x=x(isfinite(x)); n=numel(x);
if n<2 || std(x)==0, t=NaN; p=NaN; return; end
t=mean(x)/(std(x)/sqrt(n)); p=betainc(n/(n+t^2),n/2,.5);
end

function q = bh_fdr(p)
q=nan(size(p)); valid=isfinite(p); pv=p(valid); [sp,order]=sort(pv(:)); m=numel(sp);
adj=min(1,flipud(cummin(flipud(sp.*m./(1:m)')))); restored=zeros(size(pv)); restored(order)=adj;
q(valid)=restored;
end
