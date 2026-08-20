function summarize_q8_procedure_roi_long45_allqc()
% Q8 Procedure: HIHE - LILE long-channel (45) ROI summaries, all QC retained.
project_root = 'C:\aliceproject\Jingtong-fNIRS-GLM';
cohort_root = fullfile(project_root, 'results', 'q8_hihe_vs_lile_procedure', 'cohort_specification');
subject_ids = 1:30; hb_labels = {'HbO','HbR'};
individual_rows = table(); channel_rows = table();
for subject_id = subject_ids
    subject_label = sprintf('subject%03d', subject_id);
    for h = 1:numel(hb_labels)
        hb = hb_labels{h};
        file_path = fullfile(cohort_root, subject_label, sprintf('%s_%s_q8_HIHE_minus_LILE_procedure_channel_contrasts.csv', subject_label, hb));
        assert(isfile(file_path), 'Missing contrast file: %s', file_path);
        T = readtable(file_path, 'TextType', 'string'); T = T(T.channel_type == "Long", :);
        assert(height(T) == 45, '%s %s does not contain exactly 45 long channels.', subject_label, hb);
        channel_rows = [channel_rows; T]; %#ok<AGROW>
        keys = unique(T(:, {'roi','laterality'}), 'rows', 'stable');
        for k = 1:height(keys)
            idx = T.roi == keys.roi(k) & T.laterality == keys.laterality(k);
            row = table(subject_id, string(hb), keys.roi(k), keys.laterality(k), sum(idx), mean(T.contrast_beta(idx), 'omitnan'), median(T.contrast_beta(idx), 'omitnan'), sum(T.contrast_beta(idx)>0), sum(T.contrast_beta(idx)<0), sum(T.contrast_p_two_sided(idx)<.05), 'VariableNames', {'subject_id','hemoglobin','roi','laterality','n_long_channels','mean_contrast_beta','median_contrast_beta','n_positive_beta','n_negative_beta','n_channel_p_lt_005'});
            individual_rows = [individual_rows; row]; %#ok<AGROW>
        end
    end
end
writetable(individual_rows, fullfile(cohort_root, 'q8_HIHE_minus_LILE_procedure_individual_ROI_long45_allQC_n30.csv'));

channel_group = table();
for h = 1:numel(hb_labels)
    H = channel_rows(channel_rows.hemoglobin == string(hb_labels{h}), :); ids = unique(H.channel_id); rows = table(); pvals = nan(numel(ids),1);
    for c = 1:numel(ids)
        C = H(H.channel_id == ids(c), :); x = C.contrast_beta; x = x(isfinite(x)); [ts,p] = one_sample_t(x); pvals(c)=p;
        rows = [rows; table(string(hb_labels{h}), ids(c), C.channel_type(1), C.roi(1), C.laterality(1), numel(x), mean(x), std(x), std(x)/sqrt(numel(x)), ts, p, sum(x>0), sum(x<0), 'VariableNames', {'hemoglobin','channel','channel_type','roi','laterality','n_subjects','mean_contrast_beta','sd_contrast_beta','se_contrast_beta','t_statistic','p_value','n_positive_beta','n_negative_beta'})]; %#ok<AGROW>
    end
    rows.q_value_bh_fdr = bh_fdr(pvals); channel_group = [channel_group; rows]; %#ok<AGROW>
end
writetable(channel_group, fullfile(cohort_root, 'q8_HIHE_minus_LILE_procedure_group_channel_long45_allQC_n30.csv'));

roi_group = table();
for h = 1:numel(hb_labels)
    H = individual_rows(individual_rows.hemoglobin == string(hb_labels{h}), :); keys = unique(H(:, {'roi','laterality'}), 'rows', 'stable'); rows = table(); pvals = nan(height(keys),1);
    for k = 1:height(keys)
        R = H(H.roi == keys.roi(k) & H.laterality == keys.laterality(k), :); x = R.mean_contrast_beta; x=x(isfinite(x)); [ts,p] = one_sample_t(x); pvals(k)=p;
        rows = [rows; table(string(hb_labels{h}), keys.roi(k), keys.laterality(k), numel(x), mean(x), std(x), std(x)/sqrt(numel(x)), ts, p, sum(x>0), sum(x<0), 'VariableNames', {'hemoglobin','roi','laterality','n_subjects','mean_contrast_beta','sd_contrast_beta','se_contrast_beta','t_statistic','p_value','n_positive_beta','n_negative_beta'})]; %#ok<AGROW>
    end
    rows.q_value_bh_fdr = bh_fdr(pvals); roi_group = [roi_group; rows]; %#ok<AGROW>
end
writetable(roi_group, fullfile(cohort_root, 'q8_HIHE_minus_LILE_procedure_group_ROI_long45_allQC_n30.csv'));
fprintf('Q8 Procedure long-45 ROI summaries written to %s\n', cohort_root);
end
function [t,p] = one_sample_t(x)
n=numel(x); if n<2 || std(x)==0, t=NaN; p=NaN; return; end; t=mean(x)/(std(x)/sqrt(n)); p=betainc(n/(n+t^2),n/2,.5);
end
function q=bh_fdr(p)
q=nan(size(p)); v=isfinite(p); [sp,o]=sort(p(v)); m=numel(sp); a=sp.*m./(1:m)'; a=flipud(cummin(flipud(a))); ii=find(v); q(ii(o))=min(a,1);
end
