function [individual_roi, group_roi] = summarize_q10_roi_cortical_system_long45_allqc()
% Q10 ROI summaries after pooling left/right channels into cortical systems.
% Uses all 45 long channels and retains every QC status.

this_file = mfilename('fullpath');
model_root = fileparts(this_file);
q10_root = fullfile(model_root, 'q10_hiloadhe_vs_hiloadle');
cohort_root = fullfile(q10_root, 'cohort_specification');
roi_order = ["Medial Frontopolar PFC", "Lateral Frontopolar PFC", "DL PFC", ...
    "VL PFC", "Medial PFC", "PMC", "SMA", "M1"];
system_order = ["Prefrontal Cortex", "Prefrontal Cortex", "Prefrontal Cortex", ...
    "Prefrontal Cortex", "Prefrontal Cortex", "Motor Cortex", "Motor Cortex", "Motor Cortex"];
individual_roi = table();

for subject_num = 1:30
    subject_label = sprintf('subject%03d', subject_num);
    for hb_label = ["HbO", "HbR"]
        T = read_subject_contrast(q10_root, cohort_root, subject_num, hb_label);
        T = T(strcmpi(T.channel_type, 'Long'), :);
        assert(height(T) == 45, '%s %s does not contain 45 long channels.', subject_label, hb_label);
        for roi_idx = 1:numel(roi_order)
            in_roi = T.roi == roi_order(roi_idx);
            beta = T.contrast_beta(in_roi);
            individual_roi = [individual_roi; table(string(subject_label), hb_label, ...
                system_order(roi_idx), roi_order(roi_idx), sum(in_roi), ...
                mean(beta, 'omitnan'), median(beta, 'omitnan'), ...
                sum(beta > 0), sum(beta < 0), ...
                'VariableNames', {'subject_id','hemoglobin','cortical_system','roi', ...
                'n_long_channels','mean_contrast_beta','median_contrast_beta', ...
                'n_positive_beta','n_negative_beta'})]; %#ok<AGROW>
        end
    end
end

group_roi = table();
for hb_label = ["HbO", "HbR"]
    U = individual_roi(individual_roi.hemoglobin == hb_label, :);
    this_hb = table();
    for roi_idx = 1:numel(roi_order)
        V = U(U.roi == roi_order(roi_idx), :);
        beta = V.mean_contrast_beta;
        n = sum(isfinite(beta));
        mean_beta = mean(beta, 'omitnan');
        sd_beta = std(beta, 0, 'omitnan');
        se_beta = sd_beta / sqrt(n);
        t = mean_beta / se_beta;
        if n < 2 || se_beta == 0, t = NaN; end
        p = two_sided_t_p(t, n - 1);
        this_hb = [this_hb; table(hb_label, system_order(roi_idx), roi_order(roi_idx), ...
            n, mean(V.n_long_channels, 'omitnan'), mean_beta, sd_beta, se_beta, t, p, ...
            'VariableNames', {'hemoglobin','cortical_system','roi','n_subjects', ...
            'n_long_channels_per_subject','mean_roi_contrast_beta', ...
            'sd_roi_contrast_beta','group_se','group_t','group_p_two_sided'})]; %#ok<AGROW>
    end
    this_hb.group_q_bh8 = bh_fdr(this_hb.group_p_two_sided);
    group_roi = [group_roi; this_hb]; %#ok<AGROW>
end

writetable(individual_roi, fullfile(q10_root, 'q10_HE_minus_LE_individual_ROI_cortical_system_long45_allQC_n30.csv'));
writetable(group_roi, fullfile(q10_root, 'q10_HE_minus_LE_group_ROI_cortical_system_long45_allQC_n30.csv'));
fprintf('Q10 cortical-system ROI summary complete: 8 ROI groups, 30 subjects.\n');
end

function T = read_subject_contrast(q10_root, cohort_root, subject_num, hb_label)
subject_label = sprintf('subject%03d', subject_num);
if subject_num == 1
    f = fullfile(q10_root, subject_label, sprintf('%s_%s_q10_HE_minus_LE_channel_contrasts.csv', subject_label, hb_label));
else
    f = fullfile(cohort_root, subject_label, sprintf('%s_%s_q10_HE_minus_LE_channel_contrasts.csv', subject_label, hb_label));
end
T = readtable(f, 'TextType', 'string');
end

function p = two_sided_t_p(t, df)
p = betainc(df ./ (df + t.^2), df/2, 1/2);
p(~isfinite(t) | df <= 0) = NaN;
end

function q = bh_fdr(p)
q = nan(size(p)); valid = isfinite(p); pv = p(valid);
if isempty(pv), return; end
[s, order] = sort(pv(:)); m = numel(s);
qs = min(1, flipud(cummin(flipud(s .* m ./ (1:m)'))));
qv = zeros(size(pv)); qv(order) = qs; q(valid) = qv;
end
