function [individual_roi, group_roi, channel_long] = summarize_q10_roi_long45_allqc()
% Q10 ROI analysis using all 45 physical long channels, regardless of QC flag.
% Contrast is HILoadHE minus HILoadLE. This is an exploratory analysis.

this_file = mfilename('fullpath');
model_root = fileparts(this_file);
q10_root = fullfile(model_root, 'q10_hiloadhe_vs_hiloadle');
cohort_root = fullfile(q10_root, 'cohort_specification');
subject_numbers = 1:30;
individual_roi = table();
channel_long = table();

for hb_label = ["HbO", "HbR"]
    beta_by_subject = [];
    p_by_subject = [];
    channel_info = table();
    for subject_num = subject_numbers
        subject_label = sprintf('subject%03d', subject_num);
        T = read_subject_contrast(q10_root, cohort_root, subject_num, hb_label);
        long_rows = strcmpi(T.channel_type, 'Long');
        T = T(long_rows, :);
        assert(height(T) == 45, '%s %s has %d long channels, expected 45.', ...
            subject_label, hb_label, height(T));

        if isempty(channel_info)
            channel_info = T(:, {'channel_index','channel_id','roi','laterality'});
        else
            assert(isequal(channel_info.channel_index, T.channel_index) && ...
                isequal(channel_info.channel_id, T.channel_id), ...
                'Long-channel identity mismatch for %s.', subject_label);
        end
        beta_by_subject(subject_num, :) = T.contrast_beta(:)'; %#ok<AGROW>
        p_by_subject(subject_num, :) = T.contrast_p_two_sided(:)'; %#ok<AGROW>

        group_keys = unique(T(:, {'roi','laterality'}), 'rows', 'stable');
        for roi_idx = 1:height(group_keys)
            in_roi = T.roi == group_keys.roi(roi_idx) & ...
                T.laterality == group_keys.laterality(roi_idx);
            beta = T.contrast_beta(in_roi);
            p = T.contrast_p_two_sided(in_roi);
            individual_roi = [individual_roi; table(string(subject_label), hb_label, ...
                group_keys.roi(roi_idx), group_keys.laterality(roi_idx), sum(in_roi), ...
                mean(beta, 'omitnan'), median(beta, 'omitnan'), ...
                sum(beta > 0), sum(beta < 0), sum(p < 0.05), ...
                'VariableNames', {'subject_id','hemoglobin','roi','laterality', ...
                'n_long_channels','mean_contrast_beta','median_contrast_beta', ...
                'n_positive_beta','n_negative_beta','n_channel_p_lt_005'})]; %#ok<AGROW>
        end
    end

    % Third layer: long-channel group result, with BH-FDR across 45 long channels.
    n = sum(isfinite(beta_by_subject), 1);
    mean_beta = mean(beta_by_subject, 1, 'omitnan');
    sd_beta = std(beta_by_subject, 0, 1, 'omitnan');
    se_beta = sd_beta ./ sqrt(n);
    t_beta = mean_beta ./ se_beta;
    t_beta(n < 2 | se_beta == 0) = NaN;
    p_beta = two_sided_t_p(t_beta, n - 1);
    q_beta = bh_fdr(p_beta);
    channel_long = [channel_long; table(repmat(hb_label, 45, 1), n(:), ...
        mean_beta(:), sd_beta(:), se_beta(:), t_beta(:), p_beta(:), q_beta(:), ...
        sum(beta_by_subject > 0, 1)', sum(beta_by_subject < 0, 1)', ...
        sum(p_by_subject < 0.05, 1)', ...
        'VariableNames', {'hemoglobin','n_subjects','mean_contrast_beta', ...
        'sd_contrast_beta','group_se','group_t','group_p_two_sided', ...
        'group_q_bh45','n_positive_beta','n_negative_beta', ...
        'n_individual_p_lt_005'}), channel_info]; %#ok<AGROW>
end

% Second layer: each subject contributes one ROI mean beta, then group test.
group_roi = table();
for hb_label = ["HbO", "HbR"]
    U = individual_roi(individual_roi.hemoglobin == hb_label, :);
    keys = unique(U(:, {'roi','laterality'}), 'rows', 'stable');
    this_hb = table();
    for roi_idx = 1:height(keys)
        rows = U.roi == keys.roi(roi_idx) & U.laterality == keys.laterality(roi_idx);
        values = U.mean_contrast_beta(rows);
        n = sum(isfinite(values));
        m = mean(values, 'omitnan');
        sd = std(values, 0, 'omitnan');
        se = sd / sqrt(n);
        t = m / se;
        if n < 2 || se == 0
            t = NaN;
        end
        p = two_sided_t_p(t, n - 1);
        n_channels_per_subject = mean(U.n_long_channels(rows), 'omitnan');
        this_hb = [this_hb; table(hb_label, keys.roi(roi_idx), ...
            keys.laterality(roi_idx), n, n_channels_per_subject, m, sd, se, t, p, ...
            'VariableNames', {'hemoglobin','roi','laterality','n_subjects', ...
            'n_long_channels_per_subject','mean_roi_contrast_beta', ...
            'sd_roi_contrast_beta','group_se','group_t','group_p_two_sided'})]; %#ok<AGROW>
    end
    this_hb.group_q_bh_roi_cells = bh_fdr(this_hb.group_p_two_sided);
    group_roi = [group_roi; this_hb]; %#ok<AGROW>
end

writetable(individual_roi, fullfile(q10_root, 'q10_HE_minus_LE_individual_ROI_long45_allQC_n30.csv'));
writetable(group_roi, fullfile(q10_root, 'q10_HE_minus_LE_group_ROI_long45_allQC_n30.csv'));
writetable(channel_long, fullfile(q10_root, 'q10_HE_minus_LE_group_channel_long45_allQC_n30.csv'));
fprintf('Q10 ROI analysis complete: 30 subjects, all 45 long channels, QC flags retained.\n');
fprintf('Group ROI cells with BH-FDR q<.05: HbO=%d, HbR=%d\n', ...
    sum(group_roi.group_q_bh_roi_cells(group_roi.hemoglobin == "HbO") < 0.05), ...
    sum(group_roi.group_q_bh_roi_cells(group_roi.hemoglobin == "HbR") < 0.05));
end

function T = read_subject_contrast(q10_root, cohort_root, subject_num, hb_label)
subject_label = sprintf('subject%03d', subject_num);
if subject_num == 1
    result_file = fullfile(q10_root, subject_label, ...
        sprintf('%s_%s_q10_HE_minus_LE_channel_contrasts.csv', subject_label, hb_label));
else
    result_file = fullfile(cohort_root, subject_label, ...
        sprintf('%s_%s_q10_HE_minus_LE_channel_contrasts.csv', subject_label, hb_label));
end
assert(isfile(result_file), 'Missing contrast result: %s', result_file);
T = readtable(result_file, 'TextType', 'string');
end

function p = two_sided_t_p(t, df)
p = betainc(df ./ (df + t.^2), df/2, 1/2);
p(~isfinite(t) | df <= 0) = NaN;
end

function q = bh_fdr(p)
q = nan(size(p));
valid = isfinite(p);
pv = p(valid);
if isempty(pv), return; end
[sorted, order] = sort(pv(:));
m = numel(sorted);
q_sorted = sorted .* m ./ (1:m)';
q_sorted = min(1, flipud(cummin(flipud(q_sorted))));
qv = zeros(size(pv));
qv(order) = q_sorted;
q(valid) = qv;
end
