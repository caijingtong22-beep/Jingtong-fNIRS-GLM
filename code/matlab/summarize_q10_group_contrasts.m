function [group_results, overview] = summarize_q10_group_contrasts()
% Summarize Q10 HILoadHE minus HILoadLE contrasts across 30 subjects.
% Group inference is an unweighted one-sample t-test of subject-level beta
% contrasts against zero, followed by BH-FDR across the 61 channels per Hb.

this_file = mfilename('fullpath');
model_root = fileparts(this_file);
q10_root = fullfile(model_root, 'q10_hiloadhe_vs_hiloadle');
cohort_root = fullfile(q10_root, 'cohort_specification');
subject_numbers = 1:30;
group_results = table();
overview = table();

for hb_label = ["HbO", "HbR"]
    beta_matrix = [];
    p_matrix = [];
    channel_info = table();
    for subject_num = subject_numbers
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
        assert(height(T) == 61, '%s does not contain 61 channels.', result_file);
        if isempty(channel_info)
            channel_info = T(:, {'channel_index','channel_id','channel_type','roi','laterality','qc_status'});
        else
            assert(isequal(channel_info.channel_index, T.channel_index) && ...
                isequal(channel_info.channel_id, T.channel_id), ...
                'Channel identity mismatch in %s', result_file);
        end
        beta_matrix(subject_num, :) = T.contrast_beta(:)'; %#ok<AGROW>
        p_matrix(subject_num, :) = T.contrast_p_two_sided(:)'; %#ok<AGROW>
    end

    n = sum(isfinite(beta_matrix), 1);
    mean_beta = mean(beta_matrix, 1, 'omitnan');
    sd_beta = std(beta_matrix, 0, 1, 'omitnan');
    se_beta = sd_beta ./ sqrt(n);
    group_t = mean_beta ./ se_beta;
    group_t(n < 2 | se_beta == 0) = NaN;
    group_p = two_sided_t_p(group_t, n - 1);
    group_q = bh_fdr(group_p);
    n_positive = sum(beta_matrix > 0, 1);
    n_negative = sum(beta_matrix < 0, 1);
    n_individual_p = sum(p_matrix < 0.05, 1);
    n_positive_individual_p = sum((p_matrix < 0.05) & (beta_matrix > 0), 1);
    n_negative_individual_p = sum((p_matrix < 0.05) & (beta_matrix < 0), 1);

    Tgroup = [table(repmat(hb_label, 61, 1), n(:), mean_beta(:), sd_beta(:), ...
        se_beta(:), group_t(:), group_p(:), group_q(:), n_positive(:), ...
        n_negative(:), n_individual_p(:), n_positive_individual_p(:), ...
        n_negative_individual_p(:), ...
        'VariableNames', {'hemoglobin','n_subjects','mean_contrast_beta', ...
        'sd_contrast_beta','group_se','group_t','group_p_two_sided', ...
        'group_q_bh61','n_positive_beta','n_negative_beta', ...
        'n_individual_p_lt_005','n_positive_individual_p_lt_005', ...
        'n_negative_individual_p_lt_005'}), channel_info];
    group_results = [group_results; Tgroup]; %#ok<AGROW>

    fdr_sig = group_q < 0.05;
    overview = [overview; table(hb_label, sum(fdr_sig), sum(fdr_sig & mean_beta > 0), ...
        sum(fdr_sig & mean_beta < 0), median(max(n_positive, n_negative) ./ n), ...
        median(n_individual_p), ...
        'VariableNames', {'hemoglobin','n_group_q_bh61_lt_005', ...
        'n_group_positive_q_bh61_lt_005','n_group_negative_q_bh61_lt_005', ...
        'median_directional_agreement','median_n_subjects_individual_p_lt_005'})]; %#ok<AGROW>
end

writetable(group_results, fullfile(q10_root, 'q10_HE_minus_LE_group_channel_summary_n30.csv'));
writetable(overview, fullfile(q10_root, 'q10_HE_minus_LE_group_overview_n30.csv'));
fprintf('Q10 group summary complete: 30 subjects, 61 channels, HbO and HbR.\n');
disp(overview);
end

function p = two_sided_t_p(t, df)
% Equivalent to 2*tcdf(-abs(t), df), using base MATLAB betainc.
p = betainc(df ./ (df + t.^2), df/2, 1/2);
p(~isfinite(t) | df <= 0) = NaN;
end

function q = bh_fdr(p)
q = nan(size(p));
valid = isfinite(p);
pv = p(valid);
if isempty(pv)
    return;
end
[sorted, order] = sort(pv(:));
m = numel(sorted);
q_sorted = sorted .* m ./ (1:m)';
q_sorted = min(1, flipud(cummin(flipud(q_sorted))));
qv = zeros(size(pv));
qv(order) = q_sorted;
q(valid) = qv;
end
