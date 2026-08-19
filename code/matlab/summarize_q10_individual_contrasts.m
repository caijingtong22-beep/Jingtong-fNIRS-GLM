function profiles = summarize_q10_individual_contrasts()
% Create a compact per-subject Q10 HILoadHE minus HILoadLE profile.

this_file = mfilename('fullpath');
model_root = fileparts(this_file);
q10_root = fullfile(model_root, 'q10_hiloadhe_vs_hiloadle');
cohort_root = fullfile(q10_root, 'cohort_specification');
profiles = table();

for subject_num = 1:30
    subject_label = sprintf('subject%03d', subject_num);
    for hb_label = ["HbO", "HbR"]
        if subject_num == 1
            result_file = fullfile(q10_root, subject_label, ...
                sprintf('%s_%s_q10_HE_minus_LE_channel_contrasts.csv', subject_label, hb_label));
        else
            result_file = fullfile(cohort_root, subject_label, ...
                sprintf('%s_%s_q10_HE_minus_LE_channel_contrasts.csv', subject_label, hb_label));
        end
        assert(isfile(result_file), 'Missing contrast result: %s', result_file);
        T = readtable(result_file, 'TextType', 'string');
        p = T.contrast_p_two_sided;
        beta = T.contrast_beta;
        q = bh_fdr(p);
        is_long = strcmpi(T.channel_type, 'Long');
        profile = table(string(subject_label), hb_label, height(T), sum(is_long), ...
            sum(p < 0.05), sum(p < 0.05 & beta > 0), sum(p < 0.05 & beta < 0), ...
            sum(q < 0.05), sum(q < 0.05 & beta > 0), sum(q < 0.05 & beta < 0), ...
            sum(q < 0.05 & is_long), sum(q < 0.05 & is_long & beta > 0), ...
            sum(q < 0.05 & is_long & beta < 0), median(beta, 'omitnan'), ...
            median(abs(beta), 'omitnan'), ...
            'VariableNames', {'subject_id','hemoglobin','n_channels','n_long_channels', ...
            'n_p_lt_005','n_positive_p_lt_005','n_negative_p_lt_005', ...
            'n_q_bh61_lt_005','n_positive_q_bh61_lt_005','n_negative_q_bh61_lt_005', ...
            'n_long_q_bh61_lt_005','n_positive_long_q_bh61_lt_005', ...
            'n_negative_long_q_bh61_lt_005','median_contrast_beta', ...
            'median_abs_contrast_beta'});
        profiles = [profiles; profile]; %#ok<AGROW>
    end
end

writetable(profiles, fullfile(q10_root, 'q10_HE_minus_LE_individual_profile_n30.csv'));
fprintf('Q10 individual profiles complete: 30 subjects, HbO and HbR.\n');
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
