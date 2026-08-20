function [group_results, overview] = summarize_q8_procedure_group_contrasts()
% Unweighted n=30 one-sample summary of Q8 Procedure subject-level betas.
% This is a descriptive group analysis, not a formal NIRS-SPM second level.

this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(fileparts(this_file)));
cohort_root = fullfile(repo_root, 'results', 'q8_hihe_vs_lile_procedure', 'cohort_specification');
subject_numbers = 1:30; group_results = table(); overview = table();
for hb_label = ["HbO", "HbR"]
    beta_matrix = []; channel_info = table();
    for subject_num = subject_numbers
        subject_label = sprintf('subject%03d', subject_num);
        result_file = fullfile(cohort_root, subject_label, sprintf('%s_%s_q8_HIHE_minus_LILE_procedure_channel_contrasts.csv', subject_label, hb_label));
        assert(isfile(result_file), 'Missing contrast result: %s', result_file);
        T = readtable(result_file, 'TextType', 'string'); assert(height(T) == 61, 'Expected 61 channels.');
        if isempty(channel_info), channel_info = T(:, {'channel_index','channel_id','channel_type','roi','laterality','qc_status'});
        else, assert(isequal(channel_info.channel_index,T.channel_index) && isequal(channel_info.channel_id,T.channel_id), 'Channel identity mismatch.'); end
        beta_matrix(subject_num,:) = T.contrast_beta(:)'; %#ok<AGROW>
    end
    n = sum(isfinite(beta_matrix),1); mean_beta = mean(beta_matrix,1,'omitnan'); sd_beta = std(beta_matrix,0,1,'omitnan');
    se_beta = sd_beta ./ sqrt(n); group_t = mean_beta ./ se_beta; group_t(n < 2 | se_beta == 0) = NaN;
    group_p = two_sided_t_p(group_t,n-1); group_q = bh_fdr(group_p);
    n_positive = sum(beta_matrix > 0,1); n_negative = sum(beta_matrix < 0,1);
    Tgroup = [table(repmat(hb_label,61,1),n(:),mean_beta(:),sd_beta(:),se_beta(:),group_t(:),group_p(:),group_q(:),n_positive(:),n_negative(:), ...
        'VariableNames', {'hemoglobin','n_subjects','mean_contrast_beta','sd_contrast_beta','group_se','group_t','group_p_two_sided','group_q_bh61','n_positive_beta','n_negative_beta'}), channel_info];
    group_results = [group_results; Tgroup]; %#ok<AGROW>
    fdr_sig = group_q < .05;
    overview = [overview; table(hb_label,sum(fdr_sig),sum(fdr_sig & mean_beta > 0),sum(fdr_sig & mean_beta < 0),median(max(n_positive,n_negative)./n), ...
        'VariableNames', {'hemoglobin','n_group_q_bh61_lt_005','n_group_positive_q_bh61_lt_005','n_group_negative_q_bh61_lt_005','median_directional_agreement'})]; %#ok<AGROW>
end
writetable(group_results,fullfile(cohort_root,'q8_HIHE_minus_LILE_procedure_group_channel_summary_n30.csv'));
writetable(overview,fullfile(cohort_root,'q8_HIHE_minus_LILE_procedure_group_overview_n30.csv'));
disp(overview);
end

function p = two_sided_t_p(t,df)
p = betainc(df ./ (df + t.^2),df/2,1/2); p(~isfinite(t) | df <= 0) = NaN;
end

function q = bh_fdr(p)
q = nan(size(p)); valid = isfinite(p); pv = p(valid); if isempty(pv), return; end
[sorted,order] = sort(pv(:)); m = numel(sorted); q_sorted = sorted .* m ./ (1:m)';
q_sorted = min(1,flipud(cummin(flipud(q_sorted)))); qv = zeros(size(pv)); qv(order) = q_sorted; q(valid) = qv;
end
