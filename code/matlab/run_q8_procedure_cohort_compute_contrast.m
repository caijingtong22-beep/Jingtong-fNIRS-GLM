function [summary, errors] = run_q8_procedure_cohort_compute_contrast(subject_numbers)
% Q8 Procedure first-level contrast: HIHEProcedure minus LILEProcedure.
% Creates per-channel results only; no group-level inference is performed.

if nargin < 1, subject_numbers = 1; end
subject_numbers = unique(subject_numbers(:)');
assert(all(ismember(subject_numbers, 1:30)), 'Subject numbers must be in 1:30.');

this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(fileparts(this_file)));
glm_root = 'C:\\aliceproject\\GLM';
output_root = fullfile(repo_root, 'results', 'q8_hihe_vs_lile_procedure', 'cohort_specification');
raw_input_root = fullfile(glm_root, 'first_level', 'nirs_spm_inputs');
summary = table(); errors = table();

for subject_num = subject_numbers
    subject_label = sprintf('subject%03d', subject_num);
    try
        metadata_file = fullfile(raw_input_root, subject_label, sprintf('%s_channel_metadata.csv', subject_label));
        assert(isfile(metadata_file), 'Missing channel metadata: %s', metadata_file);
        metadata = readtable(metadata_file, 'TextType', 'string');
        for hb_label = ["HbO", "HbR"]
            subject_out = fullfile(output_root, subject_label);
            spm_file = fullfile(subject_out, char(hb_label), sprintf('SPM_indiv_%s.mat', hb_label));
            loaded = load(spm_file, 'SPM_nirs'); SPM_nirs = loaded.SPM_nirs;
            assert(strcmpi(string(SPM_nirs.nirs.step), 'estimation'), '%s %s is not estimated.', subject_label, hb_label);
            beta = double(SPM_nirs.nirs.beta); names = string(SPM_nirs.xX.name(:));
            idx_lile = find(contains(names, 'LILEProcedure'), 1);
            idx_hihe = find(contains(names, 'HIHEProcedure'), 1);
            assert(~isempty(idx_lile) && ~isempty(idx_hihe), 'Q8 Procedure task regressors are missing.');
            c = zeros(size(beta,1), 1); c(idx_lile) = -1; c(idx_hihe) = +1;
            n_channels = size(beta,2);
            assert(height(metadata) == n_channels, 'Metadata channel count mismatch.');
            resss = residual_sum_squares(SPM_nirs.nirs.ResSS, n_channels);
            erdf = double(SPM_nirs.xX.erdf); bcov = double(SPM_nirs.xX.Bcov);
            assert(isequal(size(bcov), [size(beta,1), size(beta,1)]), 'Bcov dimensions mismatch beta.');
            contrast_beta = c' * beta;
            contrast_variance = max((c' * bcov * c) .* (resss ./ erdf), 0);
            contrast_se = sqrt(contrast_variance); contrast_t = contrast_beta ./ contrast_se;
            contrast_p = two_sided_t_p(contrast_t, erdf); contrast_q = bh_fdr(contrast_p);
            result = table(repmat(string(subject_label),n_channels,1), repmat(hb_label,n_channels,1), ...
                (1:n_channels)', string(metadata.channel_id), string(metadata.channel_type), string(metadata.roi), ...
                string(metadata.laterality), string(metadata.qc_status), contrast_beta(:), contrast_se(:), ...
                contrast_t(:), contrast_p(:), contrast_q(:), contrast_variance(:), repmat(erdf,n_channels,1), ...
                repmat(join(string(c'), ' '),n_channels,1), ...
                'VariableNames', {'subject_id','hemoglobin','channel_index','channel_id','channel_type','roi', ...
                'laterality','qc_status','contrast_beta','contrast_se','contrast_t','contrast_p_two_sided', ...
                'contrast_q_bh61','contrast_variance','effective_df','contrast_vector'});
            writetable(result, fullfile(subject_out, sprintf('%s_%s_q8_HIHE_minus_LILE_procedure_channel_contrasts.csv', ...
                subject_label, hb_label)));
            summary = [summary; table(string(subject_label), hb_label, n_channels, sum(isfinite(contrast_beta)), ...
                sum(contrast_p < .05), sum(contrast_q < .05), min(contrast_p), median(contrast_p), ...
                'VariableNames', {'subject_id','hemoglobin','n_channels','n_finite_beta','n_p_lt_005', ...
                'n_q_bh61_lt_005','min_p_two_sided','median_p_two_sided'})]; %#ok<AGROW>
            fprintf('%s %s: %d p<.05, %d BH-FDR<.05\n', subject_label, hb_label, ...
                sum(contrast_p < .05), sum(contrast_q < .05));
        end
    catch ME
        errors = [errors; table(string(subject_label),string(ME.identifier),string(ME.message), ...
            'VariableNames', {'subject_id','error_id','error_message'})]; %#ok<AGROW>
    end
end
tag = sprintf('%03d_to_%03d', min(subject_numbers), max(subject_numbers));
writetable(summary, fullfile(output_root, sprintf('q8_procedure_contrast_%s.csv', tag)));
writetable(errors, fullfile(output_root, sprintf('q8_procedure_contrast_errors_%s.csv', tag)));
definition = table("Q8_HIHE_minus_LILE_Procedure", "HIHEProcedure - LILEProcedure", ...
    "[-1, +1, 0], located by regressor name; implemented order [LILEProcedure, HIHEProcedure, constant]", ...
    'VariableNames', {'contrast_name','description','implemented_vector'});
writetable(definition, fullfile(output_root, 'q8_procedure_contrast_definition.csv'));
assert(isempty(errors), 'Q8 Procedure contrast had %d failures.', height(errors));
end

function resss = residual_sum_squares(raw, n_channels)
raw = double(raw);
if isvector(raw), resss = raw(:)';
elseif isequal(size(raw), [n_channels,n_channels]), resss = diag(raw)';
else, error('Unexpected ResSS dimensions %dx%d.', size(raw,1), size(raw,2));
end
assert(numel(resss) == n_channels, 'ResSS channel count mismatch.');
end

function p = two_sided_t_p(t, df)
p = betainc(df ./ (df + t.^2), df/2, 1/2); p(~isfinite(t)) = NaN;
end

function q = bh_fdr(p)
q = nan(size(p)); valid = isfinite(p); pv = p(valid); if isempty(pv), return; end
[sorted, order] = sort(pv(:)); m = numel(sorted); q_sorted = sorted .* m ./ (1:m)';
q_sorted = min(1, flipud(cummin(flipud(q_sorted)))); qv = zeros(size(pv)); qv(order) = q_sorted; q(valid) = qv;
end
