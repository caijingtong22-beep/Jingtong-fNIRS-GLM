function [summary, errors] = run_q10_cohort_compute_contrast(subject_numbers)
% Compute Q10 HILoadHE minus HILoadLE contrasts after beta estimation.
% This produces individual first-level results only; it does not perform
% a group-level analysis.

if nargin < 1
    subject_numbers = 2:30;
end
subject_numbers = unique(subject_numbers(:)');
assert(all(subject_numbers >= 2 & subject_numbers <= 30), ...
    'This cohort script accepts subject numbers 2 through 30.');

this_file = mfilename('fullpath');
model_root = fileparts(this_file);
first_level = fileparts(model_root);
output_root = fullfile(model_root, 'q10_hiloadhe_vs_hiloadle', 'cohort_specification');
raw_input_root = fullfile(first_level, 'nirs_spm_inputs');

summary = table();
errors = table();
for subject_num = subject_numbers
    subject_label = sprintf('subject%03d', subject_num);
    fprintf('[%d] %s\n', subject_num, subject_label);
    % V2 holds only the trial-1 signal input.  Channel metadata is shared
    % with the raw input and is stored there for every subject.
    metadata_file = fullfile(raw_input_root, subject_label, ...
        sprintf('%s_channel_metadata.csv', subject_label));
    try
        assert(isfile(metadata_file), 'Missing channel metadata: %s', metadata_file);
        metadata = readtable(metadata_file, 'TextType', 'string');
        for hb_label = ["HbO", "HbR"]
            spm_file = fullfile(output_root, subject_label, char(hb_label), ...
                sprintf('SPM_indiv_%s.mat', hb_label));
            assert(isfile(spm_file), 'Missing estimated SPM file: %s', spm_file);
            loaded = load(spm_file, 'SPM_nirs');
            SPM_nirs = loaded.SPM_nirs;
            assert(strcmpi(string(SPM_nirs.nirs.step), 'estimation'), ...
                '%s %s is not estimated.', subject_label, hb_label);

            beta = double(SPM_nirs.nirs.beta);
            names = string(SPM_nirs.xX.name(:));
            idx_le = find(contains(names, 'HILoadLE'), 1);
            idx_he = find(contains(names, 'HILoadHE'), 1);
            assert(~isempty(idx_le) && ~isempty(idx_he), ...
                'Q10 task regressors are missing in %s %s.', subject_label, hb_label);
            c = zeros(size(beta,1), 1);
            c(idx_le) = -1;
            c(idx_he) = 1;

            n_channels = size(beta,2);
            assert(height(metadata) == n_channels, 'Metadata channel count mismatch.');
            resss = residual_sum_squares(SPM_nirs.nirs.ResSS, n_channels);
            erdf = double(SPM_nirs.xX.erdf);
            bcov = double(SPM_nirs.xX.Bcov);
            assert(isequal(size(bcov), [size(beta,1), size(beta,1)]), ...
                'Bcov dimensions mismatch beta.');
            contrast_beta = c' * beta;
            contrast_variance = max((c' * bcov * c) .* (resss ./ erdf), 0);
            contrast_se = sqrt(contrast_variance);
            contrast_t = contrast_beta ./ contrast_se;
            contrast_p = two_sided_t_p(contrast_t, erdf);
            contrast_q = bh_fdr(contrast_p);

            channel_index = (1:n_channels)';
            result = table(repmat(string(subject_label), n_channels, 1), ...
                repmat(hb_label, n_channels, 1), channel_index, ...
                string(metadata.channel_id), string(metadata.channel_type), ...
                string(metadata.roi), string(metadata.laterality), ...
                string(metadata.qc_status), contrast_beta(:), contrast_se(:), ...
                contrast_t(:), contrast_p(:), contrast_q(:), contrast_variance(:), ...
                repmat(erdf, n_channels, 1), repmat(join(string(c'), ' '), n_channels, 1), ...
                'VariableNames', {'subject_id','hemoglobin','channel_index','channel_id', ...
                'channel_type','roi','laterality','qc_status','contrast_beta', ...
                'contrast_se','contrast_t','contrast_p_two_sided','contrast_q_bh61', ...
                'contrast_variance','effective_df','contrast_vector'});
            subject_out = fullfile(output_root, subject_label);
            writetable(result, fullfile(subject_out, ...
                sprintf('%s_%s_q10_HE_minus_LE_channel_contrasts.csv', subject_label, hb_label)));
            summary = [summary; table(string(subject_label), hb_label, n_channels, ...
                sum(isfinite(contrast_beta)), sum(contrast_p < 0.05), sum(contrast_q < 0.05), ...
                min(contrast_p), median(contrast_p), ...
                'VariableNames', {'subject_id','hemoglobin','n_channels','n_finite_beta', ...
                'n_p_lt_005','n_q_bh61_lt_005','min_p_two_sided','median_p_two_sided'})]; %#ok<AGROW>
            fprintf('  %s contrast: %d p<.05, %d BH-FDR<.05\n', ...
                hb_label, sum(contrast_p < 0.05), sum(contrast_q < 0.05));
        end
    catch ME
        errors = [errors; table(string(subject_label), string(ME.identifier), string(ME.message), ...
            'VariableNames', {'subject_id','error_id','error_message'})]; %#ok<AGROW>
        fprintf(2, 'ERROR %s: %s\n', subject_label, ME.message);
    end
end

batch_label = sprintf('%03d_to_%03d', min(subject_numbers), max(subject_numbers));
writetable(summary, fullfile(output_root, sprintf('q10_contrast_%s.csv', batch_label)));
writetable(errors, fullfile(output_root, sprintf('q10_contrast_errors_%s.csv', batch_label)));
definition = table("Q10_HILoadHE_minus_HILoadLE", "HILoadHE - HILoadLE", ...
    "[-1, 1, 0] in implemented order [HILoadLE, HILoadHE, constant]", ...
    'VariableNames', {'contrast_name','description','implemented_vector'});
writetable(definition, fullfile(output_root, 'q10_contrast_definition.csv'));
fprintf('Q10 contrast batch %s complete: %d model results, %d subject errors.\n', ...
    batch_label, height(summary), height(errors));
end

function resss = residual_sum_squares(raw, n_channels)
raw = double(raw);
if isvector(raw)
    resss = raw(:)';
elseif isequal(size(raw), [n_channels, n_channels])
    resss = diag(raw)';
else
    error('Unexpected ResSS dimensions %dx%d.', size(raw,1), size(raw,2));
end
assert(numel(resss) == n_channels, 'ResSS channel count mismatch.');
end

function p = two_sided_t_p(t, df)
% Equivalent to 2*tcdf(-abs(t), df), using base MATLAB betainc.
p = betainc(df ./ (df + t.^2), df/2, 1/2);
p(~isfinite(t)) = NaN;
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
