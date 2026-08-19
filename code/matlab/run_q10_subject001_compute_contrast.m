function summary = run_q10_subject001_compute_contrast()
% Compute Q10: HILoadHE minus HILoadLE after first-level beta estimation.

this_file = mfilename('fullpath');
model_root = fileparts(this_file);
first_level = fileparts(model_root);
output_root = fullfile(model_root, 'q10_hiloadhe_vs_hiloadle', 'subject001');
metadata_file = fullfile(first_level, 'nirs_spm_inputs', ...
    'subject001', 'subject001_channel_metadata.csv');
assert(isfile(metadata_file), 'Missing channel metadata: %s', metadata_file);
metadata = readtable(metadata_file, 'TextType', 'string');

summary = table();
for hb_label = ["HbO", "HbR"]
    spm_file = fullfile(output_root, char(hb_label), sprintf('SPM_indiv_%s.mat', hb_label));
    assert(isfile(spm_file), 'Missing estimated SPM file: %s', spm_file);
    loaded = load(spm_file, 'SPM_nirs');
    SPM_nirs = loaded.SPM_nirs;
    assert(strcmpi(string(SPM_nirs.nirs.step), 'estimation'), ...
        '%s is not estimated yet.', hb_label);

    beta = double(SPM_nirs.nirs.beta);
    names = string(SPM_nirs.xX.name(:));
    idx_le = find(contains(names, 'HILoadLE'), 1);
    idx_he = find(contains(names, 'HILoadHE'), 1);
    assert(~isempty(idx_le) && ~isempty(idx_he), 'Q10 task regressors are missing.');
    c = zeros(size(beta,1), 1);
    c(idx_le) = -1;
    c(idx_he) = 1;

    resss = residual_sum_squares(SPM_nirs.nirs.ResSS, size(beta,2));
    erdf = double(SPM_nirs.xX.erdf);
    bcov = double(SPM_nirs.xX.Bcov);
    assert(isequal(size(bcov), [size(beta,1), size(beta,1)]), 'Bcov dimensions mismatch beta.');
    contrast_beta = c' * beta;
    contrast_variance = max((c' * bcov * c) .* (resss ./ erdf), 0);
    contrast_se = sqrt(contrast_variance);
    contrast_t = contrast_beta ./ contrast_se;
    contrast_p = two_sided_t_p(contrast_t, erdf);

    n_channels = size(beta,2);
    assert(height(metadata) == n_channels, 'Metadata channel count mismatch.');
    channel_index = (1:n_channels)';
    result = table(repmat("subject001", n_channels, 1), repmat(hb_label, n_channels, 1), ...
        channel_index, string(metadata.channel_id), string(metadata.channel_type), ...
        string(metadata.roi), string(metadata.laterality), string(metadata.qc_status), ...
        contrast_beta(:), contrast_se(:), contrast_t(:), contrast_p(:), ...
        contrast_variance(:), repmat(erdf, n_channels, 1), repmat(join(string(c'), ' '), n_channels, 1), ...
        'VariableNames', {'subject_id','hemoglobin','channel_index','channel_id', ...
        'channel_type','roi','laterality','qc_status','contrast_beta', ...
        'contrast_se','contrast_t','contrast_p_two_sided','contrast_variance', ...
        'effective_df','contrast_vector'});
    writetable(result, fullfile(output_root, sprintf('subject001_%s_q10_HE_minus_LE_channel_contrasts.csv', hb_label)));
    summary = [summary; table(hb_label, n_channels, sum(isfinite(contrast_beta)), ...
        sum(isfinite(contrast_t)), min(contrast_p), median(contrast_p), ...
        'VariableNames', {'hemoglobin','n_channels','n_finite_contrast_beta', ...
        'n_finite_t','min_p_two_sided','median_p_two_sided'})]; %#ok<AGROW>
end

writetable(summary, fullfile(output_root, 'subject001_q10_contrast_summary.csv'));
definition = table("Q10_HILoadHE_minus_HILoadLE", "HILoadHE - HILoadLE", ...
    "[-1, 1, 0] in the implemented order [HILoadLE, HILoadHE, constant]", ...
    'VariableNames', {'contrast_name','description','implemented_vector'});
writetable(definition, fullfile(output_root, 'subject001_q10_contrast_definition.csv'));
fprintf('Q10 subject001 contrast complete for HbO and HbR.\n');
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
% Equivalent to 2*tcdf(-abs(t),df), using base MATLAB betainc.
p = betainc(df ./ (df + t.^2), df/2, 1/2);
p(~isfinite(t)) = NaN;
end
