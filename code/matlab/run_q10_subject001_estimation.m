function summary = run_q10_subject001_estimation()
% Estimate Q10 only after run_q10_subject001_specification passes QA.

this_file = mfilename('fullpath');
model_root = fileparts(this_file);
first_level = fileparts(model_root);
glm_root = fileparts(first_level);
output_root = fullfile(model_root, 'q10_hiloadhe_vs_hiloadle', 'subject001');
nirs_file = fullfile(first_level, 'nirs_spm_inputs_v2_fulltimeline_trial1', ...
    'subject001', 'subject001_nirs_data.mat');
qa_file = fullfile(output_root, 'subject001_q10_design_qa_summary.csv');

assert(isfile(nirs_file), 'Missing NIRS input: %s', nirs_file);
assert(isfile(qa_file), 'Run Q10 specification and QA before estimation.');
qa = readtable(qa_file, 'TextType', 'string');
assert(all(qa.qa_pass), 'Q10 design QA did not pass. Estimation stopped.');

addpath(first_level);
addpath(fullfile(glm_root, 'toolboxes', 'spm'));
addpath(genpath(fullfile(glm_root, 'toolboxes', 'NIRS_SPM_v4_r1')));

input_data = load(nirs_file, 'nirs_data');
n_channels = input_data.nirs_data.nch;
clear input_data;
summary = table();
for hb_label = ["HbO", "HbR"]
    hb_root = fullfile(output_root, char(hb_label));
    spm_file = fullfile(hb_root, sprintf('SPM_indiv_%s.mat', hb_label));
    assert(isfile(spm_file), 'Missing Q10 specification: %s', spm_file);
    estimated = nirs_spm_estimation_batch_fixed(spm_file, nirs_file);
    assert(~isempty(estimated), '%s estimation returned empty output.', hb_label);
    check = load(spm_file, 'SPM_nirs');
    SPM_nirs = check.SPM_nirs;
    assert(strcmpi(string(SPM_nirs.nirs.step), 'estimation'), ...
        '%s did not reach estimation state.', hb_label);
    assert(all(size(SPM_nirs.nirs.beta) == [3, n_channels]), ...
        '%s beta dimensions are not 3 x %d.', hb_label, n_channels);
    assert(all(isfinite(SPM_nirs.nirs.beta), 'all'), '%s beta contains non-finite values.', hb_label);
    summary = [summary; table(hb_label, size(SPM_nirs.nirs.beta,1), ...
        size(SPM_nirs.nirs.beta,2), SPM_nirs.xX.erdf, ...
        'VariableNames', {'hemoglobin','n_betas','n_channels','effective_df'})]; %#ok<AGROW>
end
writetable(summary, fullfile(output_root, 'subject001_q10_estimation_summary.csv'));
fprintf('Q10 subject001 beta estimation complete for HbO and HbR.\n');
end
