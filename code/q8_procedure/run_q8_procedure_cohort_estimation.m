function [summary, errors] = run_q8_procedure_cohort_estimation(subject_numbers)
% Estimate QA-approved Q8 Procedure specifications; no contrast is computed.

if nargin < 1, subject_numbers = 1; end
subject_numbers = subject_numbers(:)';
assert(all(ismember(subject_numbers, 1:30)), 'Subject numbers must be in 1:30.');

this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(fileparts(this_file)));
glm_root = 'C:\\aliceproject\\GLM';
first_level = fullfile(glm_root, 'first_level');
spec_root = fullfile(repo_root, 'results', 'q8_hihe_vs_lile_procedure', 'cohort_specification');
assert(isfolder(spec_root), 'Missing Q8 Procedure specification directory.');
addpath(first_level);
addpath(fullfile(glm_root, 'toolboxes', 'spm'));
addpath(genpath(fullfile(glm_root, 'toolboxes', 'NIRS_SPM_v4_r1')));

summary = table(); errors = table();
for subject_num = subject_numbers
    subject_label = sprintf('subject%03d', subject_num);
    data_file = q8_data_file(subject_num, subject_label, first_level);
    data = load(data_file, 'nirs_data'); n_channels = data.nirs_data.nch; clear data;
    for hb_label = ["HbO", "HbR"]
        spm_file = fullfile(spec_root, subject_label, char(hb_label), sprintf('SPM_indiv_%s.mat', hb_label));
        timer = tic; status = "success"; message = ""; skipped_existing = false;
        beta_rows = NaN; beta_channels = NaN; erdf = NaN;
        try
            assert(isfile(spm_file), 'Missing Q8 Procedure specification: %s', spm_file);
            loaded = load(spm_file, 'SPM_nirs'); SPM_nirs = loaded.SPM_nirs;
            if isfield(SPM_nirs.nirs, 'step') && strcmpi(string(SPM_nirs.nirs.step), 'estimation')
                skipped_existing = true;
            else
                nirs_spm_estimation_batch_fixed(spm_file, data_file);
                loaded = load(spm_file, 'SPM_nirs'); SPM_nirs = loaded.SPM_nirs;
            end
            assert(strcmpi(string(SPM_nirs.nirs.step), 'estimation'), '%s did not reach estimation.', hb_label);
            beta_rows = size(SPM_nirs.nirs.beta, 1); beta_channels = size(SPM_nirs.nirs.beta, 2);
            erdf = double(SPM_nirs.xX.erdf);
            assert(isequal([beta_rows beta_channels], [3 n_channels]), ...
                '%s beta dimensions are %d x %d; expected 3 x %d.', hb_label, beta_rows, beta_channels, n_channels);
            assert(all(isfinite(SPM_nirs.nirs.beta), 'all'), '%s beta contains non-finite values.', hb_label);
            fprintf('%s %s: beta %dx%d, df %.0f%s\n', subject_label, hb_label, beta_rows, beta_channels, ...
                erdf, ternary(skipped_existing, ' (already estimated)', ''));
        catch ME
            status = "failed"; message = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
            errors = [errors; table(string(subject_label), subject_num, hb_label, string(ME.identifier), message, ...
                'VariableNames', {'subject_id','subject_num','hemoglobin','error_identifier','error_message'})]; %#ok<AGROW>
        end
        summary = [summary; table(string(subject_label), subject_num, hb_label, status, skipped_existing, ...
            beta_rows, beta_channels, erdf, toc(timer), ...
            'VariableNames', {'subject_id','subject_num','hemoglobin','status','skipped_existing', ...
            'n_betas','n_channels','effective_df','elapsed_sec'})]; %#ok<AGROW>
    end
end
tag = sprintf('%03d_to_%03d', min(subject_numbers), max(subject_numbers));
writetable(summary, fullfile(spec_root, sprintf('q8_procedure_estimation_%s.csv', tag)));
writetable(errors, fullfile(spec_root, sprintf('q8_procedure_estimation_errors_%s.csv', tag)));
assert(isempty(errors), 'Q8 Procedure estimation had %d failures.', height(errors));
end

function data_file = q8_data_file(subject_num, subject_label, first_level)
if ismember(subject_num, [1, 2, 10])
    data_file = fullfile(first_level, 'nirs_spm_inputs_v2_fulltimeline_trial1', subject_label, ...
        sprintf('%s_nirs_data.mat', subject_label));
else
    data_file = fullfile(first_level, 'nirs_spm_inputs', subject_label, sprintf('%s_nirs_data.mat', subject_label));
end
assert(isfile(data_file), 'Missing NIRS input: %s', data_file);
end

function value = ternary(condition, if_true, if_false)
if condition, value = if_true; else, value = if_false; end
end
