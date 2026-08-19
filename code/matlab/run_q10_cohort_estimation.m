function [summary, errors] = run_q10_cohort_estimation(subject_numbers)
% Estimate previously QA-approved Q10 specifications. No contrast is computed here.

if nargin < 1, subject_numbers = 2:30; end
subject_numbers = subject_numbers(:)';
assert(all(ismember(subject_numbers, 2:30)), 'This cohort runner supports subject002:subject030.');

this_file = mfilename('fullpath');
model_root = fileparts(this_file);
first_level = fileparts(model_root);
glm_root = fileparts(first_level);
v2_root = fullfile(first_level, 'nirs_spm_inputs_v2_fulltimeline_trial1');
raw_root = fullfile(first_level, 'nirs_spm_inputs');
spec_root = fullfile(model_root, 'q10_hiloadhe_vs_hiloadle', 'cohort_specification');
assert(isfolder(spec_root), 'Missing Q10 cohort specification directory.');

addpath(first_level);
addpath(fullfile(glm_root, 'toolboxes', 'spm'));
addpath(genpath(fullfile(glm_root, 'toolboxes', 'NIRS_SPM_v4_r1')));

summary = table(); errors = table();
for subject_num = subject_numbers
    subject_label = sprintf('subject%03d', subject_num);
    if ismember(subject_num, [2, 10])
        nirs_file = fullfile(v2_root, subject_label, sprintf('%s_nirs_data.mat', subject_label));
    else
        nirs_file = fullfile(raw_root, subject_label, sprintf('%s_nirs_data.mat', subject_label));
    end
    assert(isfile(nirs_file), 'Missing NIRS input: %s', nirs_file);
    data = load(nirs_file, 'nirs_data'); n_channels = data.nirs_data.nch; clear data;
    fprintf('[%d] %s\n', subject_num, subject_label);
    for hb_label = ["HbO", "HbR"]
        spm_file = fullfile(spec_root, subject_label, char(hb_label), sprintf('SPM_indiv_%s.mat', hb_label));
        started = datetime('now'); timer = tic;
        status = "success"; message = ""; skipped_existing = false;
        beta_rows = NaN; beta_channels = NaN; erdf = NaN;
        try
            assert(isfile(spm_file), 'Missing Q10 specification: %s', spm_file);
            existing = load(spm_file, 'SPM_nirs');
            if isfield(existing, 'SPM_nirs') && isfield(existing.SPM_nirs, 'nirs') && ...
                    isfield(existing.SPM_nirs.nirs, 'step') && strcmpi(string(existing.SPM_nirs.nirs.step), 'estimation')
                SPM_nirs = existing.SPM_nirs;
                skipped_existing = true;
            else
                clear existing;
                estimated = nirs_spm_estimation_batch_fixed(spm_file, nirs_file);
                assert(~isempty(estimated), '%s estimation returned empty output.', hb_label);
                clear estimated;
                loaded = load(spm_file, 'SPM_nirs'); SPM_nirs = loaded.SPM_nirs;
            end
            assert(strcmpi(string(SPM_nirs.nirs.step), 'estimation'), '%s did not reach estimation state.', hb_label);
            beta_rows = size(SPM_nirs.nirs.beta,1);
            beta_channels = size(SPM_nirs.nirs.beta,2);
            erdf = double(SPM_nirs.xX.erdf);
            assert(isequal([beta_rows, beta_channels], [3, n_channels]), ...
                '%s beta dimensions are %d x %d; expected 3 x %d.', hb_label, beta_rows, beta_channels, n_channels);
            assert(all(isfinite(SPM_nirs.nirs.beta), 'all'), '%s beta contains non-finite values.', hb_label);
            fprintf('  %s beta %dx%d, df %.0f%s\n', hb_label, beta_rows, beta_channels, erdf, ternary(skipped_existing, ' (already estimated)', ''));
        catch ME
            status = "failed"; message = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
            errors = [errors; table(string(subject_label),subject_num,hb_label,string(ME.identifier),message, ...
                'VariableNames',{'subject_id','subject_num','hemoglobin','error_identifier','error_message'})]; %#ok<AGROW>
            fprintf(2, '  %s ERROR: %s\n', hb_label, ME.message);
        end
        summary = [summary; table(string(subject_label),subject_num,hb_label,status,skipped_existing, ...
            beta_rows,beta_channels,erdf,toc(timer),string(spm_file),string(started),string(datetime('now')), ...
            'VariableNames',{'subject_id','subject_num','hemoglobin','status','skipped_existing', ...
            'n_betas','n_channels','effective_df','elapsed_sec','spm_file','started_at','finished_at'})]; %#ok<AGROW>
    end
end

tag = sprintf('%03d_to_%03d', min(subject_numbers), max(subject_numbers));
writetable(summary, fullfile(spec_root, sprintf('q10_estimation_%s.csv', tag)));
writetable(errors, fullfile(spec_root, sprintf('q10_estimation_errors_%s.csv', tag)));
assert(isempty(errors), 'Q10 estimation had %d failed subject/chromophore models.', height(errors));
fprintf('Q10 estimation batch %s complete: %d models.\n', tag, height(summary));
end

function out = ternary(condition, if_true, if_false)
if condition, out = if_true; else, out = if_false; end
end
