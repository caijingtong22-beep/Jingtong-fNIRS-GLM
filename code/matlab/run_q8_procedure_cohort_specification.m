function [summary, errors] = run_q8_procedure_cohort_specification(subject_numbers)
% Q8 Procedure specification/HRF convolution only: HIHE minus LILE.
% The Procedure events are used at their recorded full durations.  DOM and
% AOE are deliberately not removed; they are analysed in separate GLMs.
% This function does not estimate beta parameters or compute contrasts.

if nargin < 1, subject_numbers = 1; end
subject_numbers = subject_numbers(:)';
assert(all(ismember(subject_numbers, 1:30)), 'Subject numbers must be in 1:30.');

this_file = mfilename('fullpath');
model_root = fileparts(this_file);
repo_root = fileparts(fileparts(model_root));
glm_root = 'C:\\aliceproject\\GLM';
first_level = fullfile(glm_root, 'first_level');
v2_root = fullfile(first_level, 'nirs_spm_inputs_v2_fulltimeline_trial1');
raw_root = fullfile(first_level, 'nirs_spm_inputs');
event_table_file = fullfile(glm_root, 'amended_snirf', 'event_tables_full_annotations', ...
    'event_table_all_subjects_full_annotations.csv');
output_root = fullfile(repo_root, 'results', 'q8_hihe_vs_lile_procedure', 'cohort_specification');
assert(isfile(event_table_file), 'Missing full event table: %s', event_table_file);
if ~isfolder(output_root), mkdir(output_root); end

addpath(first_level);
addpath(fullfile(glm_root, 'toolboxes', 'spm'));
addpath(genpath(fullfile(glm_root, 'toolboxes', 'NIRS_SPM_v4_r1')));
full_events = readtable(event_table_file, 'TextType', 'string');
full_events.condition = string(full_events.condition);
full_events.subject_num = double(full_events.subject_id);

summary = table(); errors = table(); manifest = table();
for subject_num = subject_numbers
    subject_label = sprintf('subject%03d', subject_num);
    [nirs_file, events, policy] = q8_inputs(subject_num, subject_label, v2_root, raw_root, full_events);
    subject_root = fullfile(output_root, subject_label);
    if ~isfolder(subject_root), mkdir(subject_root); end
    try
        [names, onsets, durations, selected] = q8_procedure_conditions(events, subject_label);
        writetable(selected, fullfile(subject_root, sprintf('%s_q8_procedure_events_used.csv', subject_label)));
        save(fullfile(subject_root, sprintf('%s_q8_procedure_multiple_conditions.mat', subject_label)), ...
            'names', 'onsets', 'durations', '-v7');
        manifest = [manifest; table(string(subject_label), subject_num, string(policy), string(nirs_file), ...
            'VariableNames', {'subject_id','subject_num','selection_policy','nirs_data_file'})]; %#ok<AGROW>
        fprintf('[%d] %s (%s)\n', subject_num, subject_label, policy);
        for hb = ["hbo", "hbr"]
            if hb == "hbo", hb_label = "HbO"; else, hb_label = "HbR"; end
            hb_root = fullfile(subject_root, char(hb_label));
            if ~isfolder(hb_root), mkdir(hb_root); end
            % specification_batch HRF-convolves the two named Procedure events.
            SPM_nirs = specification_batch(nirs_file, char(hb), 'DCT, 128', ...
                'none', 0, hb_root, 0, 0, 1, names, onsets, durations);
            qa = assess_q8_design(SPM_nirs, string(subject_label), string(hb_label), hb_root);
            summary = [summary; qa]; %#ok<AGROW>
            fprintf('  %s rank %d/%d, r %.3f, VIF %.3f, pass=%d\n', ...
                hb_label, qa.rank, qa.n_columns, qa.max_abs_correlation, qa.max_vif, qa.qa_pass);
        end
    catch ME
        errors = [errors; table(string(subject_label), subject_num, string(ME.identifier), string(ME.message), ...
            'VariableNames', {'subject_id','subject_num','error_identifier','error_message'})]; %#ok<AGROW>
        fprintf(2, '  ERROR %s: %s\n', subject_label, ME.message);
    end
end

tag = sprintf('%03d_to_%03d', min(subject_numbers), max(subject_numbers));
writetable(manifest, fullfile(output_root, sprintf('q8_procedure_manifest_%s.csv', tag)));
writetable(summary, fullfile(output_root, sprintf('q8_procedure_design_qa_%s.csv', tag)));
writetable(errors, fullfile(output_root, sprintf('q8_procedure_specification_errors_%s.csv', tag)));
assert(isempty(errors), 'Q8 Procedure specification had %d failed subjects.', height(errors));
assert(all(summary.qa_pass), 'At least one Q8 Procedure design failed QA.');
end

function [nirs_file, events, policy] = q8_inputs(subject_num, subject_label, v2_root, raw_root, full_events)
if ismember(subject_num, [1, 2, 10])
    source_dir = fullfile(v2_root, subject_label);
    nirs_file = fullfile(source_dir, sprintf('%s_nirs_data.mat', subject_label));
    event_file = fullfile(source_dir, sprintf('%s_events_used_v2.csv', subject_label));
    assert(isfile(event_file), 'Missing V2 events: %s', event_file);
    events = readtable(event_file, 'TextType', 'string'); policy = "first_complete_cycle";
else
    source_dir = fullfile(raw_root, subject_label);
    nirs_file = fullfile(source_dir, sprintf('%s_nirs_data.mat', subject_label));
    events = full_events(full_events.subject_num == subject_num, :); policy = "full_recording";
end
assert(isfile(nirs_file), 'Missing NIRS input: %s', nirs_file);
events.condition = string(events.condition);
end

function [names, onsets, durations, selected] = q8_procedure_conditions(events, subject_label)
% Full Procedure intervals; no DOM/AOE subtraction or trimming is applied.
target_names = ["LILEProcedure", "HIHEProcedure"];
selected = sortrows(events(ismember(events.condition, target_names), :), 'onset_sec');
assert(all(ismember(target_names, unique(selected.condition))), ...
    '%s is missing LILEProcedure or HIHEProcedure.', subject_label);
names = cellstr(target_names); onsets = cell(1, 2); durations = cell(1, 2);
for i = 1:2
    rows = selected.condition == target_names(i);
    assert(sum(rows) == 1, '%s must have exactly one %s occurrence after selection.', ...
        subject_label, target_names(i));
    onsets{i} = double(selected.onset_sec(rows))';
    durations{i} = double(selected.duration_sec(rows))';
end
end

function qa = assess_q8_design(SPM_nirs, subject_label, hb_label, output_root)
X = full(SPM_nirs.xX.X); regressors = string(SPM_nirs.xX.name(:));
is_constant = contains(lower(regressors), 'constant');
assert(size(X,2) == 3 && sum(is_constant) == 1, ...
    'Expected LILEProcedure, HIHEProcedure, and one constant column.');
task_X = X(:, ~is_constant);
task_Z = (task_X - mean(task_X,1)) ./ std(task_X,0,1);
assert(all(isfinite(task_Z), 'all'), 'Task regressors have zero/non-finite variance.');
R = corrcoef(task_Z); max_r = abs(R(1,2));
vif = zeros(2,1);
for j = 1:2
    y = task_X(:,j); z = [ones(size(X,1),1), task_X(:,setdiff(1:2,j))];
    r2 = 1 - sum((y-z*(z\y)).^2) / sum((y-mean(y)).^2);
    vif(j) = 1 / max(1-r2, eps);
end
rank_x = rank(X); qa_pass = rank_x == size(X,2) && max_r < 0.95 && max(vif) <= 10;
writematrix(X, fullfile(output_root, sprintf('%s_%s_q8_procedure_convolved_design.csv', subject_label, hb_label)));
writetable(table((1:numel(regressors))', regressors, is_constant, ...
    'VariableNames', {'column_index','regressor','is_constant'}), ...
    fullfile(output_root, sprintf('%s_%s_q8_procedure_regressor_table.csv', subject_label, hb_label)));
writetable(table(regressors(~is_constant), vif, 'VariableNames', {'regressor','vif'}), ...
    fullfile(output_root, sprintf('%s_%s_q8_procedure_vif.csv', subject_label, hb_label)));
qa = table(subject_label, hb_label, size(X,1), size(X,2), rank_x, cond(task_Z), max_r, max(vif), qa_pass, ...
    'VariableNames', {'subject_id','hemoglobin','n_scans','n_columns','rank', ...
    'normalized_condition_number','max_abs_correlation','max_vif','qa_pass'});
end
