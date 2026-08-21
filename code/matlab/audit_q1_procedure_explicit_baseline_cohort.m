function manifest = audit_q1_procedure_explicit_baseline_cohort()
% Q1 Procedure cohort eligibility audit for explicit ExperimentalBaseline.
% Target-event policy: use the first LILEProcedure occurrence by onset.
repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
glm = 'C:\aliceproject\GLM';
event_file = fullfile(glm, 'amended_snirf', 'event_tables_full_annotations', ...
    'event_table_all_subjects_full_annotations.csv');
input_root = fullfile(glm, 'first_level', 'nirs_spm_inputs');
output_root = fullfile(repo, 'results', 'q1_lile_vs_experimentalbaseline_procedure', 'cohort_specification');
if ~isfolder(output_root), mkdir(output_root); end

events = readtable(event_file, 'TextType', 'string');
events.subject_id = compose('%03d', double(events.subject_id));
events.condition = string(events.condition);
target_label = "LILEProcedure";
manifest = table();
selected = table();
for id = 1:30
    subject_id = string(sprintf('subject%03d', id));
    event_id = string(sprintf('%03d', id));
    input_file = fullfile(input_root, char(subject_id), char(subject_id + "_nirs_data.mat"));
    sub = events(events.subject_id == event_id, :);
    baseline = sortrows(sub(sub.condition == "ExperimentalBaseline", :), 'onset_sec');
    target = sortrows(sub(sub.condition == target_label, :), 'onset_sec');
    reason = "included";
    if ~isfile(input_file)
        reason = "continuous input missing";
    elseif height(baseline) ~= 1
        reason = sprintf('ExperimentalBaseline count is %d (requires exactly 1)', height(baseline));
    elseif height(target) < 1
        reason = sprintf('%s count is 0', target_label);
    end
    included = reason == "included";
    first_target_onset = NaN; first_target_duration = NaN;
    baseline_onset = NaN; baseline_duration = NaN;
    if height(baseline) == 1
        baseline_onset = double(baseline.onset_sec); baseline_duration = double(baseline.duration_sec);
    end
    if height(target) >= 1
        first_target_onset = double(target.onset_sec(1)); first_target_duration = double(target.duration_sec(1));
    end
    manifest = [manifest; table(subject_id, included, string(reason), isfile(input_file), height(baseline), height(target), ...
        baseline_onset, baseline_duration, first_target_onset, first_target_duration, ...
        repmat("first occurrence by onset_sec (user-confirmed)",1,1), ...
        'VariableNames', {'subject_id','included','reason','input_file_exists','experimentalbaseline_count','lileprocedure_count','baseline_onset_sec','baseline_duration_sec','selected_target_onset_sec','selected_target_duration_sec','target_selection_policy'})]; %#ok<AGROW>
    if included
        b = baseline; b.role = "ExperimentalBaseline_boxcar"; b.selection_policy = "single ExperimentalBaseline event";
        t = target(1,:); t.role = "target_phase_HRF"; t.selection_policy = "first occurrence by onset_sec (user-confirmed)";
        b.input_file = string(input_file); t.input_file = string(input_file);
        selected = [selected; b; t]; %#ok<AGROW>
    end
end
writetable(manifest, fullfile(output_root, 'q1_procedure_explicit_baseline_inclusion_manifest.csv'));
writetable(selected, fullfile(output_root, 'q1_procedure_explicit_baseline_selected_events.csv'));
fprintf('Q1 Procedure explicit-baseline eligible cohort: n=%d; excluded: n=%d.\n', sum(manifest.included), sum(~manifest.included));
end
