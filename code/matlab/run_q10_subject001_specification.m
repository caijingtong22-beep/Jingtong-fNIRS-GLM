function summary = run_q10_subject001_specification()
% Q10 pilot: HILoadHE minus HILoadLE, using subject001's first complete cycle.
% This is intentionally independent of any explicit/global baseline choice.

this_file = mfilename('fullpath');
model_root = fileparts(this_file);
first_level = fileparts(model_root);
glm_root = fileparts(first_level);
source_root = fullfile(first_level, 'nirs_spm_inputs_v2_fulltimeline_trial1', 'subject001');
output_root = fullfile(model_root, 'q10_hiloadhe_vs_hiloadle', 'subject001');
nirs_file = fullfile(source_root, 'subject001_nirs_data.mat');
source_events = fullfile(source_root, 'subject001_events_used_v2.csv');

assert(isfile(nirs_file), 'Missing NIRS input: %s', nirs_file);
assert(isfile(source_events), 'Missing source events: %s', source_events);
if ~isfolder(output_root), mkdir(output_root); end

addpath(first_level);
addpath(fullfile(glm_root, 'toolboxes', 'spm'));
addpath(genpath(fullfile(glm_root, 'toolboxes', 'NIRS_SPM_v4_r1')));

events = readtable(source_events, 'TextType', 'string');
events.condition = string(events.condition);
target_names = ["HILoadLE", "HILoadHE"];
selected = events(ismember(events.condition, target_names), :);
selected = sortrows(selected, 'onset_sec');
assert(all(ismember(target_names, unique(selected.condition))), ...
    'Q10 requires both HILoadLE and HILoadHE.');

% specification_batch iterates across columns, so keep this as a 1x2 cell array.
names = cellstr(target_names);
onsets = cell(1, numel(names));
durations = cell(1, numel(names));
for i = 1:numel(names)
    rows = selected.condition == string(names{i});
    assert(sum(rows) == 1, ...
        'subject001 first-cycle input must contain exactly one %s event.', names{i});
    onsets{i} = double(selected.onset_sec(rows))';
    durations{i} = double(selected.duration_sec(rows))';
end

writetable(selected, fullfile(output_root, 'subject001_q10_events_used.csv'));
save(fullfile(output_root, 'subject001_q10_multiple_conditions.mat'), ...
    'names', 'onsets', 'durations', '-v7');

summary = table();
for hb = ["hbo", "hbr"]
    if hb == "hbo"
        hb_label = "HbO";
    else
        hb_label = "HbR";
    end
    hb_root = fullfile(output_root, char(hb_label));
    if ~isfolder(hb_root), mkdir(hb_root); end

    % NIRS-SPM applies its canonical HRF to the two task regressors here.
    % The third design column is NIRS-SPM's constant/intercept.
    SPM_nirs = specification_batch(nirs_file, char(hb), 'DCT, 128', ...
        'none', 0, hb_root, 0, 0, 1, names, onsets, durations);
    assert(isfile(fullfile(hb_root, sprintf('SPM_indiv_%s.mat', hb_label))), ...
        'Missing NIRS-SPM specification output for %s.', hb_label);
    qa = assess_q10_design(SPM_nirs, string(hb_label), hb_root);
    summary = [summary; qa]; %#ok<AGROW>
end

writetable(summary, fullfile(output_root, 'subject001_q10_design_qa_summary.csv'));
assert(all(summary.qa_pass), 'Q10 design QA failed; do not run estimation.');
fprintf('Q10 subject001 specification complete. Both HbO and HbR passed QA.\n');
end

function qa = assess_q10_design(SPM_nirs, hb_label, output_root)
X = full(SPM_nirs.xX.X);
regressors = string(SPM_nirs.xX.name(:));
is_constant = contains(lower(regressors), 'constant');
assert(size(X, 2) == 3 && sum(is_constant) == 1, ...
    'Expected exactly HILoadLE, HILoadHE, and one constant column.');
task_X = X(:, ~is_constant);
task_Z = (task_X - mean(task_X, 1)) ./ std(task_X, 0, 1);
assert(all(isfinite(task_Z), 'all'), 'Task regressors have zero or non-finite variance.');
R = corrcoef(task_Z);
max_r = abs(R(1,2));
vif = zeros(2,1);
for j = 1:2
    y = task_X(:,j);
    z = [ones(size(task_X,1),1), task_X(:, setdiff(1:2,j))];
    coef = z \ y;
    r2 = 1 - sum((y-z*coef).^2) / sum((y-mean(y)).^2);
    vif(j) = 1 / max(1-r2, eps);
end
rank_x = rank(X);
qa_pass = rank_x == size(X,2) && max_r < 0.95 && max(vif) <= 10;

writetable(table((1:numel(regressors))', regressors, is_constant, ...
    'VariableNames', {'column_index','regressor','is_constant'}), ...
    fullfile(output_root, sprintf('subject001_%s_q10_regressor_table.csv', hb_label)));
writetable(table(regressors(~is_constant), vif, ...
    'VariableNames', {'regressor','vif'}), ...
    fullfile(output_root, sprintf('subject001_%s_q10_vif.csv', hb_label)));

fig = figure('Visible', 'off', 'Color', 'w');
imagesc(X'); axis tight; colormap(gray); colorbar;
yticks(1:numel(regressors)); yticklabels(regressors);
xlabel('Scan'); ylabel('Design-matrix column');
title(sprintf('subject001 Q10 design matrix (%s): HILoadHE - HILoadLE', hb_label));
exportgraphics(fig, fullfile(output_root, sprintf('subject001_%s_q10_design_matrix.png', hb_label)), 'Resolution', 180);
close(fig);

qa = table(hb_label, size(X,1), size(X,2), rank_x, cond(task_Z), ...
    max_r, max(vif), qa_pass, ...
    'VariableNames', {'hemoglobin','n_scans','n_columns','rank', ...
    'normalized_condition_number','max_abs_correlation','max_vif','qa_pass'});
end
