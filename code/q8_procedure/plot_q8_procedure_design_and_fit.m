function outputs = plot_q8_procedure_design_and_fit(subject_num, channel_id, hb_label, phase_key)
% Plot Q8 Procedure's canonical HRF, convolved-regressor heatmap, and,
% after estimation, high-pass-filtered observed/fitted/residual traces.
%
% Examples:
%   plot_q8_procedure_design_and_fit                 % subject001, ch1, HbO
%   plot_q8_procedure_design_and_fit(1, 1, "HbR")
%
% The fit panel intentionally requires an estimated SPM_indiv file. Before
% estimation this function still exports the HRF and design figures, but
% does not invent a fitted signal from unestimated coefficients.

if nargin < 1, subject_num = 1; end
if nargin < 2, channel_id = 1; end
if nargin < 3, hb_label = "HbO"; end
if nargin < 4, phase_key = "procedure"; end
hb_label = string(hb_label);
phase_key = lower(string(phase_key));
assert(ismember(phase_key, ["procedure", "dom", "aoe"]), 'phase_key must be procedure, dom, or aoe.');
phase_display = upper(phase_key);
assert(ismember(hb_label, ["HbO", "HbR"]), 'hb_label must be HbO or HbR.');
assert(isscalar(subject_num) && subject_num >= 1 && subject_num <= 30, ...
    'subject_num must be in 1:30.');

this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(fileparts(this_file)));
glm_root = 'C:\\aliceproject\\GLM';
first_level = fullfile(glm_root, 'first_level');
subject_label = sprintf('subject%03d', subject_num);
spec_root = fullfile(repo_root, 'results', sprintf('q8_hihe_vs_lile_%s', phase_key), ...
    'cohort_specification', subject_label, char(hb_label));
spm_file = fullfile(spec_root, sprintf('SPM_indiv_%s.mat', hb_label));
assert(isfile(spm_file), 'Missing Q8 %s specification: %s', phase_display, spm_file);
output_dir = fullfile(spec_root, 'figures');
if ~isfolder(output_dir), mkdir(output_dir); end

addpath(first_level);
addpath(fullfile(glm_root, 'toolboxes', 'spm'));
addpath(genpath(fullfile(glm_root, 'toolboxes', 'NIRS_SPM_v4_r1')));
loaded = load(spm_file, 'SPM_nirs'); SPM_nirs = loaded.SPM_nirs;
X = full(SPM_nirs.xX.X);
regressor_names = string(SPM_nirs.xX.name(:));
task_columns = ~contains(lower(regressor_names), 'constant');
assert(sum(task_columns) == 2, 'Expected exactly two Q8 task regressors.');

% 1. Canonical NIRS-SPM HRF used by specification_batch.
hrf_time = (0:size(SPM_nirs.xBF.bf,1)-1)' * double(SPM_nirs.xBF.dt);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 520]);
set(fig, 'DefaultAxesColor', 'w', 'DefaultAxesXColor', 'k', ...
    'DefaultAxesYColor', 'k', 'DefaultTextColor', 'k');
plot(hrf_time, SPM_nirs.xBF.bf(:,1), 'LineWidth', 1.8, 'Color', [0 0.35 0.7]);
xlabel('Time after event onset (s)'); ylabel('HRF amplitude'); grid on;
title(sprintf('Q8 %s canonical HRF (%s)', phase_display, hb_label),'Color','k');
exportgraphics(fig, fullfile(output_dir, sprintf('%s_%s_q8_%s_canonical_hrf.png', ...
    subject_label, hb_label, phase_key)), 'Resolution', 220);
close(fig);

% 2. The task columns of X are the HRF-convolved regressors.  The constant
% is excluded so that colour is informative for condition timing/amplitude.
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1350 520]);
set(fig, 'DefaultAxesColor', 'w', 'DefaultAxesXColor', 'k', ...
    'DefaultAxesYColor', 'k', 'DefaultTextColor', 'k');
imagesc(X(:,task_columns)'); axis tight;
colormap(fig, parula); cb=colorbar; cb.Color='k';
yticks(1:sum(task_columns)); yticklabels(regressor_names(task_columns));
xlabel('Scan'); ylabel('HRF-convolved task regressor');
title(sprintf('%s Q8 %s: HRF-convolved regressor heatmap (%s)', subject_label, phase_display, hb_label),'Color','k');
exportgraphics(fig, fullfile(output_dir, sprintf('%s_%s_q8_%s_convolved_regressor_heatmap.png', ...
    subject_label, hb_label, phase_key)), 'Resolution', 220);
close(fig);

% 3. Line-view equivalent of the task-regressor panel in the reference
% figure: the two convolved task columns, deliberately excluding constant.
time_sec = (0:size(X,1)-1)' / double(data_sampling_rate(SPM_nirs));
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 600]);
set(fig, 'DefaultAxesColor', 'w', 'DefaultAxesXColor', 'k', ...
    'DefaultAxesYColor', 'k', 'DefaultTextColor', 'k');
plot(time_sec, X(:,task_columns), 'LineWidth', 1.1); grid on;
xlabel('Time (s)'); ylabel('Task regressor');
title(sprintf('%s Q8 %s: two HRF-convolved task regressors (%s)', subject_label, phase_display, hb_label),'Color','k');
legend(regressor_names(task_columns), 'Location', 'eastoutside');
line_png = fullfile(output_dir, sprintf('%s_%s_q8_%s_convolved_regressor_timeseries.png', ...
    subject_label, hb_label, phase_key));
exportgraphics(fig, line_png, 'Resolution', 220); close(fig);

outputs = struct();
outputs.hrf_png = fullfile(output_dir, sprintf('%s_%s_q8_%s_canonical_hrf.png', subject_label, hb_label, phase_key));
outputs.heatmap_png = fullfile(output_dir, sprintf('%s_%s_q8_%s_convolved_regressor_heatmap.png', subject_label, hb_label, phase_key));
outputs.timeseries_png = line_png;
outputs.fit_png = "";

% 4. Fit diagnostics are valid only after NIRS-SPM estimation has supplied
% beta and the exact DCT high-pass filtered design/data used in estimation.
if ~isfield(SPM_nirs.nirs, 'beta') || ~isfield(SPM_nirs.xX, 'xKXs')
    warning('Q8 %s has not been estimated; exported HRF and heatmap only.', phase_display);
    return;
end
n_channels = size(SPM_nirs.nirs.beta, 2);
assert(channel_id >= 1 && channel_id <= n_channels && mod(channel_id,1) == 0, ...
    'channel_id must be an integer in 1:%d.', n_channels);
if hb_label == "HbO", data_field = 'oxyData'; else, data_field = 'dxyData'; end
if ismember(subject_num, [1, 2, 10])
    data_file = fullfile(first_level, 'nirs_spm_inputs_v2_fulltimeline_trial1', ...
        subject_label, sprintf('%s_nirs_data.mat', subject_label));
else
    data_file = fullfile(first_level, 'nirs_spm_inputs', subject_label, ...
        sprintf('%s_nirs_data.mat', subject_label));
end
data = load(data_file, 'nirs_data');
y = double(data.nirs_data.(data_field)(:, channel_id));
t = double(data.nirs_data.time(:));
X_fit = full(SPM_nirs.xX.xKXs.X);
y_fit = spm_filter_HPF_LPF_WMDL(SPM_nirs.xX.K, y);
beta = double(SPM_nirs.nirs.beta(:, channel_id));
fitted = X_fit * beta; residual = y_fit - fitted;
assert(numel(y_fit) == size(X_fit,1), 'Filtered data/design length mismatch.');

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1550 880]);
set(fig, 'DefaultAxesColor', 'w', 'DefaultAxesXColor', 'k', ...
    'DefaultAxesYColor', 'k', 'DefaultTextColor', 'k');
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, y_fit, 'LineWidth', 0.7, 'Color', [0 0.35 0.7]); hold on;
plot(t, fitted, 'LineWidth', 1.2, 'Color', [0.85 0.33 0.1]); grid on;
xlabel('Time (s)'); ylabel('High-pass filtered concentration');
title(sprintf('%s, channel %d, %s: observed and GLM fitted signal', subject_label, channel_id, hb_label));
legend({'Observed (DCT high-pass filtered)', 'GLM fitted'}, 'Location', 'eastoutside');
nexttile;
plot(t, residual, 'LineWidth', 0.7, 'Color', [0 0.35 0.7]); hold on; yline(0, '--k'); grid on;
xlabel('Time (s)'); ylabel('Residual'); title('GLM residual (estimation space)');
nexttile;
plot(t, X(:,task_columns), 'LineWidth', 1.0); grid on;
xlabel('Time (s)'); ylabel('Task regressor');
title(sprintf('Q8 %s HRF-convolved task regressors', phase_display),'Color','k');
legend(regressor_names(task_columns), 'Location', 'eastoutside');
fit_png = fullfile(output_dir, sprintf('%s_channel%02d_%s_q8_%s_fit.png', ...
    subject_label, channel_id, hb_label, phase_key));
exportgraphics(fig, fit_png, 'Resolution', 220); close(fig);
outputs.fit_png = fit_png;
end

function fs = data_sampling_rate(SPM_nirs)
% SPM_nirs stores scan timing as RT.  This avoids loading raw data when only
% design visualisations are requested.
assert(isfield(SPM_nirs, 'xY') && isfield(SPM_nirs.xY, 'RT') && SPM_nirs.xY.RT > 0, ...
    'SPM_nirs does not contain a positive scan interval (xY.RT).');
fs = 1 / double(SPM_nirs.xY.RT);
end
