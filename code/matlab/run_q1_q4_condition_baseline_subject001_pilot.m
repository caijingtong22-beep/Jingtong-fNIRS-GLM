function result = run_q1_q4_condition_baseline_subject001_pilot(question_key, phase_key)
% Explicit ExperimentalBaseline pilot for one Q1-Q4 condition model.
% Target-event rule: retain only the first occurrence (smallest onset) of
% each condition/phase label. The continuous recording is never truncated.
if nargin < 1, question_key = "q1"; end
if nargin < 2, phase_key = "procedure"; end
question_key = lower(string(question_key));
phase_key = lower(string(phase_key));
assert(any(question_key == ["q1","q2","q3","q4"]), 'Pilot supports Q1-Q4.');
assert(any(phase_key == ["procedure", "dom", "aoe"]), 'This pilot supports Procedure, DOM or AOE.');

sid = "subject001";
repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
glm = 'C:\aliceproject\GLM';
first_level = fullfile(glm, 'first_level');
input_dir = fullfile(first_level, 'nirs_spm_inputs', char(sid));
input_file = fullfile(input_dir, char(sid + "_nirs_data.mat"));
event_file = fullfile(glm, 'amended_snirf', 'event_tables_full_annotations', ...
    'event_table_all_subjects_full_annotations.csv');
condition_key="lile"; if question_key=="q2", condition_key="hile"; elseif question_key=="q3", condition_key="lihe"; elseif question_key=="q4", condition_key="hihe"; end
root = fullfile(repo, 'results', char(question_key+"_"+condition_key+"_vs_experimentalbaseline_"+phase_key), ...
    'subject001_pilot');
if ~isfolder(root), mkdir(root); end

addpath(first_level);
addpath(fullfile(glm, 'toolboxes', 'spm'));
addpath(genpath(fullfile(glm, 'toolboxes', 'NIRS_SPM_v4_r1')));
assert(isfile(input_file), 'Continuous input is missing: %s', input_file);

target_label = "LILEProcedure";
if question_key=="q2", target_label="HILEProcedure"; elseif question_key=="q3", target_label="LIHEProcedure"; elseif question_key=="q4", target_label="HIHEProcedure"; end
if phase_key == "dom", target_label = "DivisionOfMesoappendixLILE"; end
if phase_key == "aoe", target_label = "ApplicationOfEndoloopsLILE"; end
if question_key=="q2" && phase_key=="dom", target_label="DivisionOfMesoappendixHILE"; end
if question_key=="q2" && phase_key=="aoe", target_label="ApplicationOfEndoloopsHILE"; end
if question_key=="q3" && phase_key=="dom", target_label="DivisionOfMesoappendixLIHE"; end
if question_key=="q3" && phase_key=="aoe", target_label="ApplicationOfEndoloopsLIHE"; end
if question_key=="q4" && phase_key=="dom", target_label="DivisionOfMesoappendixHIHE"; end
if question_key=="q4" && phase_key=="aoe", target_label="ApplicationOfEndoloopsHIHE"; end
tag = question_key + "_" + phase_key;
events = readtable(event_file, 'TextType', 'string');
events.subject_id = compose('%03d', double(events.subject_id));
events.condition = string(events.condition);
sub_events = events(events.subject_id == "001", :);
baseline_rows = sortrows(sub_events(sub_events.condition == "ExperimentalBaseline", :), 'onset_sec');
target_rows = sortrows(sub_events(sub_events.condition == target_label, :), 'onset_sec');
assert(height(baseline_rows) == 1, 'Subject001 must have exactly one ExperimentalBaseline.');
assert(height(target_rows) >= 1, 'Subject001 has no %s event.', target_label);
target_row = target_rows(1, :); % User-confirmed first-occurrence policy.

provenance = [baseline_rows; target_row];
provenance.role = ["ExperimentalBaseline_boxcar"; "target_phase_HRF"];
provenance.selection_policy = ["single ExperimentalBaseline event"; ...
    "first occurrence by onset_sec (user-confirmed)"];
provenance.input_file = repmat(string(input_file), height(provenance), 1);
writetable(provenance, fullfile(root, char(sid + "_" + tag + "_selected_events.csv")));

D = load(input_file, 'nirs_data');
nscan = size(D.nirs_data.oxyData, 1);
fs = double(D.nirs_data.fs);
scan_time_sec = (0:nscan-1)' ./ fs;
baseline = scan_time_sec >= double(baseline_rows.onset_sec) & ...
    scan_time_sec < double(baseline_rows.onset_sec + baseline_rows.duration_sec);
assert(any(baseline) && any(~baseline), 'Baseline boxcar must contain both zeros and ones.');
boxcar_audit = table(repmat(sid, nscan, 1), (1:nscan)', scan_time_sec, double(baseline), ...
    'VariableNames', {'subject_id', 'scan_index', 'scan_time_sec', 'ExperimentalBaseline_boxcar'});
writetable(boxcar_audit, fullfile(root, char(sid + "_" + tag + "_baseline_boxcar_audit.csv")));

names = {char(target_label)};
onsets = {double(target_row.onset_sec)};
durations = {double(target_row.duration_sec)};
save(fullfile(root, char(sid + "_" + tag + "_target_condition.mat")), 'names', 'onsets', 'durations');

qa = table();
S_by_hb = cell(2, 1);
hb_labels = ["HbO", "HbR"];
for h = 1:2
    hb = hb_labels(h);
    hb_short = "hbr";
    if hb == "HbO", hb_short = "hbo"; end
    hb_root = fullfile(root, char(hb));
    if ~isfolder(hb_root), mkdir(hb_root); end

    % specification_batch receives only the HRF-convolved target condition.
    specification_batch(input_file, char(hb_short), 'DCT, 128', 'none', 0, ...
        hb_root, 0, 0, 1, names, onsets, durations);
    spm_file = fullfile(hb_root, sprintf('SPM_indiv_%s.mat', hb));
    L = load(spm_file, 'SPM_nirs');
    S = add_unconvolved_baseline(L.SPM_nirs, baseline, target_label);
    S.nirs.step = 'specification';
    SPM_nirs = S; %#ok<NASGU>
    save(spm_file, 'SPM_nirs');

    [qa_row, regressor_table] = design_qa(S, sid, hb);
    qa = [qa; qa_row]; %#ok<AGROW>
    writematrix(full(S.xX.X), fullfile(hb_root, char(sid + "_" + hb + "_" + tag + "_final_design.csv")));
    writetable(regressor_table, fullfile(hb_root, char(sid + "_" + hb + "_" + tag + "_regressor_table.csv")));
    S_by_hb{h} = S;
end

task_hbo = full(S_by_hb{1}.xX.X(:, 2));
task_hbr = full(S_by_hb{2}.xX.X(:, 2));
baseline_matches = isequal(full(S_by_hb{1}.xX.X(:, 1)), full(S_by_hb{2}.xX.X(:, 1)));
constant_matches = isequal(full(S_by_hb{1}.xX.X(:, 3)), full(S_by_hb{2}.xX.X(:, 3)));
inverse_hrf_matches = max(abs(task_hbr + task_hbo)) < 1e-10;
hb_check = table(baseline_matches, constant_matches, inverse_hrf_matches, ...
    'VariableNames', {'baseline_boxcar_matches', 'constant_matches', 'hbr_target_is_negative_hbo_target'});
writetable(hb_check, fullfile(root, char(sid + "_" + tag + "_hbo_hbr_design_check.csv")));
writetable(qa, fullfile(root, char(sid + "_" + tag + "_design_qa.csv")));
assert(all(qa.qa_pass), 'Design QA failed; estimation was not run.');
assert(all(hb_check{1,:}), 'HbO/HbR design convention check failed.');

contrast_summary = table();
for h = 1:2
    hb = hb_labels(h);
    hb_root = fullfile(root, char(hb));
    spm_file = fullfile(hb_root, sprintf('SPM_indiv_%s.mat', hb));
    current = load(spm_file, 'SPM_nirs');
    if ~isfield(current.SPM_nirs.nirs, 'step') || ~strcmpi(string(current.SPM_nirs.nirs.step), 'estimation')
        nirs_spm_estimation_batch_fixed(spm_file, input_file);
    end
    L = load(spm_file, 'SPM_nirs');
    S = L.SPM_nirs;
    contrast_summary = [contrast_summary; write_contrast(S, D.nirs_data, sid, hb, root, target_label, phase_key)]; %#ok<AGROW>
    plot_pilot_figures(S, D.nirs_data, sid, hb, root, phase_key, 1);
end
writetable(contrast_summary, fullfile(root, char(sid + "_" + tag + "_contrast_summary.csv")));
result = struct('output_root', root, 'design_qa', qa, 'contrast_summary', contrast_summary);
end

function S = add_unconvolved_baseline(S, baseline, target_label)
X0 = full(S.xX.X);
old_names = string(S.xX.name(:));
constant_index = find(contains(lower(old_names), 'constant'), 1);
target_index = find(contains(old_names, target_label), 1);
assert(numel(constant_index) == 1 && numel(target_index) == 1 && size(X0,2) == 2);
S.xX.X = [double(baseline(:)), X0(:, target_index), X0(:, constant_index)];
S.xX.name = {'Sn(1) ExperimentalBaseline_boxcar (unconvolved)', ...
    char(old_names(target_index)), 'Sn(1) constant'};
S.xX.iH = [];
S.xX.iC = [1 2];
S.xX.iB = 3;
S.xX.iG = [];
S.Sess(1).C.C = double(baseline(:));
S.Sess(1).C.name = {'ExperimentalBaseline_boxcar (unconvolved)'};
S.Sess(1).col = [1 2];
if isfield(S.xX, 'xKXs'), S.xX = rmfield(S.xX, 'xKXs'); end
if isfield(S.xX, 'pKX'), S.xX = rmfield(S.xX, 'pKX'); end
if isfield(S.xX, 'Bcov'), S.xX = rmfield(S.xX, 'Bcov'); end
end

function [row, regressor_table] = design_qa(S, sid, hb)
X = full(S.xX.X);
names = string(S.xX.name(:));
assert(size(X, 2) == 3 && rank(X) == 3, 'Final design must have rank-3 and exactly three columns.');
baseline = X(:,1); task = X(:,2);
baseline_centered = baseline - mean(baseline);
task_centered = task - mean(task);
r = (baseline_centered' * task_centered) / sqrt((baseline_centered' * baseline_centered) * (task_centered' * task_centered));
Z = [ones(size(X,1),1), baseline];
r2 = 1 - sum((task - Z * (Z \ task)).^2) / sum((task - mean(task)).^2);
vif = 1 / max(1 - r2, eps);
pass = abs(r) < .95 && vif <= 10;
row = table(sid, hb, size(X,1), size(X,2), rank(X), r, vif, pass, ...
    'VariableNames', {'subject_id','hemoglobin','n_scans','n_columns','rank','task_baseline_correlation','task_vif','qa_pass'});
regressor_table = table((1:3)', names, [true; false; false], [false; true; false], [false; false; true], ...
    'VariableNames', {'column_index','regressor','is_unconvolved_baseline','is_hrf_convolved_target','is_constant'});
end

function summary = write_contrast(S, nirs_data, sid, hb, root, target_label, phase_key)
names = string(S.xX.name(:));
baseline_idx = find(contains(names, 'ExperimentalBaseline_boxcar'), 1);
target_idx = find(contains(names, target_label), 1);
assert(~isempty(baseline_idx) && ~isempty(target_idx), 'Could not match contrast regressors by name.');
c = zeros(numel(names),1); c(baseline_idx) = -1; c(target_idx) = 1;
B = double(S.nirs.beta);
nch = size(B,2);
assert(size(B,1) == 3 && nch == 61, 'Expected 3 beta rows and 61 channels.');
rss = double(S.nirs.ResSS); if ~isvector(rss), rss = diag(rss)'; end
dfree = double(S.xX.erdf);
beta = c' * B;
se = sqrt(max((c' * double(S.xX.Bcov) * c) .* (rss / dfree), 0));
t = beta ./ se;
p = betainc(dfree ./ (dfree + t.^2), dfree/2, .5);
q = bh(p);
metadata = readtable('C:\aliceproject\GLM\first_level\nirs_spm_inputs\subject001\subject001_channel_metadata.csv', 'TextType', 'string');
R = table(repmat(sid,nch,1), repmat(hb,nch,1), (1:nch)', string(metadata.channel_id), ...
    string(metadata.channel_type), string(metadata.roi), string(metadata.laterality), beta(:), se(:), t(:), p(:), q(:), ...
    repmat(join(string(c'), ' '),nch,1), 'VariableNames', {'subject_id','hemoglobin','channel_index','channel_id','channel_type','roi','laterality','contrast_beta','contrast_se','contrast_t','contrast_p_two_sided','contrast_q_bh61','contrast_vector'});
qkey="q1"; condition_key="lile"; if contains(target_label,"HILE"),qkey="q2";condition_key="hile";elseif contains(target_label,"LIHE"),qkey="q3";condition_key="lihe";elseif contains(target_label,"HIHE"),qkey="q4";condition_key="hihe";end
writetable(R, fullfile(root, char(sid + "_" + hb + "_"+qkey+"_"+condition_key+"_minus_experimentalbaseline_" + phase_key + "_channel_contrasts.csv")));
summary = table(hb, nch, dfree, mean(beta), median(beta), sum(p < .05), sum(q < .05), ...
    'VariableNames', {'hemoglobin','n_channels','effective_df','mean_contrast_beta','median_contrast_beta','n_p_lt_005','n_q_bh61_lt_005'});
end

function plot_pilot_figures(S, nirs_data, sid, hb, root, phase_key, channel_id)
out = fullfile(root, char(hb), 'figures'); if ~isfolder(out), mkdir(out); end
set(0, 'DefaultFigureVisible', 'off');
set(groot, 'defaultAxesXColor', 'k', 'defaultAxesYColor', 'k', 'defaultTextColor', 'k');
X = full(S.xX.X); names = string(S.xX.name(:)); t = (0:size(X,1)-1)' * double(S.xY.RT);
h = (0:size(S.xBF.bf,1)-1)' * double(S.xBF.dt);
f = figure('Color','w'); plot(h,S.xBF.bf(:,1),'LineWidth',1.5); grid on; xlabel('Time after onset (s)'); ylabel('HRF amplitude'); title(sprintf('Q1 %s canonical HRF (%s)',upper(phase_key),hb));
exportgraphics(f,fullfile(out,char(sid+"_"+hb+"_q1_"+phase_key+"_canonical_hrf.png")),'Resolution',220); close(f);
f = figure('Color','w','Position',[100 100 1300 500]); imagesc(X'); axis tight; colorbar; yticks(1:3); yticklabels(names); xlabel('Scan'); ylabel('Final regressor'); title(sprintf('%s Q1 %s final design (%s): baseline is unconvolved boxcar',sid,upper(phase_key),hb));
exportgraphics(f,fullfile(out,char(sid+"_"+hb+"_q1_"+phase_key+"_final_design_heatmap.png")),'Resolution',220); close(f);
f = figure('Color','w','Position',[100 100 1400 500]); plot(t,X,'LineWidth',1); grid on; xlabel('Time (s)'); ylabel('Regressor amplitude'); legend(names,'Location','eastoutside'); title(sprintf('%s Q1 %s final design (%s)',sid,upper(phase_key),hb));
exportgraphics(f,fullfile(out,char(sid+"_"+hb+"_q1_"+phase_key+"_final_design_timeseries.png")),'Resolution',220); close(f);
field = 'dxyData'; if hb == "HbO", field = 'oxyData'; end
y = double(nirs_data.(field)(:,channel_id)); yf = spm_filter_HPF_LPF_WMDL(S.xX.K,y); fit = full(S.xX.xKXs.X) * double(S.nirs.beta(:,channel_id)); res = yf - fit;
f = figure('Color','w','Position',[100 100 1500 900]); tiledlayout(3,1,'TileSpacing','compact');
nexttile; plot(t,yf); hold on; plot(t,fit,'LineWidth',1.2); grid on; legend('Observed filtered signal','GLM fit','Location','eastoutside'); title(sprintf('%s channel %d %s: observed versus GLM fit',sid,channel_id,hb));
nexttile; plot(t,res); yline(0,'--k'); grid on; title('Residual');
nexttile; plot(t,X,'LineWidth',1); grid on; legend(names,'Location','eastoutside'); title('Final regressors: baseline boxcar is unconvolved'); xlabel('Time (s)');
exportgraphics(f,fullfile(out,char(sid+"_channel"+sprintf('%02d',channel_id)+"_"+hb+"_q1_"+phase_key+"_fit.png")),'Resolution',220); close(f);
end

function q = bh(p)
q = nan(size(p)); valid = isfinite(p); pv = p(valid); [s, order] = sort(pv(:)); m = numel(s);
adjusted = min(1, flipud(cummin(flipud(s .* m ./ (1:m)')))); out = zeros(size(pv)); out(order) = adjusted; q(valid) = out;
end
