function [log_rows, outputs] = run_q10_subject001_validated_spatial_maps()
%RUN_Q10_SUBJECT001_VALIDATED_SPATIAL_MAPS Render Q10 HE-LE subject001 maps.
% Uses the audited 45-long-channel registration.  These are individual-level,
% uncorrected visualisations, not a 30-subject second-level inference.

this_file = mfilename('fullpath');
model_root = fileparts(this_file);
first_level = fileparts(model_root);
glm_root = fileparts(first_level);
nirs_spm_root = fullfile(glm_root, 'toolboxes', 'NIRS_SPM_v4_r1');
spm_root = fullfile(glm_root, 'toolboxes', 'spm');
registration_file = fullfile(first_level, 'group_outputs', 'NIRS_SPM', 'channel_registration_q10_validated', 'q10_validated_channel_registration.mat');
q10_root = fullfile(model_root, 'q10_hiloadhe_vs_hiloadle');
outputs.root = fullfile(q10_root, 'spatial_maps_subject001_validated');
outputs.figure_dir = fullfile(outputs.root, 'figures');
outputs.log_csv = fullfile(outputs.root, 'q10_subject001_spatial_map_log.csv');
outputs.summary_md = fullfile(outputs.root, 'Q10_SUBJECT001_SPATIAL_MAP_SUMMARY.md');
outputs.runtime_dir = fullfile(outputs.root, 'runtime_assets');
if ~isfolder(outputs.root), mkdir(outputs.root); end
if ~isfolder(outputs.figure_dir), mkdir(outputs.figure_dir); end
if ~isfolder(outputs.runtime_dir), mkdir(outputs.runtime_dir); end

assert(isfile(registration_file), 'Build the validated registration first: %s', registration_file);
addpath(spm_root); addpath(genpath(nirs_spm_root));
ensure_split_colormap(outputs.runtime_dir); addpath(outputs.runtime_dir);

views = ["dorsal", "frontal"];
hb_labels = ["HbO", "HbR"];
log_rows = table();
for hb = hb_labels
    hb_char = char(hb);
    source_spm = fullfile(q10_root, 'subject001', hb_char, sprintf('SPM_indiv_%s.mat', hb_char));
    assert(isfile(source_spm), 'Missing Q10 estimated SPM: %s', source_spm);
    loaded = load(source_spm, 'SPM_nirs'); SPM_nirs = loaded.SPM_nirs;
    assert(strcmpi(string(SPM_nirs.nirs.step), 'estimation'), '%s is not estimated.', hb_char);
    names = string(SPM_nirs.xX.name(:));
    idx_le = find(contains(names, 'HILoadLE'), 1); idx_he = find(contains(names, 'HILoadHE'), 1);
    assert(~isempty(idx_le) && ~isempty(idx_he), 'Cannot find Q10 task regressors in %s.', hb_char);
    c = zeros(numel(names),1); c(idx_le) = -1; c(idx_he) = 1;
    assert(numel(c) == 3, 'Q10 subject001 must have 3 regressors; found %d.', numel(c));
    for view = views
        work_dir = fullfile(outputs.root, 'subject001', hb_char, char(view));
        if ~isfolder(work_dir), mkdir(work_dir); end
        work_spm = fullfile(work_dir, sprintf('SPM_indiv_%s.mat', hb_char));
        copyfile(source_spm, work_spm, 'f');
        prefix = sprintf('sub-001_%s_Q10_HILoadHE_minus_HILoadLE_%s', hb_char, view);
        out_mat = fullfile(work_dir, [prefix '_map.mat']);
        out_png = fullfile(outputs.figure_dir, [prefix '_unthresholded.png']);
        status = "failed"; message = ""; threshold = NaN; timer = tic;
        try
            [stat_brain, act_brain, threshold] = activation_map_batch(work_spm, registration_file, ...
                'Q10_HILoadHE_minus_HILoadLE', c, 'T', char(view), 0.05, 'none', 0);
            save(out_mat, 'stat_brain', 'act_brain', 'threshold', 'c', 'names', 'source_spm', 'registration_file', '-v7.3');
            save_map_figure(stat_brain, out_png, sprintf('subject001 %s Q10: HILoadHE - HILoadLE (%s)', hb_char, view));
            status = "success";
        catch ME
            message = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
        end
        row = table("subject001", hb, view, "Q10_HILoadHE_minus_HILoadLE", ...
            string(strjoin(string(c'), ' ')), status, toc(timer), threshold, string(out_mat), string(out_png), message, ...
            'VariableNames', {'subject_id','hemoglobin','view','contrast_name','contrast_vector', ...
            'status','elapsed_sec','threshold','map_mat','figure_png','message'});
        log_rows = [log_rows; row]; %#ok<AGROW>
        writetable(log_rows, outputs.log_csv);
    end
end
write_summary(outputs.summary_md, outputs, log_rows, registration_file);
assert(all(log_rows.status == "success"), 'At least one Q10 spatial map failed; inspect %s.', outputs.log_csv);
end

function ensure_split_colormap(asset_root)
filename = fullfile(asset_root, 'Split.mat');
if ~isfile(filename), split = jet(128); save(filename, 'split'); end %#ok<NASGU>
end

function save_map_figure(stat_brain, filename, title_text)
fig = figure('Visible', 'off', 'Color', 'w'); imagesc(stat_brain); axis image off; colormap(jet); colorbar;
title(title_text, 'Interpreter', 'none'); set(fig, 'PaperPositionMode', 'auto'); print(fig, filename, '-dpng', '-r300'); close(fig);
end

function write_summary(filename, outputs, log_rows, registration_file)
fid = fopen(filename, 'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# Q10 subject001 validated spatial maps\n\n');
fprintf(fid, '- Registration: `%s`\n- Contrast: `HILoadHE - HILoadLE`\n- Implemented regressors: `[HILoadLE, HILoadHE, constant]`\n- Implemented contrast: `[-1, +1, 0]`\n- Views: dorsal and frontal\n- Threshold configuration: uncorrected `p=0.05` for individual-level visualisation\n\n', registration_file);
fprintf(fid, '## Status\n\n');
for k = 1:height(log_rows), fprintf(fid, '- %s %s: %s\n', log_rows.hemoglobin(k), log_rows.view(k), log_rows.status(k)); end
fprintf(fid, '\n## Interpretation boundary\n\nThese are first-level subject001 maps. They validate channel registration and visualisation, but are not 30-subject group-level NIRS-SPM inference.\n');
end
