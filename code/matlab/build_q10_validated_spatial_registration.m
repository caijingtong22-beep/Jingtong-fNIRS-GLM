function outputs = build_q10_validated_spatial_registration()
%BUILD_Q10_VALIDATED_SPATIAL_REGISTRATION Build an audited Q10 map registration.
% This does not alter channel_allocation.csv.  It uses the VC source/detector
% names for the 45 long channels, with one audited correction:
% channel 40: Src11-Det22 -> Src11-Det06.  The latter is the SNIRF pair and
% follows the raw-to-VC detector conversion used throughout the allocation.

this_file = mfilename('fullpath');
output_root = fileparts(this_file);
first_level = fileparts(fileparts(output_root));
glm_root = fileparts(first_level);

allocation_file = fullfile(glm_root, 'channel_allocation.csv');
snirf_file = fullfile(glm_root, 'amended_snirf', 'Subj001_Sess101_run001.snirf');
probe_file = fullfile(first_level, 'nirs_spm_inputs', 'subject001', 'subject001_probe_geometry.mat');
spm_root = fullfile(glm_root, 'toolboxes', 'spm');
nirs_spm_root = fullfile(glm_root, 'toolboxes', 'NIRS_SPM_v4_r1');
out_dir = fullfile(output_root, 'channel_registration_q10_validated');
if ~isfolder(out_dir), mkdir(out_dir); end

outputs.mapping_csv = fullfile(out_dir, 'q10_validated_channel_pair_mapping.csv');
outputs.coordinate_csv = fullfile(out_dir, 'q10_validated_channel_mni_coordinates.csv');
outputs.registration_mat = fullfile(out_dir, 'q10_validated_channel_registration.mat');
outputs.report_md = fullfile(out_dir, 'Q10_SPATIAL_REGISTRATION_REPORT.md');

assert(isfile(allocation_file), 'Missing allocation: %s', allocation_file);
assert(isfile(snirf_file), 'Missing SNIRF: %s', snirf_file);
assert(isfile(probe_file), 'Missing probe geometry: %s', probe_file);
addpath(spm_root);
addpath(genpath(nirs_spm_root));

allocation = readtable(allocation_file, 'TextType', 'string', 'VariableNamingRule', 'preserve');
channel_id = double(allocation.('SamplingSite.id'));
channel_name = string(allocation.('SamplingSite.name'));
vc_name = string(allocation.('VCSamplingSite.name'));
channel_type = string(allocation.('Channel'));
laterality = string(allocation.('Laterality'));
roi = string(allocation.('JC RoI Allocation'));
n_channels = height(allocation);
assert(n_channels == 61 && isequal(channel_id(:)', 1:61), 'Expected allocation channels 1:61.');
is_long = strcmpi(channel_type, 'Long');
assert(sum(is_long) == 45, 'Expected 45 long channels.');

source_index = nan(n_channels, 1);
detector_index = nan(n_channels, 1);
pair_text = strings(n_channels, 1);
mapping_rule = repmat("short_channel_not_rendered", n_channels, 1);
for ch = find(is_long)'
    pair_text(ch) = vc_name(ch);
    if channel_id(ch) == 40
        pair_text(ch) = "ch40:Src11-Det06";
        mapping_rule(ch) = "audited_correction_raw_Det22_to_VC_Det06";
    else
        mapping_rule(ch) = "VCSamplingSite_name";
    end
    [source_index(ch), detector_index(ch)] = parse_src_det(pair_text(ch));
end

snirf_pairs = read_snirf_pairs(snirf_file);
pair_present = false(n_channels, 1);
long_rows = find(is_long);
[pair_present(long_rows), ~] = ismember([source_index(long_rows), detector_index(long_rows)], snirf_pairs, 'rows');
assert(all(pair_present(long_rows)), 'At least one long-channel mapping pair is absent from the SNIRF measurement list.');
assert(size(unique([source_index(long_rows), detector_index(long_rows)], 'rows'), 1) == 45, ...
    'Long-channel mappings must be unique.');

loaded = load(probe_file, 'probe');
probe = loaded.probe;
required = {'sourcePos3D', 'detectorPos3D', 'sourcePos2D', 'detectorPos2D'};
for k = 1:numel(required), assert(isfield(probe, required{k}), 'Missing probe.%s.', required{k}); end

channel_mni_unprojected = ones(4, n_channels);
midpoint_2d = nan(n_channels, 2);
valid_geometry = false(n_channels, 1);
for ch = long_rows'
    src_xyz = double(probe.sourcePos3D(:, source_index(ch)));
    det_xyz = double(probe.detectorPos3D(:, detector_index(ch)));
    src_xy = double(probe.sourcePos2D(:, source_index(ch)));
    det_xy = double(probe.detectorPos2D(:, detector_index(ch)));
    channel_mni_unprojected(1:3, ch) = (src_xyz + det_xyz) ./ 2;
    midpoint_2d(ch, :) = ((src_xy + det_xy) ./ 2)';
    valid_geometry(ch) = all(isfinite(channel_mni_unprojected(1:3, ch)));
end
assert(all(valid_geometry(long_rows)), 'A long channel has invalid probe geometry.');

template_info = select_template_info(spm_root);
channel_mni_rendered = channel_mni_unprojected;
projection_used = false;
try
    channel_mni_rendered(:, long_rows) = projection_CS(channel_mni_rendered(:, long_rows));
    projection_used = true;
catch ME
    warning('projection_CS failed; retaining unprojected midpoints: %s', ME.message);
end

ch_mni_vx = inv(template_info.mat) * channel_mni_rendered(:, long_rows);
[rend, rendered_valid] = render_MNI_coordinates(ch_mni_vx, template_info);
rendered_MNI = cell(size(rendered_valid));
for view = 1:numel(rendered_valid)
    rendered_MNI{view}.rchn = -ones(n_channels, 1);
    rendered_MNI{view}.cchn = -ones(n_channels, 1);
    rendered_MNI{view}.rchn(long_rows) = rendered_valid{view}.rchn;
    rendered_MNI{view}.cchn(long_rows) = rendered_valid{view}.cchn;
    rendered_MNI{view}.ren = rend{view}.ren;
end

preproc_info = struct();
preproc_info.rend_ch_pos = rendered_MNI;
preproc_info.ch_MNI_mm = channel_mni_rendered;
preproc_info.ch_MNI_mm_unprojected = channel_mni_unprojected;
preproc_info.wT1_info = template_info;
preproc_info.channel_labels = compose("CH%d", channel_id);
preproc_info.source_index = source_index;
preproc_info.detector_index = detector_index;
preproc_info.valid_geometry = valid_geometry;
preproc_info.channel_type = channel_type;
preproc_info.laterality = laterality;
preproc_info.roi = roi;
preproc_info.mapping_rule = mapping_rule;
preproc_info.pair_present_in_subject001_snirf = pair_present;
preproc_info.source_file = struct('allocation', allocation_file, 'snirf', snirf_file, 'probe', probe_file, 'template', template_info.fname);
preproc_info.notes = [
    "Q10 spatial registration; independent of atlas selection."
    "Long-channel coordinates are source-detector midpoints from subject001 SNIRF probe geometry."
    "Channel 40 mapping corrected from unconverted Src11-Det22 to validated Src11-Det06."
    "All 45 long-channel source-detector pairs were checked against subject001 SNIRF measurement list."
    "Short channels are retained in model estimation but are not spatially rendered."
    ];
save(outputs.registration_mat, 'preproc_info');

mapping = table(channel_id, channel_name, vc_name, channel_type, laterality, roi, pair_text, ...
    source_index, detector_index, mapping_rule, pair_present, valid_geometry, ...
    'VariableNames', {'channel_id','sampling_site_name','vc_sampling_site_name','channel_type', ...
    'laterality','roi','validated_pair_name','source_index','detector_index', ...
    'mapping_rule','pair_present_in_subject001_snirf','rendered_geometry'});
writetable(mapping, outputs.mapping_csv);
coordinates = [mapping, table(channel_mni_unprojected(1,:)', channel_mni_unprojected(2,:)', channel_mni_unprojected(3,:)', ...
    channel_mni_rendered(1,:)', channel_mni_rendered(2,:)', channel_mni_rendered(3,:)', midpoint_2d(:,1), midpoint_2d(:,2), ...
    'VariableNames', {'mni_x_unprojected','mni_y_unprojected','mni_z_unprojected', ...
    'mni_x_rendered_input','mni_y_rendered_input','mni_z_rendered_input','midpoint_2d_x','midpoint_2d_y'})];
writetable(coordinates, outputs.coordinate_csv);
write_report(outputs.report_md, outputs, projection_used, template_info, mapping, snirf_pairs, allocation_file, snirf_file, probe_file);
fprintf('Q10 validated registration completed: 45/45 long-channel pairs verified.\n');
end

function [src, det] = parse_src_det(text)
token = regexp(text, 'Src([0-9]+)-Det([0-9]+)', 'tokens', 'once');
assert(~isempty(token), 'Cannot parse source/detector pair: %s', text);
src = str2double(token{1}); det = str2double(token{2});
end

function pairs = read_snirf_pairs(filename)
info = h5info(filename, '/nirs/data1');
n = numel(info.Groups); src = nan(n,1); det = nan(n,1);
for k = 1:n
    group = info.Groups(k).Name;
    src(k) = double(h5read(filename, [group '/sourceIndex']));
    det(k) = double(h5read(filename, [group '/detectorIndex']));
end
pairs = unique([src det], 'rows');
assert(size(pairs,1) == 61, 'Expected 61 unique SNIRF pairs, found %d.', size(pairs,1));
end

function info = select_template_info(spm_root)
candidates = {fullfile(spm_root, 'canonical', 'avg152T1.nii'), fullfile(spm_root, 'canonical', 'single_subj_T1.nii')};
for k = 1:numel(candidates)
    if isfile(candidates{k})
        V = spm_vol(candidates{k});
        info = struct('fname', candidates{k}, 'mat', V.mat, 'dim', V.dim);
        return;
    end
end
error('No usable SPM canonical T1 template found.');
end

function write_report(filename, outputs, projection_used, template_info, mapping, snirf_pairs, allocation_file, snirf_file, probe_file)
fid = fopen(filename, 'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# Q10 spatial registration report\n\n');
fprintf(fid, '## Purpose\n\nThis registration is for NIRS-SPM visualisation of Q10 and later small-model contrasts. It does not define the GLM or its ROI statistics.\n\n');
fprintf(fid, '## Inputs\n\n- Original allocation: `%s`\n- Subject001 SNIRF measurement list: `%s`\n- Subject001 SNIRF-derived probe geometry: `%s`\n- SPM template: `%s`\n\n', allocation_file, snirf_file, probe_file, template_info.fname);
fprintf(fid, '## Validated mapping\n\n- Total channels: %d\n- Long channels spatially rendered: %d\n- Short channels not rendered: %d\n- Unique SNIRF source-detector pairs: %d\n- Validated long-channel pairs present in SNIRF: %d / %d\n\n', height(mapping), sum(strcmpi(mapping.channel_type,'Long')), sum(strcmpi(mapping.channel_type,'Short')), size(snirf_pairs,1), sum(mapping.pair_present_in_subject001_snirf), sum(strcmpi(mapping.channel_type,'Long')));
fprintf(fid, '## Channel 40 audit correction\n\nThe original allocation has `ch40:Src11-Det22` in both source and VC fields. That pair is absent from the subject001 SNIRF measurement list. The allocation’s established conversion pattern maps raw `Det22` to VC/SNIRF `Det06`; SNIRF contains `Src11-Det06`. This registration therefore uses `ch40:Src11-Det06`, without altering the original allocation file.\n\n');
fprintf(fid, '## Coordinate method\n\nFor each validated long channel, the coordinate is the midpoint of source and detector `sourcePos3D` / `detectorPos3D` from the subject001 SNIRF-derived probe geometry. Coordinates are projected to the cortical surface with `projection_CS`: %d. NIRS-SPM then renders the channel positions on the SPM canonical template.\n\n', projection_used);
fprintf(fid, '## Limits\n\n- This is MNI-compatible montage registration, not an atlas-based cortical label assignment.\n- Short channels remain available to the GLM but are intentionally omitted from spatial rendering.\n- Q10 maps made with this file are subject-level unless a separately specified group-level model is fitted.\n');
end
