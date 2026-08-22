function run_q4_aoe_rebuild()
% Rebuild Q4 AOE only: HIHE AOE minus explicit ExperimentalBaseline.
% This is intentionally independent of the legacy Q1-Q4 cohort script,
% which overwrote the AOE target label with LILE.

repo = 'C:\aliceproject\Jingtong-fNIRS-GLM';
glm = 'C:\aliceproject\GLM';
fl = fullfile(glm, 'first_level');
root = fullfile(repo, 'results', 'q4_hihe_vs_experimentalbaseline_aoe');
cohort = fullfile(root, 'cohort_specification');
event_file = fullfile(glm, 'amended_snirf', 'event_tables_full_annotations', ...
    'event_table_all_subjects_full_annotations.csv');
manifest_file = fullfile(repo, 'results', 'q1_lile_vs_experimentalbaseline_procedure', ...
    'cohort_specification', 'q1_procedure_explicit_baseline_inclusion_manifest.csv');
target_label = "ApplicationOfEndoloopsHIHE";
phase_key = "aoe";
question_key = "q4";
condition_key = "hihe";

addpath(fl);
addpath(fullfile(glm, 'toolboxes', 'spm'));
addpath(genpath(fullfile(glm, 'toolboxes', 'NIRS_SPM_v4_r1')));
assert(isfile(event_file) && isfile(manifest_file), 'Required input/manifest missing.');
if ~isfolder(cohort), mkdir(cohort); end

% Re-run the already-established Subject001 pilot with the correct Q2 label.
run_q1_q4_condition_baseline_subject001_pilot(question_key, phase_key);

M = readtable(manifest_file, 'TextType', 'string');
M.subject_id = string(M.subject_id);
E = readtable(event_file, 'TextType', 'string');
E.subject_id = compose('%03d', double(E.subject_id));
E.condition = string(E.condition);
status = table();

for i = 1:height(M)
    sid = M.subject_id(i);
    if ~logical(M.included(i)) || sid == "subject001"
        continue
    end
    subject_root = fullfile(cohort, char(sid));
    if ~isfolder(subject_root), mkdir(subject_root); end
    input_file = fullfile(fl, 'nirs_spm_inputs', char(sid), char(sid + "_nirs_data.mat"));
    try
        D = load(input_file, 'nirs_data');
        nscan = size(D.nirs_data.oxyData, 1);
        fs = double(D.nirs_data.fs);
        event_id = extractAfter(sid, 'subject');
        sub = E(E.subject_id == event_id, :);
        baseline_row = sortrows(sub(sub.condition == "ExperimentalBaseline", :), 'onset_sec');
        target_rows = sortrows(sub(sub.condition == target_label, :), 'onset_sec');
        assert(height(baseline_row) == 1, '%s lacks unique ExperimentalBaseline.', sid);
        assert(height(target_rows) >= 1, '%s lacks %s.', sid, target_label);
        target_row = target_rows(1, :);

        scan_time = (0:nscan-1)' ./ fs;
        baseline = scan_time >= double(baseline_row.onset_sec) & ...
            scan_time < double(baseline_row.onset_sec + baseline_row.duration_sec);
        assert(any(baseline) && any(~baseline), '%s baseline boxcar has no variation.', sid);

        provenance = [baseline_row; target_row];
        provenance.role = ["ExperimentalBaseline_boxcar"; "target_phase_HRF"];
        provenance.selection_policy = ["single ExperimentalBaseline event"; "first occurrence by onset_sec"];
        provenance.input_file = repmat(string(input_file), 2, 1);
        writetable(provenance, fullfile(subject_root, char(sid + "_q4_aoe_selected_events.csv")));
        writetable(table((1:nscan)', scan_time, double(baseline), ...
            'VariableNames', {'scan_index','scan_time_sec','ExperimentalBaseline_boxcar'}), ...
            fullfile(subject_root, char(sid + "_q4_aoe_baseline_boxcar_audit.csv")));

        qa = table();
        Sset = cell(2,1);
        for h = 1:2
            hb = ["HbO","HbR"];
            hb = hb(h);
            hb_short = "hbr"; if hb == "HbO", hb_short = "hbo"; end
            hb_root = fullfile(subject_root, char(hb));
            if ~isfolder(hb_root), mkdir(hb_root); end
            specification_batch(input_file, char(hb_short), 'DCT, 128', 'none', 0, ...
                hb_root, 0, 0, 1, {char(target_label)}, {double(target_row.onset_sec)}, {double(target_row.duration_sec)});
            spm_file = fullfile(hb_root, sprintf('SPM_indiv_%s.mat', hb));
            L = load(spm_file, 'SPM_nirs');
            S = add_baseline(L.SPM_nirs, baseline, target_label);
            S.nirs.step = 'specification';
            SPM_nirs = S; %#ok<NASGU>
            save(spm_file, 'SPM_nirs');
            [qa_row, regressor_table] = design_qa(S, sid, hb);
            qa = [qa; qa_row]; %#ok<AGROW>
            writematrix(full(S.xX.X), fullfile(hb_root, char(sid + "_" + hb + "_q4_aoe_final_design.csv")));
            writetable(regressor_table, fullfile(hb_root, char(sid + "_" + hb + "_q4_aoe_regressor_table.csv")));
            Sset{h} = S;
        end
        hb_ok = max(abs(full(Sset{1}.xX.X(:,2)) + full(Sset{2}.xX.X(:,2)))) < 1e-10 && ...
            isequal(full(Sset{1}.xX.X(:,1)), full(Sset{2}.xX.X(:,1))) && ...
            isequal(full(Sset{1}.xX.X(:,3)), full(Sset{2}.xX.X(:,3)));
        writetable(qa, fullfile(subject_root, char(sid + "_q4_aoe_design_qa.csv")));
        assert(all(qa.qa_pass) && hb_ok, '%s design QA/Hb convention failed.', sid);

        for h = 1:2
            hb = ["HbO","HbR"];
            hb = hb(h);
            spm_file = fullfile(subject_root, char(hb), sprintf('SPM_indiv_%s.mat', hb));
            nirs_spm_estimation_batch_fixed(spm_file, input_file);
            L = load(spm_file, 'SPM_nirs');
            write_contrast(L.SPM_nirs, sid, hb, subject_root, fl, target_label);
        end
        status = [status; table(sid, "completed", "", 'VariableNames', {'subject_id','status','detail'})]; %#ok<AGROW>
    catch ME
        status = [status; table(sid, "failed", string(ME.message), 'VariableNames', {'subject_id','status','detail'})]; %#ok<AGROW>
    end
    writetable(status, fullfile(cohort, 'q4_aoe_rebuild_estimation_manifest.csv'));
end

assert(all(status.status == "completed"), 'At least one Q2 AOE cohort subject failed.');
summarize_q1_procedure_explicit_baseline_group(phase_key, question_key);
plot_group_figures(cohort, question_key, condition_key, phase_key);
end

function S = add_baseline(S, baseline, target_label)
X = full(S.xX.X); names = string(S.xX.name(:));
ti = find(contains(names, target_label), 1);
ci = find(contains(lower(names), 'constant'), 1);
assert(size(X,2) == 2 && ~isempty(ti) && ~isempty(ci));
S.xX.X = [double(baseline(:)), X(:,ti), X(:,ci)];
S.xX.name = {'Sn(1) ExperimentalBaseline_boxcar (unconvolved)', ...
    char("Sn(1) " + target_label + "*bf(1)"), 'Sn(1) constant'};
S.xX.iH = []; S.xX.iC = [1 2]; S.xX.iB = 3; S.xX.iG = [];
S.Sess(1).C.C = double(baseline(:));
S.Sess(1).C.name = {'ExperimentalBaseline_boxcar (unconvolved)'};
S.Sess(1).col = [1 2];
for f = {'xKXs','pKX','Bcov'}
    if isfield(S.xX, f{1}), S.xX = rmfield(S.xX, f{1}); end
end
end

function [q, rows] = design_qa(S, sid, hb)
X = full(S.xX.X); b = X(:,1); t = X(:,2);
bc = b - mean(b); tc = t - mean(t);
r = (bc' * tc) / sqrt((bc' * bc) * (tc' * tc));
Z = [ones(size(X,1),1), b];
r2 = 1 - sum((t - Z*(Z\t)).^2) / sum((t - mean(t)).^2);
vif = 1 / max(1-r2, eps);
q = table(sid,hb,size(X,1),size(X,2),rank(X),r,vif,rank(X)==3 && abs(r)<.95 && vif<=10, ...
    'VariableNames', {'subject_id','hemoglobin','n_scans','n_columns','rank','task_baseline_correlation','task_vif','qa_pass'});
n = string(S.xX.name(:));
rows = table((1:3)',n,[true;false;false],[false;true;false],[false;false;true], ...
    'VariableNames', {'column_index','regressor','is_unconvolved_baseline','is_hrf_convolved_target','is_constant'});
end

function write_contrast(S, sid, hb, out, fl, target_label)
n = string(S.xX.name(:)); c = zeros(3,1);
c(find(contains(n,'ExperimentalBaseline_boxcar'),1)) = -1;
c(find(contains(n,target_label),1)) = 1;
assert(sum(c~=0)==2);
B = double(S.nirs.beta); rss = double(S.nirs.ResSS);
if ~isvector(rss), rss = diag(rss)'; end
df = double(S.xX.erdf); beta = c' * B;
se = sqrt(max((c' * double(S.xX.Bcov) * c) .* (rss / df), 0));
t = beta ./ se; p = betainc(df ./ (df+t.^2), df/2, .5); q = bh(p);
meta = readtable(fullfile(fl,'nirs_spm_inputs',char(sid),char(sid+'_channel_metadata.csv')), 'TextType','string');
nch = size(B,2);
T = table(repmat(sid,nch,1),repmat(hb,nch,1),(1:nch)',string(meta.channel_id), ...
    string(meta.channel_type),string(meta.roi),string(meta.laterality),beta(:),se(:),t(:),p(:),q(:), ...
    repmat(join(string(c'),' '),nch,1), 'VariableNames', ...
    {'subject_id','hemoglobin','channel_index','channel_id','channel_type','roi','laterality','contrast_beta','contrast_se','contrast_t','contrast_p_two_sided','contrast_q_bh61','contrast_vector'});
writetable(T, fullfile(out, char(sid+'_'+hb+'_q4_hihe_minus_experimentalbaseline_aoe_channel_contrasts.csv')));
end

function plot_group_figures(cohort, q, condition, phase)
I = readtable(fullfile(cohort, char(q+'_'+condition+'_minus_experimentalbaseline_'+phase+'_individual_ROI_long45_allQC_n29.csv')), 'TextType','string');
G = readtable(fullfile(cohort, char(q+'_'+condition+'_minus_experimentalbaseline_'+phase+'_group_ROI_long45_allQC_n29.csv')), 'TextType','string');
out = fullfile(cohort, 'figures_long45_allQC'); if ~isfolder(out), mkdir(out); end
keys = unique(I(I.hemoglobin=="HbO", {'roi','laterality'}), 'rows', 'stable');
f = figure('Visible','off','Color','w','Position',[100 100 1800 950]); tiledlayout(2,1);
for j=1:2
    hb = ["HbO","HbR"]; hb = hb(j); U = I(I.hemoglobin==hb,:); M = nan(29,height(keys));
    ids = unique(U.subject_id,'stable');
    for k=1:height(keys)
        for s=1:numel(ids)
            v=U(U.subject_id==ids(s)&U.roi==keys.roi(k)&U.laterality==keys.laterality(k),:);
            M(s,k)=v.roi_contrast_beta;
        end
    end
    nexttile; z=max(abs(M(:))); if z==0,z=1;end; imagesc(M,[-z z]); colormap(gca,parula); cb=colorbar;cb.Color='k';
    set(gca,'XTick',1:height(keys),'XTickLabel',strcat(keys.roi," | ",keys.laterality),'XTickLabelRotation',35,'YTick',1:numel(ids),'YTickLabel',ids,'XColor','k','YColor','k');
    title(hb+' individual ROI beta','Color','k'); xlabel('ROI | laterality'); ylabel('Subject');
end
sgtitle('Q4 AOE: HIHE - ExperimentalBaseline | individual 45-long-channel ROI beta (n=29)','Color','k');
exportgraphics(f,fullfile(out,'q4_aoe_layer1_individual_roi_heatmaps_long45_allQC_n29.png'),'Resolution',300);close(f);

f=figure('Visible','off','Color','w','Position',[100 100 1800 850]); tiledlayout(2,1);
for j=1:2
    hb=["HbO","HbR"];hb=hb(j); U=G(G.hemoglobin==hb,:); nexttile;
    bar(U.mean_contrast_beta,'FaceColor',[0.2 0.45 0.7]); hold on; errorbar(1:height(U),U.mean_contrast_beta,U.se_contrast_beta,'k.','LineWidth',1);
    ix=U.p_two_sided<.05; plot(find(ix),U.mean_contrast_beta(ix),'ko','MarkerFaceColor','k'); yline(0,'k-'); grid on;
    set(gca,'XTick',1:height(U),'XTickLabel',strcat(U.roi," | ",U.laterality),'XTickLabelRotation',35,'XColor','k','YColor','k');
    title(hb+' group ROI beta; black dots = uncorrected p < .05','Color','k'); ylabel('Mean contrast beta');
end
sgtitle('Q4 AOE: HIHE - ExperimentalBaseline | group 45-long-channel ROI effects (n=29)','Color','k');
exportgraphics(f,fullfile(out,'q4_aoe_layer2_group_roi_effects_long45_allQC_n29.png'),'Resolution',300);close(f);
end

function q = bh(p)
q = nan(size(p)); valid=isfinite(p); pv=p(valid); [s,o]=sort(pv(:)); m=numel(s);
a=min(1,flipud(cummin(flipud(s.*m./(1:m)')))); z=zeros(size(pv)); z(o)=a; q(valid)=z;
end
