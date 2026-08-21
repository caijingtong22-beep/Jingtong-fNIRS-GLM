function run_q1_procedure_cohort_estimation(phase_key,question_key)
% Q1 Procedure explicit-baseline GLM for the 28 eligible subjects after
% the validated Subject001 pilot. Target policy: first LILEProcedure only.
if nargin < 1, phase_key = "procedure"; end
if nargin < 2, question_key="q1"; end
question_key=lower(string(question_key)); assert(any(question_key==["q1","q2","q3","q4"]));
phase_key = lower(string(phase_key)); assert(any(phase_key == ["procedure","dom","aoe"]));
target_label = "LILEProcedure"; if question_key=="q2",target_label="HILEProcedure";elseif question_key=="q3",target_label="LIHEProcedure";elseif question_key=="q4",target_label="HIHEProcedure";end; if phase_key == "dom", target_label="DivisionOfMesoappendixLILE"; if question_key=="q2",target_label="DivisionOfMesoappendixHILE";elseif question_key=="q3",target_label="DivisionOfMesoappendixLIHE";elseif question_key=="q4",target_label="DivisionOfMesoappendixHIHE";end;end
if phase_key=="aoe",target_label="ApplicationOfEndoloopsLILE";if question_key=="q2",target_label="ApplicationOfEndoloopsHILE";elseif question_key=="q3",target_label="ApplicationOfEndoloopsLIHE";elseif question_key=="q4",target_label="ApplicationOfEndoloopsHIHE";end;end
if phase_key == "aoe", target_label = "ApplicationOfEndoloopsLILE"; end
tag = question_key + "_" + phase_key;
repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
glm = 'C:\aliceproject\GLM'; fl = fullfile(glm, 'first_level');
condition_key="lile";if question_key=="q2",condition_key="hile";elseif question_key=="q3",condition_key="lihe";elseif question_key=="q4",condition_key="hihe";end
root = fullfile(repo,'results',char(question_key+"_"+condition_key+"_vs_experimentalbaseline_"+phase_key),'cohort_specification');
manifest_file = fullfile(repo, 'results', 'q1_lile_vs_experimentalbaseline_procedure', 'cohort_specification', 'q1_procedure_explicit_baseline_inclusion_manifest.csv');
event_file = fullfile(glm, 'amended_snirf', 'event_tables_full_annotations', 'event_table_all_subjects_full_annotations.csv');
addpath(fl); addpath(fullfile(glm, 'toolboxes', 'spm')); addpath(genpath(fullfile(glm, 'toolboxes', 'NIRS_SPM_v4_r1')));
M = readtable(manifest_file, 'TextType', 'string'); M.subject_id = string(M.subject_id);
E = readtable(event_file, 'TextType', 'string'); E.subject_id = compose('%03d', double(E.subject_id)); E.condition = string(E.condition);
status = table();
fprintf('Q2/Q1 cohort start: phase=%s question=%s root=%s\n',phase_key,question_key,root);
if ~isfolder(root), mkdir(root); end
for i = 1:height(M)
    sid = M.subject_id(i);
    if ~logical(M.included(i)) || sid == "subject001", continue; end
    subject_root = fullfile(root, char(sid)); if ~isfolder(subject_root), mkdir(subject_root); end
    fprintf('Starting %s\n',sid);
    input_file = fullfile(fl, 'nirs_spm_inputs', char(sid), char(sid + "_nirs_data.mat"));
    try
        D = load(input_file, 'nirs_data'); nscan = size(D.nirs_data.oxyData,1); fs = double(D.nirs_data.fs);
        event_id = extractAfter(sid, 'subject'); sub = E(E.subject_id == event_id,:);
        baseline_row = sortrows(sub(sub.condition == "ExperimentalBaseline",:),'onset_sec');
        target_row = sortrows(sub(sub.condition == target_label,:),'onset_sec');
        assert(height(baseline_row) == 1 && height(target_row) >= 1, 'Event availability changed after audit.');
        target_row = target_row(1,:);
        scan_time = (0:nscan-1)' ./ fs;
        baseline = scan_time >= double(baseline_row.onset_sec) & scan_time < double(baseline_row.onset_sec + baseline_row.duration_sec);
        assert(any(baseline) && any(~baseline), 'Baseline boxcar has no variation.');
        provenance = [baseline_row; target_row]; provenance.role = ["ExperimentalBaseline_boxcar";"target_phase_HRF"];
        provenance.selection_policy = ["single ExperimentalBaseline event";"first occurrence by onset_sec (user-confirmed)"];
        provenance.input_file = repmat(string(input_file),2,1);
        writetable(provenance, fullfile(subject_root, char(sid + "_" + tag + "_selected_events.csv")));
        writetable(table((1:nscan)',scan_time,double(baseline),'VariableNames',{'scan_index','scan_time_sec','ExperimentalBaseline_boxcar'}), fullfile(subject_root, char(sid + "_q1_procedure_baseline_boxcar_audit.csv")));
        qa = table(); Sset = cell(2,1); hbs = ["HbO","HbR"];
        for h = 1:2
            hb = hbs(h); short = "hbr"; if hb == "HbO", short = "hbo"; end
            hb_root = fullfile(subject_root,char(hb)); if ~isfolder(hb_root), mkdir(hb_root); end
            specification_batch(input_file,char(short),'DCT, 128','none',0,hb_root,0,0,1,{char(target_label)},{double(target_row.onset_sec)},{double(target_row.duration_sec)});
            spm_file = fullfile(hb_root,sprintf('SPM_indiv_%s.mat',hb)); L = load(spm_file,'SPM_nirs'); S = add_baseline(L.SPM_nirs,baseline,target_label);
            S.nirs.step = 'specification'; SPM_nirs = S; save(spm_file,'SPM_nirs');
            [q,row] = qa_design(S,sid,hb); qa=[qa;q]; %#ok<AGROW>
            writematrix(full(S.xX.X),fullfile(hb_root,char(sid+"_"+hb+"_"+tag+"_final_design.csv")));
            writetable(row,fullfile(hb_root,char(sid+"_"+hb+"_q1_procedure_regressor_table.csv"))); Sset{h}=S;
        end
        hbr_check = max(abs(full(Sset{1}.xX.X(:,2)) + full(Sset{2}.xX.X(:,2)))) < 1e-10 && isequal(full(Sset{1}.xX.X(:,1)),full(Sset{2}.xX.X(:,1))) && isequal(full(Sset{1}.xX.X(:,3)),full(Sset{2}.xX.X(:,3)));
        writetable(qa,fullfile(subject_root,char(sid+"_"+tag+"_design_qa.csv"))); assert(all(qa.qa_pass) && hbr_check,'Design QA failed.');
        for h = 1:2
            hb=hbs(h); hb_root=fullfile(subject_root,char(hb)); spm_file=fullfile(hb_root,sprintf('SPM_indiv_%s.mat',hb));
            nirs_spm_estimation_batch_fixed(spm_file,input_file); L=load(spm_file,'SPM_nirs'); write_contrast(L.SPM_nirs,sid,hb,subject_root,fl,target_label,phase_key);
        end
        status=[status;table(sid,"completed","",'VariableNames',{'subject_id','status','detail'})]; %#ok<AGROW>
    catch ME
        status=[status;table(sid,"failed",string(ME.message),'VariableNames',{'subject_id','status','detail'})]; %#ok<AGROW>
    end
    writetable(status,fullfile(root,'q1_procedure_estimation_manifest_remaining_subjects.csv'));
end
writetable(status,fullfile(root,char(tag+"_estimation_manifest_remaining_subjects.csv")));
end

function S=add_baseline(S,baseline,target_label)
X=full(S.xX.X); names=string(S.xX.name(:)); ti=find(contains(names,target_label),1); ci=find(contains(lower(names),'constant'),1);
assert(size(X,2)==2 && ~isempty(ti) && ~isempty(ci)); S.xX.X=[double(baseline(:)),X(:,ti),X(:,ci)];
S.xX.name={'Sn(1) ExperimentalBaseline_boxcar (unconvolved)',char("Sn(1) "+target_label+"*bf(1)"),'Sn(1) constant'};
S.xX.iH=[]; S.xX.iC=[1 2]; S.xX.iB=3; S.xX.iG=[]; S.Sess(1).C.C=double(baseline(:)); S.Sess(1).C.name={'ExperimentalBaseline_boxcar (unconvolved)'}; S.Sess(1).col=[1 2];
for f={'xKXs','pKX','Bcov'}, if isfield(S.xX,f{1}), S.xX=rmfield(S.xX,f{1}); end, end
end

function [q,row]=qa_design(S,sid,hb)
X=full(S.xX.X); b=X(:,1); t=X(:,2); bc=b-mean(b); tc=t-mean(t); r=(bc'*tc)/sqrt((bc'*bc)*(tc'*tc)); Z=[ones(size(X,1),1),b]; r2=1-sum((t-Z*(Z\t)).^2)/sum((t-mean(t)).^2); vif=1/max(1-r2,eps);
q=table(sid,hb,size(X,1),size(X,2),rank(X),r,vif,rank(X)==3 && abs(r)<.95 && vif<=10,'VariableNames',{'subject_id','hemoglobin','n_scans','n_columns','rank','task_baseline_correlation','task_vif','qa_pass'});
n=string(S.xX.name(:)); row=table((1:3)',n,[true;false;false],[false;true;false],[false;false;true],'VariableNames',{'column_index','regressor','is_unconvolved_baseline','is_hrf_convolved_target','is_constant'});
end

function write_contrast(S,sid,hb,out,fl,target_label,phase_key)
n=string(S.xX.name(:)); c=zeros(3,1); c(find(contains(n,'ExperimentalBaseline_boxcar'),1))=-1; c(find(contains(n,target_label),1))=1; assert(sum(c~=0)==2);
B=double(S.nirs.beta); rss=double(S.nirs.ResSS); if ~isvector(rss),rss=diag(rss)';end; df=double(S.xX.erdf); beta=c'*B; se=sqrt(max((c'*double(S.xX.Bcov)*c).*(rss/df),0)); t=beta./se; p=betainc(df./(df+t.^2),df/2,.5); q=bh(p);
meta=readtable(fullfile(fl,'nirs_spm_inputs',char(sid),char(sid+'_channel_metadata.csv')),'TextType','string'); nch=size(B,2);
T=table(repmat(sid,nch,1),repmat(hb,nch,1),(1:nch)',string(meta.channel_id),string(meta.channel_type),string(meta.roi),string(meta.laterality),beta(:),se(:),t(:),p(:),q(:),repmat(join(string(c'),' '),nch,1),'VariableNames',{'subject_id','hemoglobin','channel_index','channel_id','channel_type','roi','laterality','contrast_beta','contrast_se','contrast_t','contrast_p_two_sided','contrast_q_bh61','contrast_vector'});
qkey="q1";condition_key="lile";if contains(target_label,"HILE"),qkey="q2";condition_key="hile";elseif contains(target_label,"LIHE"),qkey="q3";condition_key="lihe";elseif contains(target_label,"HIHE"),qkey="q4";condition_key="hihe";end
writetable(T,fullfile(out,char(sid+'_'+hb+'_'+qkey+'_'+condition_key+'_minus_experimentalbaseline_'+phase_key+'_channel_contrasts.csv')));
end

function q=bh(p)
q=nan(size(p)); v=isfinite(p); pv=p(v); [s,o]=sort(pv(:)); m=numel(s); a=min(1,flipud(cummin(flipud(s.*m./(1:m)')))); z=zeros(size(pv)); z(o)=a; q(v)=z;
end
