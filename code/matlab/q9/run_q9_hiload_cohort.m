function [audit,status] = run_q9_hiload_cohort(subject_numbers)
%RUN_Q9_HILOAD_COHORT Rebuild Q9 using the independent HILoad annotation.
% Model: [ExperimentalBaseline boxcar, HILoad HRF, constant].
% Contrast: HILoad - ExperimentalBaseline = [-1, +1, 0].
% This deliberately does NOT merge HILE and HIHE condition events.
% It writes new outputs and does not overwrite the earlier phase-based Q9.

if nargin < 1, subject_numbers = 1:30; end
subject_numbers = subject_numbers(:)';
assert(all(ismember(subject_numbers,1:30)), 'Subject numbers must be in 1:30.');

sid_all = compose('%03d',1:30)';
repo = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
glm = 'C:\aliceproject\GLM';
fl = fullfile(glm,'first_level');
event_file = fullfile(glm,'amended_snirf','event_tables_full_annotations', ...
    'event_table_all_subjects_full_annotations.csv');
root = fullfile(repo,'results','q9_hiload_vs_experimentalbaseline','cohort_specification');
if ~isfolder(root), mkdir(root); end

addpath(fl);
addpath(fullfile(glm,'toolboxes','spm'));
addpath(genpath(fullfile(glm,'toolboxes','NIRS_SPM_v4_r1')));
E = readtable(event_file,'TextType','string');
E.subject_id = string(compose('%03d',double(E.subject_id)));
E.condition = string(E.condition);

input_files = strings(30,1); input_exists = false(30,1);
for i = 1:30
    input_files(i) = fullfile(fl,'nirs_spm_inputs',char("subject"+sid_all(i)), ...
        char("subject"+sid_all(i)+"_nirs_data.mat"));
    input_exists(i) = isfile(input_files(i));
end

audit = table();
for i = 1:30
    sid = sid_all(i); sub = E(E.subject_id==sid,:);
    nbase = sum(sub.condition=="ExperimentalBaseline");
    nhiload = sum(sub.condition=="HILoad");
    included = nbase==1 && nhiload>=1 && input_exists(i);
    reason = "included";
    if ~included
        why = strings(0,1);
        if nbase~=1, why(end+1) = "ExperimentalBaseline count is not 1"; end
        if nhiload<1, why(end+1) = "HILoad annotation missing"; end
        if ~input_exists(i), why(end+1) = "continuous nirs_data.mat missing"; end
        reason = strjoin(why,'; ');
    end
    audit = [audit; table(sid,i,nbase,nhiload,input_exists(i),included,input_files(i),reason, ...
        'VariableNames',{'subject_id','subject_num','experimentalbaseline_n','hiload_n','input_exists','included','nirs_data_file','exclusion_reason'})]; %#ok<AGROW>
end
writetable(audit,fullfile(root,'q9_hiload_event_availability_audit_n30.csv'));
assert(sum(audit.included)==29 && any(audit.subject_id=="003" & ~audit.included), ...
    'Expected 29 included subjects with subject003 excluded for missing ExperimentalBaseline.');

status = table(); target_name = 'HILoad_HRF'; hbs = ["HbO","HbR"];
for subject_num = subject_numbers
    sid = string(compose('%03d',subject_num));
    a = audit(audit.subject_id==sid,:);
    if ~a.included
        status = [status; table(sid,"excluded",a.exclusion_reason, ...
            'VariableNames',{'subject_id','status','detail'})]; %#ok<AGROW>
        continue
    end
    sr = fullfile(root,char("subject"+sid)); if ~isfolder(sr), mkdir(sr); end
    try
        input_file = char(a.nirs_data_file);
        D = load(input_file,'nirs_data');
        nscan = size(D.nirs_data.oxyData,1); fs = double(D.nirs_data.fs);
        scan_time = (0:nscan-1)'/fs;
        sub = E(E.subject_id==sid,:);
        base = sortrows(sub(sub.condition=="ExperimentalBaseline",:),'onset_sec');
        target = sortrows(sub(sub.condition=="HILoad",:),'onset_sec');
        baseline = scan_time>=double(base.onset_sec) & scan_time<double(base.onset_sec+base.duration_sec);
        assert(any(baseline) && any(~baseline) && height(target)>=1, 'Invalid Q9 events for %s.',sid);

        selected = [base;target];
        selected.role = [repmat("ExperimentalBaseline_boxcar",height(base),1); repmat("HILoad_HRF",height(target),1)];
        selected.selection_policy = ["single ExperimentalBaseline event"; repmat("independent timeline HILoad annotation",height(target),1)];
        selected.input_file = repmat(string(input_file),height(selected),1);
        writetable(selected,fullfile(sr,char(sid+'_q9_hiload_selected_events.csv')));
        writetable(table(repmat(sid,nscan,1),(1:nscan)',scan_time,double(baseline), ...
            'VariableNames',{'subject_id','scan_index','scan_time_sec','ExperimentalBaseline_boxcar'}), ...
            fullfile(sr,char(sid+'_q9_hiload_baseline_boxcar_audit.csv')));

        names = {target_name}; onsets = {double(target.onset_sec(:)')}; durations = {double(target.duration_sec(:)')};
        save(fullfile(sr,char(sid+'_q9_hiload_task_definition.mat')),'names','onsets','durations');
        qa = table(); designs = cell(2,1);
        for ih = 1:2
            hb = hbs(ih); short = 'hbr'; if hb=="HbO", short='hbo'; end
            hr = fullfile(sr,char(hb)); if ~isfolder(hr), mkdir(hr); end
            specification_batch(input_file,short,'DCT, 128','none',0,char(hr),0,0,1,names,onsets,durations);
            sf = fullfile(hr,sprintf('SPM_indiv_%s.mat',hb)); L = load(sf,'SPM_nirs');
            S = add_baseline(L.SPM_nirs,baseline,target_name); S.nirs.step='specification'; S = rebuild_filtered_design(S);
            SPM_nirs = S; save(sf,'SPM_nirs');
            [q,row,pair] = design_qa(S,sid,hb,target_name);
            qa = [qa;q]; %#ok<AGROW>
            writematrix(full(S.xX.X),fullfile(hr,char(sid+'_Hb'+extractAfter(hb,2)+'_q9_hiload_final_design.csv')));
            writetable(row,fullfile(hr,char(sid+'_q9_hiload_regressor_table.csv')));
            writetable(pair,fullfile(hr,char(sid+'_q9_hiload_pairwise_correlation_vif.csv')));
            designs{ih} = S;
        end
        hcheck = isequal(designs{1}.xX.X(:,1),designs{2}.xX.X(:,1)) && ...
            isequal(designs{1}.xX.X(:,3),designs{2}.xX.X(:,3)) && ...
            max(abs(full(designs{1}.xX.X(:,2))+full(designs{2}.xX.X(:,2))))<1e-10;
        assert(all(qa.qa_pass) && hcheck, 'Q9 HILoad design QA failed for %s.',sid);
        writetable(qa,fullfile(sr,char(sid+'_q9_hiload_design_qa.csv')));

        for ih = 1:2
            hb = hbs(ih); sf = fullfile(sr,char(hb),sprintf('SPM_indiv_%s.mat',hb));
            nirs_spm_estimation_batch_fixed(sf,input_file);
            L = load(sf,'SPM_nirs'); write_contrast(L.SPM_nirs,sid,hb,sr,fl,target_name);
        end
        status = [status;table(sid,"completed","",'VariableNames',{'subject_id','status','detail'})]; %#ok<AGROW>
    catch ME
        status = [status;table(sid,"failed",string(ME.message),'VariableNames',{'subject_id','status','detail'})]; %#ok<AGROW>
    end
    writetable(status,fullfile(root,'q9_hiload_estimation_manifest_progress.csv'));
end
writetable(status,fullfile(root,'q9_hiload_estimation_manifest_n30.csv'));
assert(~any(status.status=="failed"),'One or more Q9 HILoad subjects failed; inspect the manifest.');
end

function S=add_baseline(S,baseline,target_name)
X=full(S.xX.X); n=string(S.xX.name(:)); ti=find(contains(n,target_name),1); ci=find(contains(lower(n),'constant'),1);
assert(size(X,2)==2 && ~isempty(ti) && ~isempty(ci));
S.xX.X=[double(baseline(:)),X(:,ti),X(:,ci)];
S.xX.name={'Sn(1) ExperimentalBaseline_boxcar (unconvolved)',char(n(ti)),'Sn(1) constant'};
S.xX.iH=[]; S.xX.iC=[1 2]; S.xX.iB=3; S.xX.iG=[]; S.Sess(1).C.C=double(baseline(:)); S.Sess(1).C.name={'ExperimentalBaseline_boxcar (unconvolved)'}; S.Sess(1).col=[1 2];
for f={'xKXs','pKX','Bcov'}, if isfield(S.xX,f{1}), S.xX=rmfield(S.xX,f{1}); end, end
end
function S=rebuild_filtered_design(S)
S.xX.xKXs=spm_sp('Set',spm_filter_HPF_LPF_WMDL(S.xX.K,S.xX.X)); S.xX.xKXs.X=full(S.xX.xKXs.X); S.xX.pKX=spm_sp('x-',S.xX.xKXs);
end
function [q,row,pair]=design_qa(S,sid,hb,target_name)
X=full(S.xX.X); n=string(S.xX.name(:)); b=X(:,1); t=X(:,2); bc=b-mean(b); tc=t-mean(t); r=(bc'*tc)/sqrt((bc'*bc)*(tc'*tc)); Z=[ones(size(b)),b]; r2=1-sum((t-Z*(Z\t)).^2)/sum((t-mean(t)).^2); vif=1/max(1-r2,eps);
q=table(sid,hb,size(X,1),size(X,2),rank(X),r,vif,rank(X)==3&&abs(r)<.95&&vif<=10,'VariableNames',{'subject_id','hemoglobin','n_scans','n_columns','rank','task_baseline_correlation','task_vif','qa_pass'});
row=table((1:3)',n,[true;false;false],[false;true;false],[false;false;true],'VariableNames',{'column_index','regressor','is_unconvolved_baseline','is_hrf_convolved_task','is_constant'});
pair=table("ExperimentalBaseline_boxcar",string(target_name),r,vif,abs(r)<.95&&vif<=10,'VariableNames',{'regressor_1','regressor_2','correlation','vif_for_task','qa_pass'});
end
function write_contrast(S,sid,hb,sr,fl,target_name)
n=string(S.xX.name(:)); ib=find(contains(n,'ExperimentalBaseline_boxcar'),1); it=find(contains(n,target_name),1); c=zeros(3,1); c(ib)=-1; c(it)=1;
B=double(S.nirs.beta); rss=double(S.nirs.ResSS); if ~isvector(rss), rss=diag(rss)'; end; df=double(S.xX.erdf); beta=c'*B; se=sqrt(max((c'*double(S.xX.Bcov)*c).*(rss/df),0)); tt=beta./se; p=betainc(df./(df+tt.^2),df/2,.5); q=bh(p);
subject_dir=char("subject"+sid); meta=readtable(fullfile(fl,'nirs_spm_inputs',subject_dir,[subject_dir '_channel_metadata.csv']),'TextType','string'); nch=size(B,2);
T=table(repmat(sid,nch,1),repmat(hb,nch,1),(1:nch)',string(meta.channel_id),string(meta.channel_type),string(meta.roi),string(meta.laterality),beta(:),se(:),tt(:),p(:),q(:),repmat(join(string(c'),' '),nch,1),'VariableNames',{'subject_id','hemoglobin','channel_index','channel_id','channel_type','roi','laterality','contrast_beta','contrast_se','contrast_t','contrast_p_two_sided','contrast_q_bh61','contrast_vector'});
writetable(T,fullfile(sr,char(sid+'_'+hb+'_q9_hiload_minus_experimentalbaseline_channel_contrasts.csv')));
end
function q=bh(p)
q=nan(size(p)); v=isfinite(p); pv=p(v); [s,o]=sort(pv(:)); m=numel(s); a=min(1,flipud(cummin(flipud(s.*m./(1:m)')))); z=zeros(size(pv)); z(o)=a; q(v)=z;
end
