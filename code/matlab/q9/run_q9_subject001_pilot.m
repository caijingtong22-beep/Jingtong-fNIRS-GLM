function results = run_q9_subject001_pilot(phase_key)
% Q9 Subject001 pilot: high intrinsic load minus ExperimentalBaseline.
% The three phases are separate models.  HILE and HIHE are merged into one
% HRF-convolved task regressor; ExperimentalBaseline is an un-convolved
% scan-aligned boxcar covariate.
if nargin < 1 || isempty(phase_key)
    phases = ["procedure","dom","aoe"];
else
    phases = lower(string(phase_key));
end
assert(all(ismember(phases,["procedure","dom","aoe"])), ...
    'phase_key must be procedure, dom, or aoe.');

sid = "subject001";
repo = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
glm = 'C:\aliceproject\GLM';
fl = fullfile(glm,'first_level');
input_file = fullfile(fl,'nirs_spm_inputs',char(sid),char(sid+'_nirs_data.mat'));
event_file = fullfile(glm,'amended_snirf','event_tables_full_annotations', ...
    'event_table_all_subjects_full_annotations.csv');
assert(isfile(input_file) && isfile(event_file), 'Missing Q9 pilot input.');
addpath(fl); addpath(fullfile(glm,'toolboxes','spm')); 
addpath(genpath(fullfile(glm,'toolboxes','NIRS_SPM_v4_r1')));

E = readtable(event_file,'TextType','string');
E.subject_id = compose('%03d',double(E.subject_id));
E.condition = string(E.condition);
D = load(input_file,'nirs_data');
nscan = size(D.nirs_data.oxyData,1); fs = double(D.nirs_data.fs);
scan_time = (0:nscan-1)'/fs;
results = struct([]);

for ip = 1:numel(phases)
    phase = phases(ip);
    [target_labels, phase_label] = q9_labels(phase);
    root = fullfile(repo,'results','q9_high_intrinsic_vs_experimentalbaseline_'+phase_label, ...
        'subject001_pilot');
    if ~isfolder(root), mkdir(root); end
    sub = E(E.subject_id=="001",:);
    base = sortrows(sub(sub.condition=="ExperimentalBaseline",:),'onset_sec');
    target = sortrows(sub(ismember(sub.condition,target_labels),:),'onset_sec');
    assert(height(base)==1,'Subject001 must have exactly one ExperimentalBaseline.');
    assert(height(target)>=1,'Subject001 has no Q9 %s target events.',phase);
    assert(all(ismember(target.condition,target_labels)));

    role = [repmat("ExperimentalBaseline_boxcar",height(base),1); ...
        repmat("HILoad_phase_HRF",height(target),1)];
    selected = [base;target];
    selected.role = role;
    selected.selection_policy = ["single ExperimentalBaseline event"; ...
        repmat("all HILE/HIHE events for this phase, merged into one task regressor",height(target),1)];
    selected.input_file = repmat(string(input_file),height(selected),1);
    writetable(selected,fullfile(root,char(sid+"_q9_"+phase+"_selected_events.csv")));
    save(fullfile(root,char(sid+"_q9_"+phase+"_target_events.mat")), ...
        'target_labels','phase_label');

    baseline = scan_time >= double(base.onset_sec) & ...
        scan_time < double(base.onset_sec + base.duration_sec);
    assert(any(baseline) && any(~baseline),'Baseline boxcar has no variation.');
    audit = table(repmat(sid,nscan,1),(1:nscan)',scan_time,double(baseline), ...
        'VariableNames',{'subject_id','scan_index','scan_time_sec','ExperimentalBaseline_boxcar'});
    writetable(audit,fullfile(root,char(sid+"_q9_"+phase+"_baseline_boxcar_audit.csv")));

    names = {['HILoad_' char(phase_label) '_HRF']};
    onsets = {double(target.onset_sec(:)')};
    durations = {double(target.duration_sec(:)')};
    save(fullfile(root,char(sid+"_q9_"+phase+"_merged_task_definition.mat")), ...
        'names','onsets','durations','target_labels');

    qa = table(); S_by_hb = cell(2,1); hbs = ["HbO","HbR"];
    for ih = 1:2
        hb = hbs(ih); short = 'hbr'; if hb=="HbO", short='hbo'; end
        hb_root = fullfile(root,char(hb)); if ~isfolder(hb_root), mkdir(hb_root); end
        specification_batch(char(input_file),short,'DCT, 128','none',0,char(hb_root),0,0,1, ...
            names,onsets,durations);
        spm_file = fullfile(hb_root,sprintf('SPM_indiv_%s.mat',hb));
        L = load(spm_file,'SPM_nirs'); S = add_baseline(L.SPM_nirs,baseline,names{1});
        S.nirs.step = 'specification';
        S = rebuild_filtered_design(S); % before estimation, never append post-fit
        SPM_nirs = S; %#ok<NASGU>
        save(spm_file,'SPM_nirs');
        [qrow,rtable,qa_pair] = design_qa(S,sid,hb,root);
        qa = [qa;qrow]; %#ok<AGROW>
        writematrix(full(S.xX.X),fullfile(hb_root,char(sid+"_"+hb+"_q9_"+phase+"_final_design.csv")));
        writetable(rtable,fullfile(hb_root,char(sid+"_"+hb+"_q9_"+phase+"_regressor_table.csv")));
        writetable(qa_pair,fullfile(hb_root,char(sid+"_"+hb+"_q9_"+phase+"_pairwise_correlation_vif.csv")));
        S_by_hb{ih}=S;
    end
    hcheck = table(isequal(S_by_hb{1}.xX.X(:,1),S_by_hb{2}.xX.X(:,1)), ...
        isequal(S_by_hb{1}.xX.X(:,3),S_by_hb{2}.xX.X(:,3)), ...
        max(abs(full(S_by_hb{1}.xX.X(:,2))+full(S_by_hb{2}.xX.X(:,2))))<1e-10, ...
        'VariableNames',{'baseline_boxcar_matches','constant_matches','hbr_target_is_negative_hbo_target'});
    writetable(hcheck,fullfile(root,char(sid+"_q9_"+phase+"_hbo_hbr_design_check.csv")));
    writetable(qa,fullfile(root,char(sid+"_q9_"+phase+"_design_qa.csv")));
    assert(all(qa.qa_pass) && all(hcheck{1,:}), ...
        'Q9 %s Subject001 design QA failed; estimation was not run.',phase);

    contrasts = table();
    for ih = 1:2
        hb = hbs(ih); hb_root=fullfile(root,char(hb));
        spm_file=fullfile(hb_root,sprintf('SPM_indiv_%s.mat',hb));
        nirs_spm_estimation_batch_fixed(spm_file,input_file);
        S=load(spm_file,'SPM_nirs'); S=S.SPM_nirs;
        contrasts=[contrasts;write_contrast(S,sid,hb,root,fl,phase_label)]; %#ok<AGROW>
        plot_q9_figures(S,D.nirs_data,sid,hb,root,phase,1);
    end
    writetable(contrasts,fullfile(root,char(sid+"_q9_"+phase+"_contrast_summary.csv")));
    results(ip).phase = phase; results(ip).output_root = root; results(ip).design_qa = qa;
    results(ip).contrast_summary = contrasts;
end
end

function [labels,phase_label] = q9_labels(phase)
phase_label=phase;
if phase=="procedure", labels=["HILEProcedure","HIHEProcedure"];
elseif phase=="dom", labels=["DivisionOfMesoappendixHILE","DivisionOfMesoappendixHIHE"];
else, labels=["ApplicationOfEndoloopsHILE","ApplicationOfEndoloopsHIHE"]; end
end

function S=add_baseline(S,baseline,target_name)
X=full(S.xX.X); n=string(S.xX.name(:)); ti=find(contains(n,target_name),1); ci=find(contains(lower(n),'constant'),1);
assert(size(X,2)==2 && ~isempty(ti) && ~isempty(ci));
S.xX.X=[double(baseline(:)),X(:,ti),X(:,ci)];
S.xX.name={'Sn(1) ExperimentalBaseline_boxcar (unconvolved)',char(n(ti)),'Sn(1) constant'};
S.xX.iH=[]; S.xX.iC=[1 2]; S.xX.iB=3; S.xX.iG=[];
S.Sess(1).C.C=double(baseline(:)); S.Sess(1).C.name={'ExperimentalBaseline_boxcar (unconvolved)'}; S.Sess(1).col=[1 2];
for f={'xKXs','pKX','Bcov'}, if isfield(S.xX,f{1}), S.xX=rmfield(S.xX,f{1}); end, end
end

function S=rebuild_filtered_design(S)
% Recreate the filtered design and pseudoinverse while still at specification.
S.xX.xKXs=spm_sp('Set',spm_filter_HPF_LPF_WMDL(S.xX.K,S.xX.X));
S.xX.xKXs.X=full(S.xX.xKXs.X); S.xX.pKX=spm_sp('x-',S.xX.xKXs);
end

function [q,row,pair]=design_qa(S,sid,hb,root)
X=full(S.xX.X); n=string(S.xX.name(:)); assert(isequal(size(X,2),3)&&rank(X)==3);
b=X(:,1); t=X(:,2); bc=b-mean(b); tc=t-mean(t); r=(bc'*tc)/sqrt((bc'*bc)*(tc'*tc)); Z=[ones(size(b)),b]; r2=1-sum((t-Z*(Z\t)).^2)/sum((t-mean(t)).^2); vif=1/max(1-r2,eps);
pass=rank(X)==3 && abs(r)<.95 && vif<=10;
q=table(sid,hb,size(X,1),size(X,2),rank(X),r,vif,pass,'VariableNames',{'subject_id','hemoglobin','n_scans','n_columns','rank','task_baseline_correlation','task_vif','qa_pass'});
row=table((1:3)',n,[true;false;false],[false;true;false],[false;false;true],... 
 'VariableNames',{'column_index','regressor','is_unconvolved_baseline','is_hrf_convolved_task','is_constant'});
pair=table("ExperimentalBaseline_boxcar","HILoad_phase_HRF",r,vif,abs(r)<.95&&vif<=10,'VariableNames',{'regressor_1','regressor_2','correlation','vif_for_task','qa_pass'});
end

function summary=write_contrast(S,sid,hb,root,fl,phase_label)
n=string(S.xX.name(:)); ib=find(contains(n,'ExperimentalBaseline_boxcar'),1); it=find(contains(n,'HILoad_'+phase_label+'_HRF'),1); c=zeros(3,1); c(ib)=-1; c(it)=1;
B=double(S.nirs.beta); rss=double(S.nirs.ResSS); if ~isvector(rss),rss=diag(rss)';end; df=double(S.xX.erdf); beta=c'*B; se=sqrt(max((c'*double(S.xX.Bcov)*c).*(rss/df),0)); t=beta./se; p=betainc(df./(df+t.^2),df/2,.5); q=bh(p);
meta=readtable(fullfile(fl,'nirs_spm_inputs',char(sid),char(sid+'_channel_metadata.csv')),'TextType','string'); nch=size(B,2);
T=table(repmat(sid,nch,1),repmat(hb,nch,1),(1:nch)',string(meta.channel_id),string(meta.channel_type),string(meta.roi),string(meta.laterality),beta(:),se(:),t(:),p(:),q(:),repmat(join(string(c'),' '),nch,1),'VariableNames',{'subject_id','hemoglobin','channel_index','channel_id','channel_type','roi','laterality','contrast_beta','contrast_se','contrast_t','contrast_p_two_sided','contrast_q_bh61','contrast_vector'});
writetable(T,fullfile(root,char(sid+'_'+hb+'_q9_high_intrinsic_minus_experimentalbaseline_'+phase_label+'_channel_contrasts.csv')));
summary=table(hb,nch,df,mean(beta),median(beta),sum(p<.05),sum(q<.05),'VariableNames',{'hemoglobin','n_channels','effective_df','mean_contrast_beta','median_contrast_beta','n_p_lt_005','n_q_bh61_lt_005'});
end

function plot_q9_figures(S,d,sid,hb,root,phase,ch)
out=fullfile(root,char(hb),'figures'); if ~isfolder(out),mkdir(out);end; set(0,'DefaultFigureVisible','off');
X=full(S.xX.X); n=string(S.xX.name(:)); t=(0:size(X,1)-1)'*double(S.xY.RT); h=(0:size(S.xBF.bf,1)-1)'*double(S.xBF.dt);
f=figure('Color','w');plot(h,S.xBF.bf(:,1),'LineWidth',1.5);grid on;xlabel('Time after onset (s)');ylabel('HRF amplitude');title(sprintf('Q9 %s canonical HRF (%s)',upper(phase),hb));exportgraphics(f,fullfile(out,char(sid+'_'+hb+'_q9_'+phase+'_canonical_hrf.png')),'Resolution',220);close(f);
f=figure('Color','w','Position',[100 100 1300 500]);imagesc(X');axis tight;colorbar;yticks(1:3);yticklabels(n);xlabel('Scan');ylabel('Final regressor');title(sprintf('%s Q9 %s final design: baseline is unconvolved boxcar (%s)',sid,upper(phase),hb));exportgraphics(f,fullfile(out,char(sid+'_'+hb+'_q9_'+phase+'_final_design_heatmap.png')),'Resolution',220);close(f);
f=figure('Color','w','Position',[100 100 1400 500]);plot(t,X,'LineWidth',1);grid on;xlabel('Time (s)');ylabel('Regressor amplitude');legend(n,'Location','eastoutside');title(sprintf('%s Q9 %s design (%s)',sid,upper(phase),hb));exportgraphics(f,fullfile(out,char(sid+'_'+hb+'_q9_'+phase+'_final_design_timeseries.png')),'Resolution',220);close(f);
field='dxyData';if hb=="HbO",field='oxyData';end;y=double(d.(field)(:,ch));yf=spm_filter_HPF_LPF_WMDL(S.xX.K,y);fit=full(S.xX.xKXs.X)*double(S.nirs.beta(:,ch));res=yf-fit;
f=figure('Color','w','Position',[100 100 1500 900]);tiledlayout(3,1,'TileSpacing','compact');nexttile;plot(t,yf);hold on;plot(t,fit,'LineWidth',1.2);grid on;legend('Observed filtered signal','GLM fit','Location','eastoutside');title(sprintf('%s channel %d %s: observed versus GLM fit',sid,ch,hb));nexttile;plot(t,res);yline(0,'--k');grid on;title('Residual');nexttile;plot(t,X,'LineWidth',1);grid on;legend(n,'Location','eastoutside');title('Final regressors: baseline = unconvolved boxcar; task = HRF-convolved');xlabel('Time (s)');exportgraphics(f,fullfile(out,char(sid+'_channel'+sprintf('%02d',ch)+'_'+hb+'_q9_'+phase+'_fit.png')),'Resolution',220);close(f);
end

function q=bh(p)
q=nan(size(p));v=isfinite(p);pv=p(v);[s,o]=sort(pv(:));m=numel(s);a=min(1,flipud(cummin(flipud(s.*m./(1:m)'))));z=zeros(size(pv));z(o)=a;q(v)=z;
end
