function [summary,errors]=run_q5_aoe_cohort_specification(subject_numbers)
% Q5 AOE four-cell factorial design specification.
if nargin<1,subject_numbers=1;end
subject_numbers=subject_numbers(:)';assert(all(ismember(subject_numbers,1:30)));
repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));glm='C:\aliceproject\GLM';fl=fullfile(glm,'first_level');v2=fullfile(fl,'nirs_spm_inputs_v2_fulltimeline_trial1');raw=fullfile(fl,'nirs_spm_inputs');root=fullfile(repo,'results','q5_intrinsic_aoe','cohort_specification');if ~isfolder(root),mkdir(root);end
addpath(fl);addpath(fullfile(glm,'toolboxes','spm'));addpath(genpath(fullfile(glm,'toolboxes','NIRS_SPM_v4_r1')));
E=readtable(fullfile(glm,'amended_snirf','event_tables_full_annotations','event_table_all_subjects_full_annotations.csv'),'TextType','string');E.condition=string(E.condition);E.sn=double(E.subject_id);summary=table();errors=table();
src=["ApplicationOfEndoloopsLILE","ApplicationOfEndoloopsHILE","ApplicationOfEndoloopsLIHE","ApplicationOfEndoloopsHIHE"];labs={'LILEAOE','HILEAOE','LIHEAOE','HIHEAOE'};
for sn=subject_numbers
 sid=sprintf('subject%03d',sn);sr=fullfile(root,sid);if ~isfolder(sr),mkdir(sr);end
 try
  if ismember(sn,[1 2 10]),d=fullfile(v2,sid);nf=fullfile(d,sprintf('%s_nirs_data.mat',sid));ev=readtable(fullfile(d,sprintf('%s_events_used_v2.csv',sid)),'TextType','string');else,d=fullfile(raw,sid);nf=fullfile(d,sprintf('%s_nirs_data.mat',sid));ev=E(E.sn==sn,:);end
  ev.condition=string(ev.condition);S=sortrows(ev(ismember(ev.condition,src),:),'onset_sec');assert(all(ismember(src,unique(S.condition))));on=cell(1,4);du=cell(1,4);for i=1:4,ix=S.condition==src(i);assert(sum(ix)==1);on{i}=double(S.onset_sec(ix))';du{i}=double(S.duration_sec(ix))';end
  writetable(S,fullfile(sr,sprintf('%s_q5_aoe_events_used.csv',sid)));save(fullfile(sr,sprintf('%s_q5_aoe_multiple_conditions.mat',sid)),'labs','on','du');
  for h=["hbo","hbr"]
   if h=="hbo",hb="HbO";else,hb="HbR";end;hr=fullfile(sr,char(hb));if ~isfolder(hr),mkdir(hr);end
   A=specification_batch(nf,char(h),'DCT, 128','none',0,hr,0,0,1,labs,on,du);X=full(A.xX.X);rn=string(A.xX.name(:));ic=contains(lower(rn),'constant');assert(size(X,2)==5&&sum(ic)==1);Z=(X(:,~ic)-mean(X(:,~ic)))./std(X(:,~ic),0,1);R=corrcoef(Z);mr=max(abs(R(~eye(4))));v=zeros(4,1);for j=1:4,y=X(:,j);z=[ones(size(X,1),1),X(:,setdiff(1:4,j))];v(j)=1/max(1-(1-sum((y-z*(z\y)).^2)/sum((y-mean(y)).^2)),eps);end;rk=rank(X);pass=rk==5&&mr<.95&&max(v)<=10
   writematrix(X,fullfile(hr,sprintf('%s_%s_q5_aoe_convolved_design.csv',sid,hb)));writetable(table(rn(~ic),v,'VariableNames',{'regressor','vif'}),fullfile(hr,sprintf('%s_%s_q5_aoe_vif.csv',sid,hb)));summary=[summary;table(string(sid),hb,rk,mr,max(v),pass,'VariableNames',{'subject_id','hemoglobin','rank','max_abs_correlation','max_vif','qa_pass'})];fprintf('%s %s rank %d/5 r %.3f VIF %.3f\n',sid,hb,rk,mr,max(v));
  end
 catch ME,errors=[errors;table(string(sid),string(ME.message),'VariableNames',{'subject_id','error_message'})];end
end
tag=sprintf('%03d_to_%03d',min(subject_numbers),max(subject_numbers));writetable(summary,fullfile(root,sprintf('q5_aoe_design_qa_%s.csv',tag)));writetable(errors,fullfile(root,sprintf('q5_aoe_specification_errors_%s.csv',tag)));assert(isempty(errors));assert(all(summary.qa_pass));
end
