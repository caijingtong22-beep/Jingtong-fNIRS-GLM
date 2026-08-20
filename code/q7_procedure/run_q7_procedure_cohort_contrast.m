function summary=run_q7_procedure_cohort_contrast(subject_numbers)
% Q7 Procedure interaction contrast using validated Q5 Procedure GLM.
if nargin<1,subject_numbers=2:30;end;subject_numbers=subject_numbers(:)';assert(all(ismember(subject_numbers,2:30)));
repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));glm='C:\aliceproject\GLM';fl=fullfile(glm,'first_level');srcroot=fullfile(repo,'results','q5_intrinsic_procedure','cohort_specification');root=fullfile(repo,'results','q7_interaction_procedure','cohort_specification');addpath(fl);addpath(fullfile(glm,'toolboxes','spm'));addpath(genpath(fullfile(glm,'toolboxes','NIRS_SPM_v4_r1')));summary=table();
for sn=subject_numbers
 sid=sprintf('subject%03d',sn);srcsid=fullfile(srcroot,sid);outsid=fullfile(root,sid);if ~isfolder(outsid),mkdir(outsid);end
 for hb=["HbO","HbR"]
  sf=fullfile(srcsid,char(hb),sprintf('SPM_indiv_%s.mat',hb));assert(isfile(sf),'Missing validated Q5 Procedure SPM: %s',sf);L=load(sf,'SPM_nirs');S=L.SPM_nirs;N=string(S.xX.name(:));B=double(S.nirs.beta);nch=size(B,2);targets=["LILEProcedure","HILEProcedure","LIHEProcedure","HIHEProcedure"];k=zeros(4,1);for i=1:4,k(i)=find(contains(N,targets(i)),1);end;c=zeros(size(B,1),1);c(k)=[1 -1 -1 1];rss=double(S.nirs.ResSS);if ~isvector(rss),rss=diag(rss)';end;df=double(S.xX.erdf);b=c'*B;se=sqrt(max((c'*double(S.xX.Bcov)*c).*(rss/df),0));t=b./se;p=betainc(df./(df+t.^2),df/2,.5);q=bh(p);meta=readtable(fullfile(fl,'nirs_spm_inputs',sid,sprintf('%s_channel_metadata.csv',sid)),'TextType','string');R=table(repmat(string(sid),nch,1),repmat(hb,nch,1),(1:nch)',string(meta.channel_id),string(meta.channel_type),string(meta.roi),string(meta.laterality),b(:),se(:),t(:),p(:),q(:),repmat(join(string(c'),' '),nch,1),'VariableNames',{'subject_id','hemoglobin','channel_index','channel_id','channel_type','roi','laterality','contrast_beta','contrast_se','contrast_t','contrast_p_two_sided','contrast_q_bh61','contrast_vector'});writetable(R,fullfile(outsid,sprintf('%s_%s_q7_interaction_procedure_channel_contrasts.csv',sid,hb)));summary=[summary;table(string(sid),hb,df,sum(p<.05),sum(q<.05),'VariableNames',{'subject_id','hemoglobin','effective_df','n_p_lt_005','n_q_bh61_lt_005'})];
 end
end
writetable(summary,fullfile(root,'q7_procedure_estimation_contrast_n29_reusing_q5_glm.csv'));disp(summary);
end
function q=bh(p),q=nan(size(p));v=isfinite(p);pv=p(v);[s,o]=sort(pv(:));m=numel(s);z=min(1,flipud(cummin(flipud(s.*m./(1:m)'))));zv=zeros(size(pv));zv(o)=z;q(v)=zv;end
