function [summary,errors]=run_q8_dom_cohort_compute_contrast(subject_numbers,phase_key)
% First-level Q8 DOM contrast: HIHEDOM minus LILEDOM.
if nargin<1,subject_numbers=1;end;if nargin<2,phase_key="dom";end;phase_key=lower(string(phase_key));assert(ismember(phase_key,["dom","aoe"]));if phase_key=="dom",low="LILEDOM";high="HIHEDOM";else,low="LILEAOE";high="HIHEAOE";end;subject_numbers=unique(subject_numbers(:)');assert(all(ismember(subject_numbers,1:30)));
repo_root=fileparts(fileparts(fileparts(mfilename('fullpath'))));glm_root='C:\aliceproject\GLM';root=fullfile(repo_root,'results',sprintf('q8_hihe_vs_lile_%s',phase_key),'cohort_specification');meta_root=fullfile(glm_root,'first_level','nirs_spm_inputs');summary=table();errors=table();
for sn=subject_numbers
 sid=sprintf('subject%03d',sn);
 try
  meta=readtable(fullfile(meta_root,sid,sprintf('%s_channel_metadata.csv',sid)),'TextType','string');
  for hb=["HbO","HbR"]
   L=load(fullfile(root,sid,char(hb),sprintf('SPM_indiv_%s.mat',hb)),'SPM_nirs');S=L.SPM_nirs;assert(strcmpi(string(S.nirs.step),'estimation'),'Not estimated.');beta=double(S.nirs.beta);names=string(S.xX.name(:));il=find(contains(names,low),1);ih=find(contains(names,high),1);assert(~isempty(il)&&~isempty(ih),'Q8 task regressors are missing.');c=zeros(size(beta,1),1);c(il)=-1;c(ih)=1;nch=size(beta,2);assert(height(meta)==nch,'Metadata channel count mismatch.');rss=resss(S.nirs.ResSS,nch);df=double(S.xX.erdf);bcov=double(S.xX.Bcov);cb=c'*beta;cv=max((c'*bcov*c).*(rss./df),0);se=sqrt(cv);ts=cb./se;p=tp(ts,df);q=bh(p);
   R=table(repmat(string(sid),nch,1),repmat(hb,nch,1),(1:nch)',string(meta.channel_id),string(meta.channel_type),string(meta.roi),string(meta.laterality),string(meta.qc_status),cb(:),se(:),ts(:),p(:),q(:),cv(:),repmat(df,nch,1),repmat(join(string(c'),' '),nch,1),'VariableNames',{'subject_id','hemoglobin','channel_index','channel_id','channel_type','roi','laterality','qc_status','contrast_beta','contrast_se','contrast_t','contrast_p_two_sided','contrast_q_bh61','contrast_variance','effective_df','contrast_vector'});
   writetable(R,fullfile(root,sid,sprintf('%s_%s_q8_HIHE_minus_LILE_%s_channel_contrasts.csv',sid,hb,phase_key)));summary=[summary;table(string(sid),hb,nch,sum(isfinite(cb)),sum(p<.05),sum(q<.05),min(p),median(p),'VariableNames',{'subject_id','hemoglobin','n_channels','n_finite_beta','n_p_lt_005','n_q_bh61_lt_005','min_p_two_sided','median_p_two_sided'})]; %#ok<AGROW>
   fprintf('%s %s: %d p<.05, %d BH-FDR<.05\n',sid,hb,sum(p<.05),sum(q<.05));
  end
 catch ME
  errors=[errors;table(string(sid),string(ME.identifier),string(getReport(ME,'extended','hyperlinks','off')),'VariableNames',{'subject_id','error_id','error_message'})]; %#ok<AGROW>
 end
end
tag=sprintf('%03d_to_%03d',min(subject_numbers),max(subject_numbers));writetable(summary,fullfile(root,sprintf('q8_%s_contrast_%s.csv',phase_key,tag)));writetable(errors,fullfile(root,sprintf('q8_%s_contrast_errors_%s.csv',phase_key,tag)));writetable(table("Q8_HIHE_minus_LILE_"+upper(phase_key),high+" - "+low,"[-1, +1, 0], located by regressor name",'VariableNames',{'contrast_name','description','implemented_vector'}),fullfile(root,sprintf('q8_%s_contrast_definition.csv',phase_key)));assert(isempty(errors),'Q8 %s contrast had %d failures.',upper(phase_key),height(errors));
end
function r=resss(v,n)
v=double(v);if isvector(v),r=v(:)';elseif isequal(size(v),[n n]),r=diag(v)';else,error('Unexpected ResSS dimensions.');end
end
function p=tp(t,df),p=betainc(df./(df+t.^2),df/2,.5);p(~isfinite(t))=NaN;end
function q=bh(p)
q=nan(size(p));v=isfinite(p);pv=p(v);if isempty(pv),return;end;[s,o]=sort(pv(:));m=numel(s);z=s.*m./(1:m)';z=min(1,flipud(cummin(flipud(z))));qv=zeros(size(pv));qv(o)=z;q(v)=qv;
end
