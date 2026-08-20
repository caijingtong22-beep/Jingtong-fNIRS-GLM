function summarize_q8_dom_roi_long45_allqc(phase_key)
% Q8 DOM: long-45 broad ROI/laterality summaries, retaining all QC flags.
if nargin<1,phase_key="dom";end;phase_key=lower(string(phase_key));repo_root=fileparts(fileparts(fileparts(mfilename('fullpath'))));root=fullfile(repo_root,'results',sprintf('q8_hihe_vs_lile_%s',phase_key),'cohort_specification');subs=1:30;hbs=["HbO","HbR"];I=table();C=table();
for sn=subs
 sid=sprintf('subject%03d',sn);
 for hb=hbs
  T=readtable(fullfile(root,sid,sprintf('%s_%s_q8_HIHE_minus_LILE_%s_channel_contrasts.csv',sid,hb,phase_key)),'TextType','string');T=T(T.channel_type=="Long",:);assert(height(T)==45,'%s %s lacks 45 long channels.',sid,hb);C=[C;T]; %#ok<AGROW>
  K=unique(T(:,{'roi','laterality'}),'rows','stable');
  for k=1:height(K),ix=T.roi==K.roi(k)&T.laterality==K.laterality(k);I=[I;table(sn,hb,K.roi(k),K.laterality(k),sum(ix),mean(T.contrast_beta(ix),'omitnan'),median(T.contrast_beta(ix),'omitnan'),sum(T.contrast_beta(ix)>0),sum(T.contrast_beta(ix)<0),sum(T.contrast_p_two_sided(ix)<.05),'VariableNames',{'subject_id','hemoglobin','roi','laterality','n_long_channels','mean_contrast_beta','median_contrast_beta','n_positive_beta','n_negative_beta','n_channel_p_lt_005'})]; %#ok<AGROW>
  end
 end
end
writetable(I,fullfile(root,sprintf('q8_HIHE_minus_LILE_%s_individual_ROI_long45_allQC_n30.csv',phase_key)));CG=table();RG=table();
for hb=hbs
 H=C(C.hemoglobin==hb,:);ids=unique(H.channel_id);R=table();pvals=nan(numel(ids),1);
 for j=1:numel(ids),A=H(H.channel_id==ids(j),:);x=A.contrast_beta;x=x(isfinite(x));[t,p]=tt(x);pvals(j)=p;R=[R;table(hb,ids(j),A.channel_type(1),A.roi(1),A.laterality(1),numel(x),mean(x),std(x),std(x)/sqrt(numel(x)),t,p,sum(x>0),sum(x<0),'VariableNames',{'hemoglobin','channel','channel_type','roi','laterality','n_subjects','mean_contrast_beta','sd_contrast_beta','se_contrast_beta','t_statistic','p_value','n_positive_beta','n_negative_beta'})]; %#ok<AGROW>
 end;R.q_value_bh_fdr=bh(pvals);CG=[CG;R]; %#ok<AGROW>
 H=I(I.hemoglobin==hb,:);K=unique(H(:,{'roi','laterality'}),'rows','stable');R=table();pvals=nan(height(K),1);
 for k=1:height(K),A=H(H.roi==K.roi(k)&H.laterality==K.laterality(k),:);x=A.mean_contrast_beta;x=x(isfinite(x));[t,p]=tt(x);pvals(k)=p;R=[R;table(hb,K.roi(k),K.laterality(k),numel(x),mean(x),std(x),std(x)/sqrt(numel(x)),t,p,sum(x>0),sum(x<0),'VariableNames',{'hemoglobin','roi','laterality','n_subjects','mean_contrast_beta','sd_contrast_beta','se_contrast_beta','t_statistic','p_value','n_positive_beta','n_negative_beta'})]; %#ok<AGROW>
 end;R.q_value_bh_fdr=bh(pvals);RG=[RG;R]; %#ok<AGROW>
end
writetable(CG,fullfile(root,sprintf('q8_HIHE_minus_LILE_%s_group_channel_long45_allQC_n30.csv',phase_key)));writetable(RG,fullfile(root,sprintf('q8_HIHE_minus_LILE_%s_group_ROI_long45_allQC_n30.csv',phase_key)));fprintf('Wrote Q8 %s long-45 ROI summaries to %s\n',upper(phase_key),root);
end
function [t,p]=tt(x),n=numel(x);if n<2||std(x)==0,t=NaN;p=NaN;return;end;t=mean(x)/(std(x)/sqrt(n));p=betainc(n/(n+t^2),n/2,.5);end
function q=bh(p),q=nan(size(p));v=isfinite(p);pv=p(v);[s,o]=sort(pv(:));m=numel(s);z=min(1,flipud(cummin(flipud(s.*m./(1:m)'))));qv=zeros(size(pv));qv(o)=z;q(v)=qv;end
