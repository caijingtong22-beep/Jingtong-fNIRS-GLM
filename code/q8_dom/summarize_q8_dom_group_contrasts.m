function [group_results,overview]=summarize_q8_dom_group_contrasts(phase_key)
% Unweighted n=30 one-sample summary of Q8 DOM subject-level contrasts.
if nargin<1,phase_key="dom";end;phase_key=lower(string(phase_key));repo_root=fileparts(fileparts(fileparts(mfilename('fullpath'))));root=fullfile(repo_root,'results',sprintf('q8_hihe_vs_lile_%s',phase_key),'cohort_specification');subs=1:30;group_results=table();overview=table();
for hb=["HbO","HbR"]
 B=[];info=table();
 for sn=subs
  sid=sprintf('subject%03d',sn);T=readtable(fullfile(root,sid,sprintf('%s_%s_q8_HIHE_minus_LILE_%s_channel_contrasts.csv',sid,hb,phase_key)),'TextType','string');assert(height(T)==61,'Expected 61 channels.');
  if isempty(info),info=T(:,{'channel_index','channel_id','channel_type','roi','laterality','qc_status'});else,assert(isequal(info.channel_index,T.channel_index)&&isequal(info.channel_id,T.channel_id),'Channel identity mismatch.');end;B(sn,:)=T.contrast_beta(:)'; %#ok<AGROW>
 end
 n=sum(isfinite(B),1);mu=mean(B,1,'omitnan');sd=std(B,0,1,'omitnan');se=sd./sqrt(n);ts=mu./se;ts(n<2|se==0)=NaN;p=tp(ts,n-1);q=bh(p);np=sum(B>0,1);nn=sum(B<0,1);
 G=[table(repmat(hb,61,1),n(:),mu(:),sd(:),se(:),ts(:),p(:),q(:),np(:),nn(:),'VariableNames',{'hemoglobin','n_subjects','mean_contrast_beta','sd_contrast_beta','group_se','group_t','group_p_two_sided','group_q_bh61','n_positive_beta','n_negative_beta'}),info];group_results=[group_results;G]; %#ok<AGROW>
 sig=q<.05;overview=[overview;table(hb,sum(sig),sum(sig&mu>0),sum(sig&mu<0),median(max(np,nn)./n),'VariableNames',{'hemoglobin','n_group_q_bh61_lt_005','n_group_positive_q_bh61_lt_005','n_group_negative_q_bh61_lt_005','median_directional_agreement'})]; %#ok<AGROW>
end
writetable(group_results,fullfile(root,sprintf('q8_HIHE_minus_LILE_%s_group_channel_summary_n30.csv',phase_key)));writetable(overview,fullfile(root,sprintf('q8_HIHE_minus_LILE_%s_group_overview_n30.csv',phase_key)));disp(overview);
end
function p=tp(t,df),p=betainc(df./(df+t.^2),df/2,.5);p(~isfinite(t)|df<=0)=NaN;end
function q=bh(p)
q=nan(size(p));v=isfinite(p);pv=p(v);if isempty(pv),return;end;[s,o]=sort(pv(:));m=numel(s);z=min(1,flipud(cummin(flipud(s.*m./(1:m)'))));qv=zeros(size(pv));qv(o)=z;q(v)=qv;
end
