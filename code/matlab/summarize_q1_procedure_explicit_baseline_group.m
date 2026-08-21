function summarize_q1_procedure_explicit_baseline_group(phase_key,question_key)
% n=29 Q1 Procedure channel and 45-long-channel ROI summaries.
if nargin<1, phase_key="procedure"; end; if nargin<2,question_key="q1";end; phase_key=lower(string(phase_key));question_key=lower(string(question_key));condition_key="lile";if question_key=="q2",condition_key="hile";elseif question_key=="q3",condition_key="lihe";elseif question_key=="q4",condition_key="hihe";end;tag=question_key+"_"+phase_key;
repo=fileparts(fileparts(fileparts(mfilename('fullpath')))); root=fullfile(repo,'results',char(question_key+"_"+condition_key+"_vs_experimentalbaseline_"+phase_key)); cohort=fullfile(root,'cohort_specification');
all=table();
for id=1:30
    if id==3, continue; end
    sid=string(sprintf('subject%03d',id)); src=fullfile(cohort,char(sid)); if id==1, src=fullfile(root,'subject001_pilot'); end
    for hb=["HbO","HbR"]
        f=fullfile(src,char(sid+'_'+hb+'_'+question_key+'_'+condition_key+'_minus_experimentalbaseline_'+phase_key+'_channel_contrasts.csv'));
        if ~isfile(f)
            legacy=fullfile(src,char(sid+'_'+hb+'_q1_lile_minus_experimentalbaseline_'+phase_key+'_channel_contrasts.csv'));
            if isfile(legacy), f=legacy; end
        end
        assert(isfile(f),'Missing contrast file: %s',f); T=readtable(f,'TextType','string'); T.subject_id=string(T.subject_id); T.hemoglobin=string(T.hemoglobin); T.channel_type=string(T.channel_type); T.roi=string(T.roi); T.laterality=string(T.laterality); all=[all;T]; %#ok<AGROW>
    end
end
assert(numel(unique(all.subject_id))==29 && height(all)==29*2*61,'Incomplete cohort contrasts.');
writetable(all,fullfile(cohort,char(tag+'_all_subject_channel_contrasts_n29.csv')));
channel=table();
for hb=["HbO","HbR"]
 for ch=1:61
    A=all(all.hemoglobin==hb & all.channel_index==ch,:); b=double(A.contrast_beta); [meanb,sdb,seb,t,p]=stats(b);
    channel=[channel;table(hb,ch,string(A.channel_id(1)),string(A.channel_type(1)),string(A.roi(1)),string(A.laterality(1)),numel(b),meanb,sdb,seb,t,p,'VariableNames',{'hemoglobin','channel_index','channel_id','channel_type','roi','laterality','n_subjects','mean_contrast_beta','sd_contrast_beta','se_contrast_beta','t_one_sample','p_two_sided'})]; %#ok<AGROW>
 end
end
channel.q_bh61=nan(height(channel),1); for hb=["HbO","HbR"], k=channel.hemoglobin==hb; channel.q_bh61(k)=bh(channel.p_two_sided(k)); end
writetable(channel,fullfile(cohort,char(question_key+'_'+condition_key+'_minus_experimentalbaseline_'+phase_key+'_group_channel_summary_n29.csv')));

long=all(all.channel_type=="Long",:); assert(numel(unique(long.channel_index))==45,'Expected 45 long channels.');
individual=groupsummary(long,{'subject_id','hemoglobin','roi','laterality'},'mean','contrast_beta');
individual.Properties.VariableNames{'mean_contrast_beta'}='roi_contrast_beta';
writetable(individual,fullfile(cohort,char(question_key+'_'+condition_key+'_minus_experimentalbaseline_'+phase_key+'_individual_ROI_long45_allQC_n29.csv')));
roi=table(); keys=unique(individual(:,{'hemoglobin','roi','laterality'}),'rows');
for i=1:height(keys)
 A=individual(individual.hemoglobin==keys.hemoglobin(i) & individual.roi==keys.roi(i) & individual.laterality==keys.laterality(i),:); b=double(A.roi_contrast_beta); [meanb,sdb,seb,t,p]=stats(b);
 roi=[roi;table(keys.hemoglobin(i),keys.roi(i),keys.laterality(i),numel(b),meanb,sdb,seb,t,p,'VariableNames',{'hemoglobin','roi','laterality','n_subjects','mean_contrast_beta','sd_contrast_beta','se_contrast_beta','t_one_sample','p_two_sided'})]; %#ok<AGROW>
end
roi.q_bh=nan(height(roi),1); for hb=["HbO","HbR"],k=roi.hemoglobin==hb;roi.q_bh(k)=bh(roi.p_two_sided(k));end
writetable(roi,fullfile(cohort,char(question_key+'_'+condition_key+'_minus_experimentalbaseline_'+phase_key+'_group_ROI_long45_allQC_n29.csv')));
end

function [m,s,se,t,p]=stats(x)
n=numel(x);m=mean(x);s=std(x,0);se=s/sqrt(n);t=m/se;df=n-1;p=betainc(df/(df+t.^2),df/2,.5);
end
function q=bh(p)
q=nan(size(p));v=isfinite(p);pv=p(v);[s,o]=sort(pv(:));m=numel(s);a=min(1,flipud(cummin(flipud(s.*m./(1:m)'))));z=zeros(size(pv));z(o)=a;q(v)=z;
end
