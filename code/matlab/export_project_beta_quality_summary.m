function export_project_beta_quality_summary()
%EXPORT_PROJECT_BETA_QUALITY_SUMMARY Create a compact Q1--Q10 group-beta audit.
% Existing group-level CSVs are read only; no GLM is re-estimated.

repo=fileparts(fileparts(fileparts(mfilename('fullpath')))); results=fullfile(repo,'results');
out=fullfile(results,'group_outputs','project_beta_quality_summary'); if ~isfolder(out),mkdir(out);end
rows=table('Size',[0 12], ...
 'VariableTypes',{'string','string','string','double','double','double','double','double','double','double','double','double'}, ...
 'VariableNames',{'question','stage','hemoglobin','n_subjects','n_long_channels','n_nominal_channel_p_lt_005','n_fdr_channel_q_lt_005','min_p','mean_beta','median_abs_beta','max_abs_beta','n_nonfinite_beta'});
questions=["Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10"];
stages=["procedure","dom","aoe"];
for q=1:9
    for s=1:3
        f=group_file(results,q,s);
        d=readtable(f,'TextType','string');
        if ismember('channel_type',d.Properties.VariableNames), d=d(strcmpi(d.channel_type,'Long'),:); end
        hb=string(d.hemoglobin); for h=unique(hb,'stable')'
            x=d(hb==h,:); p=getvar(x,{'p_two_sided','group_p_two_sided','p_value'}); qv=getvar(x,{'q_bh','q_bh61','group_q_bh45','group_q_bh61','q_value_bh_fdr'}); b=x.mean_contrast_beta;
            rows(end+1,:)={questions(q),upper(stages(s)),h,x.n_subjects(1),height(x),sum(p<.05),sum(qv<.05),min(p),mean(b,'omitnan'),median(abs(b),'omitnan'),max(abs(b)),sum(~isfinite(b))}; %#ok<AGROW>
        end
    end
end
f=fullfile(results,'q10','q10_HE_minus_LE_group_channel_long45_allQC_n30.csv'); d=readtable(f,'TextType','string'); hb=string(d.hemoglobin);
for h=unique(hb,'stable')'
    x=d(hb==h,:); p=x.group_p_two_sided; qv=x.group_q_bh45; b=x.mean_contrast_beta;
    rows(end+1,:)={"Q10","Bleeding",h,x.n_subjects(1),height(x),sum(p<.05),sum(qv<.05),min(p),mean(b,'omitnan'),median(abs(b),'omitnan'),max(abs(b)),sum(~isfinite(b))}; %#ok<AGROW>
end
writetable(rows,fullfile(out,'q1_q10_group_beta_quality_summary.csv'));
fprintf('Wrote %d rows to %s\n',height(rows),out);
end

function f=group_file(results,q,s)
stage=["procedure","dom","aoe"]; names=["q1_lile_vs_experimentalbaseline_","q2_hile_vs_experimentalbaseline_","q3_lihe_vs_experimentalbaseline_","q4_hihe_vs_experimentalbaseline_","q5_intrinsic_","q6_extraneous_","q7_interaction_","q8_hihe_vs_lile_","q9_high_intrinsic_vs_experimentalbaseline_"];
stems=["q1_lile_minus_experimentalbaseline_","q2_hile_minus_experimentalbaseline_","q3_lihe_minus_experimentalbaseline_","q4_hihe_minus_experimentalbaseline_","q5_high_intrinsic_minus_low_intrinsic_","q6_high_extraneous_minus_low_extraneous_","q7_interaction_","q8_HIHE_minus_LILE_","q9_high_intrinsic_minus_experimentalbaseline_"];
base=fullfile(results,names(q)+stage(s),'cohort_specification');
if q==9, base=fullfile(base,'group_analysis_n29'); end
if q==8
    f=fullfile(base,stems(q)+stage(s)+"_group_channel_long45_allQC_n30.csv");
elseif q==9
    f=fullfile(base,stems(q)+stage(s)+"_group_channel_summary_61ch_n29.csv");
else
    n=30;
    if q<=4 || q==9, n=29; end
    f=fullfile(base,stems(q)+stage(s)+"_group_channel_summary_n"+string(n)+".csv");
end
end

function v=getvar(t,names)
for k=1:numel(names)
    n=names{k};
    if ismember(n,t.Properties.VariableNames), v=t.(n); return; end
end
error('No requested variable found.');
end
