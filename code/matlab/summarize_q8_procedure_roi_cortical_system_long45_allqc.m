function summarize_q8_procedure_roi_cortical_system_long45_allqc()
% Q8 Procedure: HIHE - LILE long-45 cortical-system ROI summary.
root = fullfile('C:\aliceproject\Jingtong-fNIRS-GLM', 'results', 'q8_hihe_vs_lile_procedure', 'cohort_specification');
roi_order = ["Medial Frontopolar PFC","Lateral Frontopolar PFC","DL PFC","VL PFC","Medial PFC","PMC","SMA","M1"];
systems = [repmat("Prefrontal",5,1); repmat("Motor",3,1)];
T = readtable(fullfile(root, 'q8_HIHE_minus_LILE_procedure_individual_ROI_long45_allQC_n30.csv'), 'TextType', 'string');
T.cortical_system = strings(height(T),1);
for i=1:numel(roi_order), T.cortical_system(T.roi==roi_order(i))=systems(i); end
T=T(T.cortical_system~="",:);
% Collapse left/right/medial subdivisions to one equal-weighted ROI mean per subject.
hbs=unique(T.hemoglobin,'stable'); subjects=unique(T.subject_id); individual=table();
for h=1:numel(hbs)
    for s=1:numel(subjects)
        for i=1:numel(roi_order)
            R=T(T.hemoglobin==hbs(h) & T.subject_id==subjects(s) & T.roi==roi_order(i),:);
            if ~isempty(R)
                individual=[individual; table(subjects(s),hbs(h),systems(i),roi_order(i),sum(R.n_long_channels),mean(R.mean_contrast_beta,'omitnan'),'VariableNames',{'subject_id','hemoglobin','cortical_system','roi','n_long_channels','mean_contrast_beta'})]; %#ok<AGROW>
            end
        end
    end
end
T=individual;
writetable(T, fullfile(root, 'q8_HIHE_minus_LILE_procedure_individual_ROI_cortical_system_long45_allQC_n30.csv'));
group=table();
for h=1:numel(hbs)
    H=T(T.hemoglobin==hbs(h),:); rows=table(); pvals=nan(numel(roi_order),1);
    for i=1:numel(roi_order)
        R=H(H.roi==roi_order(i),:); x=R.mean_contrast_beta; x=x(isfinite(x)); [ts,p]=one_sample_t(x); pvals(i)=p;
        rows=[rows; table(hbs(h),systems(i),roi_order(i),numel(x),mean(x),std(x),std(x)/sqrt(numel(x)),ts,p,sum(x>0),sum(x<0),'VariableNames',{'hemoglobin','cortical_system','roi','n_subjects','mean_contrast_beta','sd_contrast_beta','se_contrast_beta','t_statistic','p_value','n_positive_beta','n_negative_beta'})]; %#ok<AGROW>
    end
    rows.q_value_bh_fdr=bh_fdr(pvals); group=[group;rows]; %#ok<AGROW>
end
writetable(group, fullfile(root, 'q8_HIHE_minus_LILE_procedure_group_ROI_cortical_system_long45_allQC_n30.csv'));
fprintf('Q8 Procedure cortical-system ROI summary written to %s\n', root);
end
function [t,p]=one_sample_t(x)
n=numel(x); if n<2 || std(x)==0,t=NaN;p=NaN;return;end; t=mean(x)/(std(x)/sqrt(n));p=betainc(n/(n+t^2),n/2,.5);
end
function q=bh_fdr(p)
q=nan(size(p));v=isfinite(p);[sp,o]=sort(p(v));m=numel(sp);a=sp.*m./(1:m)';a=flipud(cummin(flipud(a)));ii=find(v);q(ii(o))=min(a,1);
end
