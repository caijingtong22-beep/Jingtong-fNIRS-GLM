function out = run_q9_procedure_group_analysis(phase_key)
% Q9 Procedure group summaries and figures, explicit-baseline n=29.
% Includes Subject001 pilot plus 28 eligible cohort subjects; Subject003 is excluded.
if nargin<1, phase_key="procedure"; end; phase_key=lower(string(phase_key)); assert(any(phase_key==["procedure","dom","aoe"]));
project = 'C:\aliceproject\Jingtong-fNIRS-GLM';
cohort = fullfile(project,'results',char("q9_high_intrinsic_vs_experimentalbaseline_"+phase_key),'cohort_specification');
pilot = fullfile(project,'results',char("q9_high_intrinsic_vs_experimentalbaseline_"+phase_key),'subject001_pilot');
figroot = fullfile(cohort,'group_analysis_n29'); if ~isfolder(figroot),mkdir(figroot);end
prefix = "q9_high_intrinsic_minus_experimentalbaseline_"+phase_key; figprefix="q9_"+phase_key;
hbs=["HbO","HbR"]; all_rows=table(); subject_ids=[];
for hb=hbs
    files=dir(fullfile(cohort,'subject*')); files=files([files.isdir]); picked=dir(fullfile(pilot,sprintf('*_%s_q9_high_intrinsic_minus_experimentalbaseline_%s_channel_contrasts.csv',char(hb),char(phase_key))));
    allfiles=picked;
    for di=1:numel(files)
        ff=dir(fullfile(cohort,files(di).name,sprintf('*_%s_q9_high_intrinsic_minus_experimentalbaseline_%s_channel_contrasts.csv',char(hb),char(phase_key))));
        allfiles=[allfiles;ff]; %#ok<AGROW>
    end
    files=allfiles;
    assert(numel(files)==29,'Expected 29 %s contrast files, found %d.',hb,numel(files));
    for k=1:numel(files)
        T=readtable(fullfile(files(k).folder,files(k).name),'TextType','string'); assert(height(T)==61,'Expected 61 channels: %s',files(k).name);
        all_rows=[all_rows;T]; %#ok<AGROW>
        subject_ids=unique([subject_ids;string(T.subject_id)]); %#ok<AGROW>
    end
end
subject_ids=unique(subject_ids,'stable'); n= numel(subject_ids); assert(n==29);
all_rows.subject_id=regexprep(string(all_rows.subject_id),'^subject',''); subject_ids=regexprep(string(subject_ids),'^subject',''); subject_ids=unique(subject_ids,'stable'); subject_ids=subject_ids(:); n=numel(subject_ids); assert(n==29);
all_rows.hemoglobin=string(all_rows.hemoglobin); all_rows.channel_id=string(all_rows.channel_id); all_rows.channel_index=double(all_rows.channel_index); all_rows.channel_type=string(all_rows.channel_type); all_rows.roi=string(all_rows.roi); all_rows.laterality=string(all_rows.laterality);

% All 61 channels.
channel_group=table();
for hb=hbs
    H=all_rows(all_rows.hemoglobin==hb,:); [~,first_idx]=unique(H.channel_index,'first'); first_idx=sort(first_idx); keys=H(first_idx,{'channel_id','channel_index','channel_type','roi','laterality'}); assert(height(keys)==61); rows=table(); pvals=nan(height(keys),1);
    for k=1:height(keys)
        x=H.contrast_beta(H.channel_index==keys.channel_index(k)); x=x(isfinite(x)); [tt,p]=one_sample_t(x); pvals(k)=p;
        rows=[rows;table(hb,keys.channel_index(k),keys.channel_id(k),keys.channel_type(k),keys.roi(k),keys.laterality(k),numel(x),mean(x),std(x),std(x)/sqrt(numel(x)),tt,p,sum(x>0),sum(x<0),'VariableNames',{'hemoglobin','channel_index','channel_id','channel_type','roi','laterality','n_subjects','mean_contrast_beta','sd_contrast_beta','se_contrast_beta','t_statistic','p_value','n_positive_beta','n_negative_beta'})]; %#ok<AGROW>
    end
    rows.q_value_bh_fdr=bh_fdr(pvals); channel_group=[channel_group;rows]; %#ok<AGROW>
end
writetable(channel_group,fullfile(figroot,char(prefix+"_group_channel_summary_61ch_n29.csv")));

% Individual ROI means and group ROI tests; all 45 physical long channels retained.
ind=table(); roi_group=table(); keys=unique(all_rows(all_rows.hemoglobin=="HbO" & all_rows.channel_type=="Long",{'roi','laterality'}),'rows');
for is=1:numel(subject_ids)
  sid=subject_ids(is);
  for hb=hbs
    sid=string(sid); hb=string(hb);
    T=all_rows(all_rows.subject_id==sid & all_rows.hemoglobin==hb & all_rows.channel_type=="Long",:); assert(height(T)==45);
    for k=1:height(keys)
      x=T.contrast_beta(T.roi==keys.roi(k)&T.laterality==keys.laterality(k)); x=x(isfinite(x));
      ind=[ind;table(sid,hb,keys.roi(k),keys.laterality(k),numel(x),mean(x),median(x),sum(x>0),sum(x<0),'VariableNames',{'subject_id','hemoglobin','roi','laterality','n_long_channels','mean_contrast_beta','median_contrast_beta','n_positive_beta','n_negative_beta'})]; %#ok<AGROW>
    end
  end
end
for hb=hbs
  H=ind(ind.hemoglobin==hb,:); pvals=nan(height(keys),1); rows=table();
  for k=1:height(keys)
    x=H.mean_contrast_beta(H.roi==keys.roi(k)&H.laterality==keys.laterality(k));x=x(isfinite(x));[tt,p]=one_sample_t(x);pvals(k)=p;
    rows=[rows;table(hb,keys.roi(k),keys.laterality(k),numel(x),mean(x),std(x),std(x)/sqrt(numel(x)),tt,p,sum(x>0),sum(x<0),'VariableNames',{'hemoglobin','roi','laterality','n_subjects','mean_contrast_beta','sd_contrast_beta','se_contrast_beta','t_statistic','p_value','n_positive_beta','n_negative_beta'})]; %#ok<AGROW>
  end
  rows.q_value_bh_fdr=bh_fdr(pvals);roi_group=[roi_group;rows]; %#ok<AGROW>
end
writetable(ind,fullfile(figroot,char(prefix+"_individual_ROI_long45_allQC_n29.csv")));
writetable(roi_group,fullfile(figroot,char(prefix+"_group_ROI_long45_allQC_n29.csv")));

% Cortical-system aggregation used by Q5-Q8: eight established systems/ROIs.
systems=["Medial Frontopolar PFC","Lateral Frontopolar PFC","DL PFC","VL PFC","Medial PFC","PMC","SMA","M1"];
cs_ind=table(); cs_group=table();
for is=1:numel(subject_ids)
 sid=subject_ids(is);
 for hb=hbs
  T=all_rows(all_rows.subject_id==sid&all_rows.hemoglobin==hb&all_rows.channel_type=="Long",:);
  for k=1:numel(systems)
   x=T.contrast_beta(T.roi==systems(k));x=x(isfinite(x));cs_ind=[cs_ind;table(sid,hb,systems(k),numel(x),mean(x),median(x),'VariableNames',{'subject_id','hemoglobin','cortical_system','n_long_channels','mean_contrast_beta','median_contrast_beta'})]; %#ok<AGROW>
  end
 end
end
for hb=hbs
 H=cs_ind(cs_ind.hemoglobin==hb,:);pvals=nan(numel(systems),1);rows=table();
 for k=1:numel(systems)
  x=H.mean_contrast_beta(H.cortical_system==systems(k));[tt,p]=one_sample_t(x);pvals(k)=p;rows=[rows;table(hb,systems(k),numel(x),mean(x),std(x),std(x)/sqrt(numel(x)),tt,p,sum(x>0),sum(x<0),'VariableNames',{'hemoglobin','cortical_system','n_subjects','mean_contrast_beta','sd_contrast_beta','se_contrast_beta','t_statistic','p_value','n_positive_beta','n_negative_beta'})]; %#ok<AGROW>
 end
 rows.q_value_bh_fdr=bh_fdr(pvals);cs_group=[cs_group;rows];
end
writetable(cs_ind,fullfile(figroot,char(prefix+"_individual_cortical_system_long45_n29.csv")));
writetable(cs_group,fullfile(figroot,char(prefix+"_group_cortical_system_long45_n29.csv")));
plot_figures(ind,roi_group,channel_group,cs_ind,cs_group,figroot,n,keys,systems,hbs,phase_key);
out=struct('output_root',figroot,'n_subjects',n,'channel_summary',channel_group,'roi_summary',roi_group,'cortical_summary',cs_group);
fprintf('Q9 %s group analysis complete: n=%d, outputs=%s\n',upper(char(phase_key)),n,figroot);
end

function [t,p]=one_sample_t(x)
x=x(isfinite(x)); n=numel(x); if n<2||std(x)==0,t=NaN;p=NaN;return;end;t=mean(x)/(std(x)/sqrt(n));p=betainc(n/(n+t^2),n/2,.5);
end
function q=bh_fdr(p)
q=nan(size(p));v=isfinite(p);pv=p(v);[s,o]=sort(pv(:));m=numel(s);a=min(1,flipud(cummin(flipud(s.*m./(1:m)'))));z=zeros(size(pv));z(o)=a;ii=find(v);q(ii)=z;
end
function plot_figures(ind,roi,chan,csind,csg,root,n,keys,systems,hbs,phase_key)
out=fullfile(root,'figures');if ~isfolder(out),mkdir(out);end;set(0,'DefaultFigureVisible','off'); figprefix="q9_"+phase_key;
f=figure('Color','w','Position',[100 100 1800 1050]);tiledlayout(2,1,'TileSpacing','compact');
for j=1:2,H=ind(ind.hemoglobin==hbs(j),:);M=nan(n,height(keys));for k=1:height(keys),V=H(H.roi==keys.roi(k)&H.laterality==keys.laterality(k),:);[~,o]=sort(V.subject_id);M(:,k)=V.mean_contrast_beta(o);end;lim=max(abs(M(:)),[],'omitnan');if ~isfinite(lim)||lim==0,lim=1;end;nexttile;imagesc(M,[-lim lim]);colormap(gca,bwr(256));colorbar;set(gca,'XTick',1:height(keys),'XTickLabel',strcat(keys.roi," | ",keys.laterality),'XTickLabelRotation',35,'YTick',1:n,'YTickLabel',compose('%03d',sort(str2double(unique(H.subject_id)))),'FontSize',9,'XColor','k','YColor','k');xlabel('ROI | laterality');ylabel('Subject');title(sprintf('Q9 %s individual ROI beta | %s | n=%d',upper(char(phase_key)),hbs(j),n));end
sgtitle(sprintf('Q9 %s: individual ROI profiles | high intrinsic - ExperimentalBaseline | n=%d',upper(char(phase_key)),n));exportgraphics(f,fullfile(out,char(figprefix+'_individual_ROI_heatmaps_long45_allQC_n29.png')),'Resolution',220);close(f);
f=figure('Color','w','Position',[100 100 1800 1000]);tiledlayout(1,2);for j=1:2,U=roi(roi.hemoglobin==hbs(j),:);[~,o]=sort(U.mean_contrast_beta);U=U(o,:);y=1:height(U);nexttile;errorbar(U.mean_contrast_beta,y,U.se_contrast_beta,'horizontal','o','LineWidth',1.2);hold on;xline(0,'k-');sig=U.q_value_bh_fdr<.05;plot(U.mean_contrast_beta(sig),y(sig),'rp','MarkerFaceColor','r');set(gca,'YDir','reverse','YTick',y,'YTickLabel',strcat(U.roi," | ",U.laterality),'XColor','k','YColor','k');xlabel('Mean beta +/- SE');ylabel('ROI | laterality');title(sprintf('%s group ROI effects | n=%d',hbs(j),n));grid on;end;sgtitle(sprintf('Q9 %s group ROI effects | BH-FDR q<.05 marker | n=%d',upper(char(phase_key)),n));exportgraphics(f,fullfile(out,char(figprefix+'_group_ROI_effects_long45_allQC_n29.png')),'Resolution',220);close(f);
f=figure('Color','w','Position',[100 100 1900 700]);tiledlayout(2,1);for j=1:2,U=chan(chan.hemoglobin==hbs(j)&chan.channel_type=="Long",:);[~,o]=sort(U.channel_index);U=U(o,:);v=U.mean_contrast_beta';lim=max(abs(v),[],'omitnan');if ~isfinite(lim)||lim==0,lim=1;end;nexttile;imagesc(v,[-lim lim]);colormap(gca,bwr(256));colorbar;set(gca,'YTick',1,'YTickLabel',hbs(j),'XTick',1:height(U),'XTickLabel',compose('ch%d',U.channel_index),'XTickLabelRotation',45,'XColor','k','YColor','k');xlabel('Long channel');title(sprintf('%s group mean beta; dot = uncorrected p<.05 | n=%d',hbs(j),n));hold on;sig=U.p_value<.05;plot(find(sig),ones(sum(sig),1),'ko','MarkerFaceColor','k');end;sgtitle(sprintf('Q9 %s 45 long-channel group effects | n=%d',upper(char(phase_key)),n));exportgraphics(f,fullfile(out,char(figprefix+'_group_long_channel_effects_long45_allQC_n29.png')),'Resolution',220);close(f);
f=figure('Color','w','Position',[100 100 1600 850]);tiledlayout(1,2);for j=1:2,U=csg(csg.hemoglobin==hbs(j),:);[~,o]=ismember(systems,U.cortical_system);U=U(o,:);y=1:height(U);nexttile;errorbar(U.mean_contrast_beta,y,U.se_contrast_beta,'horizontal','o','LineWidth',1.2);hold on;xline(0,'k-');sig=U.q_value_bh_fdr<.05;plot(U.mean_contrast_beta(sig),y(sig),'rp','MarkerFaceColor','r');set(gca,'YDir','reverse','YTick',y,'YTickLabel',U.cortical_system,'XColor','k','YColor','k');xlabel('Mean beta +/- SE');ylabel('Cortical system');title(sprintf('%s cortical-system effects | n=%d',hbs(j),n));grid on;end;sgtitle(sprintf('Q9 %s cortical-system group effects | BH-FDR q<.05 marker | n=%d',upper(char(phase_key)),n));exportgraphics(f,fullfile(out,char(figprefix+'_group_cortical_system_effects_long45_n29.png')),'Resolution',220);close(f);
end
function c=bwr(n),h=floor(n/2);c=[[linspace(.1,1,h)' linspace(.3,1,h)' ones(h,1)];[ones(n-h,1) linspace(1,.15,n-h)' linspace(1,.15,n-h)']];end
