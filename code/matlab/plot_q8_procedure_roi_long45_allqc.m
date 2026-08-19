function output_dir = plot_q8_procedure_roi_long45_allqc()
% Three Q10-style figure layers for Q8 Procedure, HIHE - LILE.
root = fullfile('C:\aliceproject\Jingtong-fNIRS-GLM', 'results', 'q8_hihe_vs_lile_procedure', 'cohort_specification');
output_dir = fullfile(root, 'figures_long45_allQC'); if ~isfolder(output_dir), mkdir(output_dir); end
I=readtable(fullfile(root,'q8_HIHE_minus_LILE_procedure_individual_ROI_long45_allQC_n30.csv'),'TextType','string');
G=readtable(fullfile(root,'q8_HIHE_minus_LILE_procedure_group_ROI_long45_allQC_n30.csv'),'TextType','string');
C=readtable(fullfile(root,'q8_HIHE_minus_LILE_procedure_group_channel_long45_allQC_n30.csv'),'TextType','string');
hbs=["HbO","HbR"]; keys=unique(I(I.hemoglobin=="HbO",{'roi','laterality'}),'rows','stable'); labels=strcat(keys.roi," | ",keys.laterality);

fig=figure('Color','w','Position',[100 100 1800 1050]); set(fig,'DefaultTextColor','k'); tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
for hi=1:2
    U=I(I.hemoglobin==hbs(hi),:); M=nan(30,height(keys));
    for r=1:height(keys), V=U(U.roi==keys.roi(r) & U.laterality==keys.laterality(r),:); [~,o]=sort(V.subject_id); M(:,r)=V.mean_contrast_beta(o); end
    nexttile; lim=max(abs(M(:)),[],'omitnan'); if ~isfinite(lim)||lim==0,lim=1;end; imagesc(M,[-lim lim]); colormap(gca,blue_white_red(256)); cb=colorbar; cb.Color='k';
    set(gca,'XTick',1:height(keys),'XTickLabel',labels,'XTickLabelRotation',35,'YTick',1:30,'YTickLabel',compose('%03d',1:30),'FontSize',9,'Color','w','XColor','k','YColor','k'); xlabel('Broad ROI | laterality','Color','k');ylabel('Subject','Color','k');title(sprintf('%s: subject-level ROI mean beta (HIHE - LILE)',hbs(hi)),'Color','k'); force_black_axis_text(gca);
end
sgtitle('Q8 Procedure layer 1: individual ROI profiles | all 45 long channels; QC flags retained','Color','k'); exportgraphics(fig,fullfile(output_dir,'q8_procedure_layer1_individual_roi_heatmaps_long45_allQC.png'),'Resolution',220);close(fig);

fig=figure('Color','w','Position',[100 100 1800 1050]); set(fig,'DefaultTextColor','k'); tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for hi=1:2
    U=G(G.hemoglobin==hbs(hi),:); [~,o]=sort(U.mean_contrast_beta);U=U(o,:); labs=strcat(U.roi," | ",U.laterality); nexttile;hold on;y=1:height(U);
    errorbar(U.mean_contrast_beta,y,U.se_contrast_beta,'horizontal','o','LineWidth',1.2,'MarkerFaceColor',[.1 .3 .55],'Color',[.1 .3 .55]);xline(0,'-','Color',[.4 .4 .4]);
    sig=U.q_value_bh_fdr<.05; plot(U.mean_contrast_beta(sig),y(sig),'rp','MarkerFaceColor','r','MarkerSize',9);set(gca,'YDir','reverse','YTick',y,'YTickLabel',labs,'FontSize',10,'Color','w','XColor','k','YColor','k');xlabel('Group mean beta (HIHE - LILE) +/- SE','Color','k');ylabel('ROI | laterality','Color','k');title(sprintf('%s: 30-subject ROI effects',hbs(hi)),'Color','k');force_black_axis_text(gca);grid on;
end
sgtitle('Q8 Procedure layer 2: group ROI effects | red pentagram = BH-FDR q < .05','Color','k');exportgraphics(fig,fullfile(output_dir,'q8_procedure_layer2_group_roi_effects_long45_allQC.png'),'Resolution',220);close(fig);

fig=figure('Color','w','Position',[100 100 1900 700]);set(fig,'DefaultTextColor','k');tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
for hi=1:2
    U=C(C.hemoglobin==hbs(hi),:);[~,o]=sort(U.channel);U=U(o,:); nexttile;hold on;v=U.mean_contrast_beta';lim=max(abs(v),[],'omitnan');if ~isfinite(lim)||lim==0,lim=1;end;imagesc(v,[-lim lim]);colormap(gca,blue_white_red(256));cb=colorbar;cb.Color='k';
    set(gca,'YTick',1,'YTickLabel',hbs(hi),'XTick',1:height(U),'XTickLabel',compose('ch%d',U.channel),'XTickLabelRotation',45,'FontSize',9,'Color','w','XColor','k','YColor','k');xlabel('Long channel','Color','k');ylabel('Hemoglobin','Color','k');title(sprintf('%s: group mean channel beta; black dot = uncorrected p < .05',hbs(hi)),'Color','k');force_black_axis_text(gca); sig=U.p_value<.05;plot(find(sig),ones(sum(sig),1),'ko','MarkerFaceColor','k','MarkerSize',5);
end
sgtitle('Q8 Procedure layer 3: 45 long-channel group effects | HIHE - LILE; QC flags retained','Color','k');exportgraphics(fig,fullfile(output_dir,'q8_procedure_layer3_long_channel_effects_long45_allQC.png'),'Resolution',220);close(fig);
fprintf('Saved Q8 Procedure long45 figures to %s\n',output_dir);
end
function cmap=blue_white_red(n)
half=floor(n/2);cmap=[[linspace(.1,1,half)' linspace(.3,1,half)' ones(half,1)];[ones(n-half,1) linspace(1,.15,n-half)' linspace(1,.15,n-half)']];
end
function force_black_axis_text(ax)
% Explicitly override MATLAB theme defaults for every visible text element.
ax.XColor = 'k'; ax.YColor = 'k';
ax.XLabel.Color = 'k'; ax.YLabel.Color = 'k'; ax.Title.Color = 'k';
end
