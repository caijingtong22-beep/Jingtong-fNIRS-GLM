function output_dir = plot_q10_cortical_system_long45_allqc()
% Publication-style Q10 figures grouped by prefrontal versus motor cortex.

this_file = mfilename('fullpath');
model_root = fileparts(this_file);
q10_root = fullfile(model_root, 'q10_hiloadhe_vs_hiloadle');
output_dir = fullfile(q10_root, 'figures_cortical_system_long45_allQC');
if ~isfolder(output_dir), mkdir(output_dir); end

I = readtable(fullfile(q10_root, 'q10_HE_minus_LE_individual_ROI_cortical_system_long45_allQC_n30.csv'), 'TextType', 'string');
G = readtable(fullfile(q10_root, 'q10_HE_minus_LE_group_ROI_cortical_system_long45_allQC_n30.csv'), 'TextType', 'string');
C = readtable(fullfile(q10_root, 'q10_HE_minus_LE_group_channel_long45_allQC_n30.csv'), 'TextType', 'string');
roi_order = ["Medial Frontopolar PFC", "Lateral Frontopolar PFC", "DL PFC", ...
    "VL PFC", "Medial PFC", "PMC", "SMA", "M1"];
roi_short = ["mFPFC", "lFPFC", "dlPFC", "vlPFC", "mPFC", "PMC", "SMA", "M1"];
systems = ["Prefrontal Cortex", "Prefrontal Cortex", "Prefrontal Cortex", ...
    "Prefrontal Cortex", "Prefrontal Cortex", "Motor Cortex", "Motor Cortex", "Motor Cortex"];
hbs = ["HbO", "HbR"];

%% Layer 1: individual heatmaps
fig = figure('Color', 'w', 'Position', [100 100 1500 1050]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for hi = 1:2
    hb = hbs(hi); U = I(I.hemoglobin == hb, :); M = nan(30, 8);
    for r = 1:8
        V = U(U.roi == roi_order(r), :); [~, ord] = sort(V.subject_id);
        M(:,r) = V.mean_contrast_beta(ord);
    end
    nexttile; lim = max(abs(M(:)), [], 'omitnan');
    imagesc(M, [-lim lim]); colormap(gca, blue_white_red(256));
    cb = colorbar; cb.Label.String = 'ROI mean beta: HE - LE';
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 11, ...
        'XTick', 1:8, 'XTickLabel', roi_short, 'YTick', 1:30, ...
        'YTickLabel', compose('%03d', 1:30));
    ylabel('Subject'); title(sprintf('%s: each cell is one subject''s ROI mean beta', hb), 'Color', 'k');
    xline(5.5, 'k-', 'LineWidth', 1.5);
    text(3, -1.6, 'Prefrontal Cortex', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(6.9, -1.6, 'Motor Cortex', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
sgtitle({'Q10 layer 1: individual ROI profiles (all 45 long channels; QC retained)', ...
    'Color: red = HE > LE; blue = LE > HE; white = approximately no difference'}, 'Color', 'k', 'FontWeight', 'bold');
exportgraphics(fig, fullfile(output_dir, 'q10_layer1_individual_roi_cortical_system_long45_allQC.png'), 'Resolution', 220);
close(fig);

%% Layer 2: bar chart as cortical-system grouping
fig = figure('Color', 'w', 'Position', [100 100 1600 720]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
bar_colors = {[0.86 0.18 0.20], [0.12 0.35 0.72]};
for hi = 1:2
    hb = hbs(hi); U = G(G.hemoglobin == hb, :); [~, ord] = ismember(roi_order, U.roi); U = U(ord,:);
    nexttile; hold on;
    bar(1:8, U.mean_roi_contrast_beta, 0.68, 'FaceColor', bar_colors{hi}, 'EdgeColor', 'none');
    errorbar(1:8, U.mean_roi_contrast_beta, U.group_se, 'k.', 'LineWidth', 1.1, 'CapSize', 8);
    yline(0, 'k-', 'LineWidth', 1);
    xline(5.5, 'k-', 'LineWidth', 1.2);
    ylim([-1.6 1.6]); xlim([0.35 8.65]);
    set(gca, 'Color','w','XColor','k','YColor','k','FontSize',11, ...
        'XTick',1:8,'XTickLabel',roi_short);
    ylabel('Group mean beta (HE - LE)');
    title(sprintf('%s: bar = mean of 30 subjects; whisker = SEM', hb), 'Color', 'k');
    text(3, 1.72, 'Prefrontal Cortex', 'HorizontalAlignment','center','FontWeight','bold');
    text(6.9, 1.72, 'Motor Cortex', 'HorizontalAlignment','center','FontWeight','bold');
end
sgtitle({'Q10 layer 2: group ROI effects (all 45 long channels; QC retained)', ...
    'No ROI survives BH-FDR across the 8 ROI groups for either HbO or HbR'}, 'Color', 'k', 'FontWeight', 'bold');
exportgraphics(fig, fullfile(output_dir, 'q10_layer2_group_roi_cortical_system_long45_allQC.png'), 'Resolution', 220);
close(fig);

%% Layer 3: channel-level heatmaps sorted by cortical system and ROI
fig = figure('Color', 'w', 'Position', [100 100 1700 800]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
for hi = 1:2
    hb = hbs(hi); U = C(C.hemoglobin == hb, :);
    sort_key = zeros(height(U), 1);
    for r = 1:8, sort_key(U.roi == roi_order(r)) = r; end
    [~, ord] = sortrows([sort_key, U.channel_id]); U = U(ord,:); sort_key = sort_key(ord);
    lim = max(abs(U.mean_contrast_beta), [], 'omitnan');
    nexttile; imagesc(U.mean_contrast_beta', [-lim lim]); colormap(gca, blue_white_red(256));
    cb = colorbar; cb.Label.String = 'Group mean channel beta: HE - LE';
    set(gca,'Color','w','XColor','k','YColor','k','FontSize',10, ...
        'XTick',1:height(U),'XTickLabel',compose('ch%d',U.channel_id), ...
        'XTickLabelRotation',45,'YTick',1,'YTickLabel',hb);
    sig = U.group_p_two_sided < 0.05;
    hold on; plot(find(sig), ones(sum(sig),1), 'ko', 'MarkerFaceColor','k','MarkerSize',6);
    boundaries = find(diff(sort_key) ~= 0);
    for b = boundaries', xline(b + 0.5, 'k-', 'LineWidth', 1); end
    title(sprintf('%s: black dot = uncorrected p < .05 (none survives BH-FDR)', hb), 'Color','k');
end
sgtitle({'Q10 layer 3: long-channel group effects grouped by cortical system / ROI', ...
    'Color: red = HE > LE; blue = LE > HE; black dot = uncorrected p < .05; no dot is FDR-significant'}, ...
    'Color','k','FontWeight','bold');
exportgraphics(fig, fullfile(output_dir, 'q10_layer3_long_channels_cortical_system_long45_allQC.png'), 'Resolution', 220);
close(fig);

fprintf('Saved revised cortical-system figures to %s\n', output_dir);
end

function cmap = blue_white_red(n)
half = floor(n/2);
cmap = [[linspace(0.1,1,half)' linspace(0.3,1,half)' ones(half,1)]; ...
        [ones(n-half,1) linspace(1,0.15,n-half)' linspace(1,0.15,n-half)']];
end
