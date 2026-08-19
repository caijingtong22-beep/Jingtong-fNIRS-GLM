function output_dir = plot_q10_roi_long45_allqc()
% Plot the three Q10 analysis layers for all 45 long channels (all QC flags kept).

this_file = mfilename('fullpath');
model_root = fileparts(this_file);
q10_root = fullfile(model_root, 'q10_hiloadhe_vs_hiloadle');
output_dir = fullfile(q10_root, 'figures_long45_allQC');
if ~isfolder(output_dir), mkdir(output_dir); end

I = readtable(fullfile(q10_root, 'q10_HE_minus_LE_individual_ROI_long45_allQC_n30.csv'), 'TextType', 'string');
G = readtable(fullfile(q10_root, 'q10_HE_minus_LE_group_ROI_long45_allQC_n30.csv'), 'TextType', 'string');
C = readtable(fullfile(q10_root, 'q10_HE_minus_LE_group_channel_long45_allQC_n30.csv'), 'TextType', 'string');
hbs = ["HbO", "HbR"];
roi_keys = unique(I(I.hemoglobin == "HbO", {'roi','laterality'}), 'rows', 'stable');
roi_labels = strcat(roi_keys.roi, " | ", roi_keys.laterality);

%% Layer 1: individual ROI mean-beta heatmaps
fig = figure('Color', 'w', 'Position', [100 100 1800 1050]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for hi = 1:2
    hb = hbs(hi);
    U = I(I.hemoglobin == hb, :);
    M = nan(30, height(roi_keys));
    for r = 1:height(roi_keys)
        V = U(U.roi == roi_keys.roi(r) & U.laterality == roi_keys.laterality(r), :);
        [~, order] = sort(V.subject_id);
        M(:, r) = V.mean_contrast_beta(order);
    end
    nexttile;
    lim = max(abs(M(:)), [], 'omitnan');
    imagesc(M, [-lim lim]); axis tight;
    colormap(gca, blue_white_red(256)); colorbar;
    set(gca, 'XTick', 1:height(roi_keys), 'XTickLabel', roi_labels, ...
        'XTickLabelRotation', 35, 'YTick', 1:30, ...
        'YTickLabel', compose('%03d', 1:30), 'FontSize', 9, ...
        'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    xlabel('Broad ROI | laterality'); ylabel('Subject');
    title(sprintf('%s: subject-level ROI mean beta (HE - LE)', hb));
end
sgtitle('Q10 layer 1: individual ROI profiles | all 45 long channels; QC flags retained');
exportgraphics(fig, fullfile(output_dir, 'q10_layer1_individual_roi_heatmaps_long45_allQC.png'), 'Resolution', 220);
close(fig);

%% Layer 2: group ROI mean plus/minus SE
fig = figure('Color', 'w', 'Position', [100 100 1800 1050]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for hi = 1:2
    hb = hbs(hi);
    U = G(G.hemoglobin == hb, :);
    labels = strcat(U.roi, " | ", U.laterality);
    [~, order] = sort(U.mean_roi_contrast_beta);
    U = U(order, :); labels = labels(order);
    nexttile; hold on;
    y = 1:height(U);
    errorbar(U.mean_roi_contrast_beta, y, U.group_se, 'horizontal', 'o', ...
        'LineWidth', 1.2, 'MarkerFaceColor', [0.1 0.3 0.55], 'Color', [0.1 0.3 0.55]);
    xline(0, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    set(gca, 'YDir', 'reverse', 'YTick', y, 'YTickLabel', labels, 'FontSize', 10, ...
        'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.82 0.82 0.82]);
    xlabel('Group mean beta (HE - LE) +/- SE'); ylabel('ROI | laterality');
    title(sprintf('%s: 30-subject ROI effects (no BH-FDR result)', hb)); grid on;
end
sgtitle('Q10 layer 2: group ROI effects | all 45 long channels; QC flags retained');
exportgraphics(fig, fullfile(output_dir, 'q10_layer2_group_roi_effects_long45_allQC.png'), 'Resolution', 220);
close(fig);

%% Layer 3: long-channel group beta heatmap and uncorrected p markers
fig = figure('Color', 'w', 'Position', [100 100 1900 700]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for hi = 1:2
    hb = hbs(hi);
    U = C(C.hemoglobin == hb, :);
    [~, order] = sort(U.channel_id); U = U(order, :);
    nexttile; hold on;
    values = U.mean_contrast_beta';
    lim = max(abs(values), [], 'omitnan');
    imagesc(values, [-lim lim]); colormap(gca, blue_white_red(256));
    set(gca, 'YTick', 1, 'YTickLabel', hb, 'XTick', 1:height(U), ...
        'XTickLabel', compose('ch%d', U.channel_id), 'XTickLabelRotation', 45, 'FontSize', 9, ...
        'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    xlabel('Long channel'); ylabel('Hemoglobin');
    title(sprintf('%s: group mean channel beta (HE - LE); black dot = uncorrected p < .05; none survives BH-FDR', hb));
    sig = U.group_p_two_sided < 0.05;
    plot(find(sig), ones(sum(sig), 1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    colorbar;
end
sgtitle('Q10 layer 3: 45 long-channel group effects | QC flags retained');
exportgraphics(fig, fullfile(output_dir, 'q10_layer3_long_channel_effects_long45_allQC.png'), 'Resolution', 220);
close(fig);

fprintf('Saved Q10 long45 figures to %s\n', output_dir);
end

function cmap = blue_white_red(n)
if nargin < 1, n = 256; end
half = floor(n/2);
blue = [linspace(0.1, 1, half)' linspace(0.3, 1, half)' ones(half,1)];
red = [ones(n-half,1) linspace(1, 0.15, n-half)' linspace(1, 0.15, n-half)'];
cmap = [blue; red];
end
