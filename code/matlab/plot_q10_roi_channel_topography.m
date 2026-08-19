function plot_q10_roi_channel_topography()
%PLOT_Q10_ROI_CHANNEL_TOPOGRAPHY
% Q10 group beta schematic on a two-dimensional, allocation-based montage.
% It is a companion figure to NIRS-SPM maps, not an atlas-derived cortical map.

    this_file = mfilename('fullpath');
    q10_dir = fullfile(fileparts(this_file), 'q10_hiloadhe_vs_hiloadle');
    registration_dir = fullfile(q10_dir, '..', '..', 'group_outputs', 'NIRS_SPM', ...
        'channel_registration_q10_validated');
    result_csv = fullfile(q10_dir, 'q10_HE_minus_LE_group_channel_long45_allQC_n30.csv');
    coord_csv = fullfile(registration_dir, 'q10_validated_channel_mni_coordinates.csv');
    output_dir = fullfile(q10_dir, 'figures');

    assert(isfile(result_csv), 'Missing Q10 group channel table: %s', result_csv);
    assert(isfile(coord_csv), 'Missing validated coordinate mapping: %s', coord_csv);
    if ~isfolder(output_dir), mkdir(output_dir); end

    result = readtable(result_csv, 'TextType', 'string');
    coord = readtable(coord_csv, 'TextType', 'string');
    result.channel_id = double(result.channel_id);
    coord.channel_id = double(coord.channel_id);

    long_coord = coord(strcmpi(coord.channel_type, 'Long') & coord.rendered_geometry == 1, :);
    assert(height(long_coord) == 45, 'Expected 45 rendered long channels.');

    [is_found, loc] = ismember(result.channel_id, long_coord.channel_id);
    assert(all(is_found), 'At least one Q10 result channel is absent from the validated long-channel mapping.');
    data = result;
    data.x = long_coord.mni_x_unprojected(loc);
    data.y = long_coord.mni_y_unprojected(loc);
    data.roi = long_coord.roi(loc);

    hbo = data(strcmpi(data.hemoglobin, 'HbO'), :);
    hbr = data(strcmpi(data.hemoglobin, 'HbR'), :);
    assert(height(hbo) == 45 && height(hbr) == 45, 'Expected 45 HbO and 45 HbR rows.');

    beta_limit = max(abs([hbo.mean_contrast_beta; hbr.mean_contrast_beta]));
    beta_limit = max(1, ceil(beta_limit * 2) / 2);

    figure_handle = figure('Color', 'w', 'Position', [80 80 1600 840], 'Renderer', 'painters');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    draw_panel(hbo, beta_limit, 'HbO');
    draw_panel(hbr, beta_limit, 'HbR');
    main_title = sgtitle({'Q10: group mean contrast beta (HILoadHE - HILoadLE)', ...
        '45 long channels; n = 30; dots are measured channels; background is a display interpolation'}, ...
        'FontWeight', 'bold', 'FontSize', 16, 'Color', 'k');

    annotation(figure_handle, 'textbox', [0.12 0.018 0.76 0.045], 'String', ...
        'Same-colour channel rings and ROI labels denote allocation-based channel groups only; they do not represent atlas-derived cortical boundaries.', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'k');

    png_file = fullfile(output_dir, 'q10_group_roi_channel_topography_long45_n30.png');
    pdf_file = fullfile(output_dir, 'q10_group_roi_channel_topography_long45_n30.pdf');
    exportgraphics(figure_handle, png_file, 'Resolution', 300);
    exportgraphics(figure_handle, pdf_file, 'ContentType', 'vector');
    close(figure_handle);

    summary = table(string(png_file), string(pdf_file), beta_limit, ...
        'VariableNames', {'png_file', 'pdf_file', 'symmetric_beta_colour_limit'});
    writetable(summary, fullfile(output_dir, 'q10_group_roi_channel_topography_manifest.csv'));
    fprintf('Wrote Q10 ROI/channel topography:\n%s\n%s\n', png_file, pdf_file);
end

function draw_panel(data, beta_limit, chromophore)
    ax = nexttile;
    hold(ax, 'on');
    axis(ax, 'equal');
    axis(ax, 'off');

    % Coordinates are MNI x (left-right) and MNI y (posterior-anterior).
    % The head outline is schematic; it does not claim individual anatomy.
    all_x = data.x;
    all_y = data.y;
    x_center = mean([min(all_x), max(all_x)]);
    y_center = mean([min(all_y), max(all_y)]) - 8;
    radius_x = (max(all_x) - min(all_x)) / 2 + 18;
    radius_y = (max(all_y) - min(all_y)) / 2 + 22;
    theta = linspace(0, 2*pi, 300);
    plot(ax, x_center + radius_x*cos(theta), y_center + radius_y*sin(theta), 'k-', 'LineWidth', 1.5);
    plot(ax, [x_center-7 x_center x_center+7], [y_center+radius_y-1 y_center+radius_y+10 y_center+radius_y-1], 'k-', 'LineWidth', 1.5);
    plot(ax, x_center-radius_x+[-2 -10 -2], y_center+[-10 0 10], 'k-', 'LineWidth', 1.2);
    plot(ax, x_center+radius_x+[2 10 2], y_center+[-10 0 10], 'k-', 'LineWidth', 1.2);

    % Gentle interpolated background, clipped visually by the head outline.
    grid_x = linspace(x_center-radius_x, x_center+radius_x, 160);
    grid_y = linspace(y_center-radius_y, y_center+radius_y, 160);
    [X, Y] = meshgrid(grid_x, grid_y);
    values = zeros(size(X)); weights = zeros(size(X));
    for idx = 1:height(data)
        distance2 = (X-data.x(idx)).^2 + (Y-data.y(idx)).^2;
        weight = exp(-distance2/(2*13^2));
        values = values + weight * data.mean_contrast_beta(idx);
        weights = weights + weight;
    end
    Z = values ./ max(weights, eps);
    in_head = ((X-x_center)/radius_x).^2 + ((Y-y_center)/radius_y).^2 <= 1;
    Z(~in_head) = NaN;
    contourf(ax, X, Y, Z, 18, 'LineStyle', 'none');
    colormap(ax, blue_white_red(256));
    clim(ax, [-beta_limit beta_limit]);

    % Broad ROI names are text annotations only.  No ROI outline is drawn,
    % so this panel does not imply a precise anatomical parcel or nesting.
    scatter(ax, data.x, data.y, 86, data.mean_contrast_beta, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
    roi_names = unique(data.roi, 'stable');
    for r = 1:numel(roi_names)
        subset = data(data.roi == roi_names(r), :);
        centre = [mean(subset.x), mean(subset.y)];
        roi_colour = roi_colour_for_name(roi_names(r));
        % The thin, unfilled ring encodes only ROI membership. The central
        % marker fill remains the contrast beta and is not overwritten.
        scatter(ax, subset.x, subset.y, 108, 'o', 'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', roi_colour, 'LineWidth', 1.6);
        text(ax, centre(1), centre(2), short_roi_name(roi_names(r)), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'Color', roi_colour, ...
            'BackgroundColor', 'w', 'Margin', 1, 'FontSize', 9);
    end
    for idx = 1:height(data)
        text(ax, data.x(idx), data.y(idx)-2.8, sprintf('%d', data.channel_id(idx)), ...
            'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold', 'Color', 'k');
    end
    title(ax, sprintf('%s: mean beta, HE - LE', chromophore), 'FontWeight', 'bold', 'Color', 'k');
    xlim(ax, x_center + [-radius_x-12 radius_x+12]);
    ylim(ax, y_center + [-radius_y-8 radius_y+19]);
    cb = colorbar(ax);
    cb.Label.String = 'Mean contrast beta (\muM)';
    cb.Label.FontSize = 10;
    cb.Label.Color = 'k';
    cb.Color = 'k';
    text(ax, x_center, y_center-radius_y-4, 'Posterior', 'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'k');
    text(ax, x_center, y_center+radius_y+13, 'Anterior / nose', 'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'k');
end

function label = short_roi_name(roi)
    switch string(roi)
        case "Medial Frontopolar PFC", label = 'mFPC';
        case "Lateral Frontopolar PFC", label = 'lFPFC';
        case {"Dorsolateral PFC", "DL PFC"}, label = 'DL PFC';
        case {"Ventrolateral PFC", "VL PFC"}, label = 'VL PFC';
        case "Medial PFC", label = 'mPFC';
        case {"Premotor Cortex", "PMC"}, label = 'PMC';
        case {"Supplementary Motor Area", "SMA"}, label = 'SMA';
        case {"Primary Motor Cortex", "M1"}, label = 'M1';
        otherwise, label = char(roi);
    end
end

function colour = roi_colour_for_name(roi)
    % Stable, high-contrast ROI palette. These colours are used only for
    % channel membership rings and ROI text, never for beta magnitude.
    switch string(roi)
        case "Medial Frontopolar PFC", colour = [0.00 0.35 0.70]; % blue
        case "Lateral Frontopolar PFC", colour = [0.85 0.33 0.10]; % vermilion
        case {"Dorsolateral PFC", "DL PFC"}, colour = [0.45 0.15 0.65]; % purple
        case {"Ventrolateral PFC", "VL PFC"}, colour = [0.10 0.55 0.45]; % green
        case "Medial PFC", colour = [0.75 0.15 0.35]; % magenta-red
        case {"Premotor Cortex", "PMC"}, colour = [0.05 0.55 0.75]; % cyan-blue
        case {"Supplementary Motor Area", "SMA"}, colour = [0.55 0.40 0.00]; % ochre
        case {"Primary Motor Cortex", "M1"}, colour = [0.30 0.30 0.30]; % charcoal
        otherwise, colour = [0.20 0.20 0.20];
    end
end

function cmap = blue_white_red(n)
    half = floor(n/2);
    blue = [linspace(0.08,1,half)' linspace(0.25,1,half)' ones(half,1)];
    red_count = n-half;
    red = [ones(red_count,1) linspace(1,0.10,red_count)' linspace(1,0.10,red_count)'];
    cmap = [blue; red];
end
