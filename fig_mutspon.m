<<<<<<< HEAD
function fig_mutspon()

clear all;
close all;

% Define the mutation rates we want to include
mut_rates = [1e-1, 1e-2, 1e-3, 1e-4];

%
%
% % Create a new figure for the combined plots
% fig1 = figure('Position', [100, 100, 1000, 800]);
% tlo1 = tiledlayout(4, 2, 'TileSpacing','tight','Padding','none');
% 
% % For each mutation rate
% for i = 1:length(mut_rates)
%     % Define folder path
%     folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/mutation_spon/mr_%.1g', mut_rates(i));
%     numericsFile = fullfile(folderName, 'numerics.mat');
% 
%     % Load the necessary data files
%     numerics_data = load(numericsFile);
% 
%     % Extract relevant data
%     numerics = numerics_data.numerics;
% 
%     % Get time vector
%     filename = 'Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat';
%     info = load(filename, 'info').info;
%     time_init = info.time;
%     max_time = size(numerics.tumor_num, 2);
%     time_vector = (time_init + (1:max_time)) / 10;
% 
%     % Left column: Plot tumor cell numbers
%     ax1 = nexttile(2 * i - 1);
% 
% % Plot tumor cell number with its components
% plot(time_vector, numerics.hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
% hold on;
% plot(time_vector, numerics.normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normoxic');
% plot(time_vector, numerics.tumor_num, 'b-', 'DisplayName', 'Tumor');
% hold off;
% 
% % Add mutation rate label at top-right in LaTeX
% text(0.75, 0.85, sprintf('$\\mu = 10^{-%d}$', i), ...
%     'Units', 'normalized', ...
%     'HorizontalAlignment', 'right', ...
%     'VerticalAlignment', 'top', ...
%     'FontSize', 10, ...
%     'Interpreter', 'latex');
% 
% if i == 4
%     xlabel('t', 'FontSize', 10);
% end
% ylabel('Cell Number', 'FontSize', 10);
% 
% % Set axis limits based on mutation rate
% xlim([time_init / 10, time_vector(end)]);
% 
% % Set tick labels
% xticks_vals = linspace(time_init / 10, time_vector(end), 6);
% xticks_labels = arrayfun(@(x) sprintf('%.2g', x), xticks_vals, 'UniformOutput', false);
% set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, 'FontSize', 8);
% 
% if i == 1
%     title(sprintf('Tumor Cell Number'), 'FontSize', 10);
% end
% if i == 4
%     legend('FontSize', 8, 'Location', 'best');
% end
% 
%     % Right column: Plot damage accumulation
%     ax2 = nexttile(2 * i);
% 
%     % Plot damage accumulation with error bars
%     plot(time_vector, numerics.dam_accum(1, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     hold on;
%     plot(time_vector, numerics.dam_accum(1, :) - numerics.dam_accum(2, :), 'b--', ...
%         'DisplayName', 'mean - std');
%     plot(time_vector, numerics.dam_accum(1, :) + numerics.dam_accum(2, :), 'b-.', ...
%         'DisplayName', 'mean + std');
% 
%     plot(time_vector, numerics.death_thres(1, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     plot(time_vector, numerics.death_thres(1, :) - numerics.death_thres(2, :), 'r--', ...
%         'DisplayName', 'mean - std');
%     plot(time_vector, numerics.death_thres(1, :) + numerics.death_thres(2, :), 'r-.', ...
%         'DisplayName', 'mean + std');
%     hold off;
% 
%     if i == 4
%         xlabel('t', 'FontSize', 10);
%     end
%     ylabel('Averaged Value', 'FontSize', 10);
% 
%     xlim([time_init / 10, time_vector(end)]);
%     % Set tick labels
%     xticks_vals = linspace(time_init / 10, time_vector(end), 6);
%     xticks_labels = arrayfun(@(x) sprintf('%.2g', x), xticks_vals, 'UniformOutput', false);
%     set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, 'FontSize', 8);
% 
%     if i == 1
%         title(sprintf('Average Damage Accumulation & \n Average Death Threshold'), 'FontSize', 10);
%     end
% 
%     if i == 4
%         legend('FontSize', 8, 'Location', 'northeast');
%     end
% 
%     % Add a text label for pr value (as a subfigure caption)
%     num_annotations = 4; % 或需要的总数
%     spacing = 0.95 / num_annotations;
%     annotation('textbox', [0, 0.945 - (i - 1) * spacing, 0.04, 0.04], ...
%         'String', sprintf('(%c)', 'a' + i - 1), ... % (a), (b), (c)
%         'FitBoxToText', 'on', ...
%         'EdgeColor', 'none', ...
%         'FontSize', 12);
% 
% end
% 
% % Set uniform background color
% set(fig1, 'Color', 'white');
% 
% % Save the combined figure
% savefolder = 'Figures/vascular/(epsilon = 1e-2)/mutation_spon';
% saveas(fig1, fullfile(savefolder, 'mutspon1.fig'), 'fig');
% saveas(fig1, fullfile(savefolder, 'mutspon1.svg'), 'svg');

%

% Create a new figure for the combined plots
fig2 = figure('Position', [100, 100, 1000, 900]);
tlo2 = tiledlayout(4, 2, 'TileSpacing','tight','Padding','none');

% For each mutation rate
for i = 1:length(mut_rates)
    % Define folder path
    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/mutation_spon/mr_%.1g', mut_rates(i));
    numericsFile = fullfile(folderName, 'numerics.mat');
    paramsFile = fullfile(folderName, 'params.mat');

    % Load the necessary data files
    numerics_data = load(numericsFile);
    params_data = load(paramsFile);

    % Extract relevant data
    numerics = numerics_data.numerics;
    params = params_data.params;

    % Get time vector
    filename = 'Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat';
    info = load(filename, 'info').info;
    time_init = info.time;
    max_time = size(numerics.tumor_num, 2);
    time_vector = (time_init + (1:max_time)) / 10;

    % Plot oxygen consumption histogram (left column)
    ax1 = nexttile(2 * i - 1);
    [X, Y] = meshgrid(linspace(0.5 * params.rhoo, 4 * params.rhoo, params.nbins), time_vector);
    s1 = surf(X, Y, numerics.oxygen_consumption_hist);

    % Improve appearance with better colormap and lighting
    colormap(gca, cool); % Use parula colormap (brighter than hot)
    s1.FaceAlpha = 0.8;    % Make surface more visible
    s1.EdgeColor = 'none'; % Remove grid lines for cleaner look

    xlabel('Oxygen Consumption', 'FontSize', 10);
    ylabel('Time', 'FontSize', 10);
    zlabel('Frequency', 'FontSize', 10);

    % Set X axis properties
    xlim([0.5 * params.rhoo, 4 * params.rhoo]);
    if i ~= 4
        ylim([time_init / 10, time_vector(end)]);
    else
        ylim([time_init / 10, 19.4]);
    end
    xticks_vals = 0.5 * params.rhoo : 0.35 * params.rhoo : 4 * params.rhoo;
    xticks_labels = arrayfun(@(x) sprintf('%.2f', x), xticks_vals, 'UniformOutput', false);
    if i ~= 4
        yticks_vals = linspace(time_init / 10, time_vector(end), 6);
        yticks_labels = arrayfun(@(x) sprintf('%.2g', x), yticks_vals, 'UniformOutput', false);
        set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, ...
            'YTick', yticks_vals, 'YTickLabel', yticks_labels, ...
            'FontSize', 8);
    else
        yticks_vals = linspace(time_init / 10, 19.4, 7);
        yticks_labels = arrayfun(@(x) sprintf('%.3g', x), yticks_vals, 'UniformOutput', false);
        set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, ...
            'YTick', yticks_vals, 'YTickLabel', yticks_labels, ...
            'FontSize', 8);
    end

    if i == 1
        title(sprintf('Oxygen Consumption'), 'FontSize', 12);
    end

    % Add mutation rate label at top-right in LaTeX
    text(0.75, 0.9, sprintf('$\\mu = 10^{-%d}$', i), ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 10, ...
    'Interpreter', 'latex');


    % Enhance lighting for better 3D appearance
    lighting phong;
    camlight('headlight');
    material shiny;
    view([-15, 50]);
    % colorbar;

    % Plot proliferation rate histogram (right column)
    ax2 = nexttile(2 * i);
    [X, Y] = meshgrid(linspace(0.5 * log(2) / params.age, 4 * log(2) / params.age, params.nbins), time_vector);
    s2 = surf(X, Y, numerics.proliferation_rate_hist);

    % Apply same visual improvements
    colormap(gca, cool);
    s2.FaceAlpha = 0.8;
    s2.EdgeColor = 'none';

    xlabel('Proliferation Rate', 'FontSize', 10);
    ylabel('Time', 'FontSize', 10);
    zlabel('Frequency', 'FontSize', 10);

    % Set X axis properties
    xlim([0.5 * log(2) / params.age, 4 * log(2) / params.age]);
    if i ~= 4
        ylim([time_init / 10, time_vector(end)]);
    else
        ylim([time_init / 10, 19.4]);
    end
    xticks_vals = 0.5 * log(2) / params.age : 0.35 * log(2) / params.age : 4 * log(2) / params.age;
    xticks_labels = arrayfun(@(x) sprintf('%.2f', x), xticks_vals, 'UniformOutput', false);
    if i ~= 4
        yticks_vals = linspace(time_init / 10, time_vector(end), 6);
        yticks_labels = arrayfun(@(x) sprintf('%.2g', x), yticks_vals, 'UniformOutput', false);
        set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, ...
            'YTick', yticks_vals, 'YTickLabel', yticks_labels, ...
            'FontSize', 8);
    else
        yticks_vals = linspace(time_init / 10, 19.4, 7);
        yticks_labels = arrayfun(@(x) sprintf('%.3g', x), yticks_vals, 'UniformOutput', false);
        set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, ...
            'YTick', yticks_vals, 'YTickLabel', yticks_labels, ...
            'FontSize', 8);
    end

    if i == 1
        title(sprintf('Proliferation Rate'), 'FontSize', 12);
    end

    lighting phong;
    camlight('headlight');
    material shiny;
    view([-15, 50]);
    % colorbar;

    % Add a text label for pr value (as a subfigure caption)
    num_annotations = 4; % 或需要的总数
    spacing = 0.95 / num_annotations;
    annotation('textbox', [0.01, 0.99 - i * spacing, 0.05, 0.05], ...
        'String', sprintf('(%c)', 'a' + i - 1), ... % (a), (b), (c)
        'FitBoxToText', 'on', ...
        'EdgeColor', 'none', ...
        'FontSize', 12);

end

% Set uniform background color
set(fig2, 'Color', 'white');

% Save the combined figure
savefolder = 'Figures/vascular/(epsilon = 1e-2)/mutation_spon';
saveas(fig2, fullfile(savefolder, 'mutspon2.fig'), 'fig');
saveas(fig2, fullfile(savefolder, 'mutspon2.svg'), 'svg');

%
% % 创建新图像
fig3 = figure('Position', [100, 100, 900, 1200]);
tlo3 = tiledlayout(4,3, 'TileSpacing','none','Padding','none');

% 遍历每个pr值
for col_idx = 1:3
    mut_value = mut_rates(col_idx);

    % 构建文件路径
    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/mutation_spon/mr_%.1g', mut_value);
    filePath = fullfile(folderName, 'vessel1.fig');

    orig_fig = openfig(filePath, 'invisible');
    all_axes = findobj(orig_fig, 'Type', 'axes');

    % 遍历每一行
    for row_idx = 1:4
        figure(fig3);
        new_subplot = nexttile(tlo3, (row_idx - 1) * 3 + col_idx);

        % 查找对应的原始子图
        found = false;
        for ax_idx = 1:length(all_axes)
            ax = all_axes(ax_idx);
            pos = get(ax, 'Position');

            if pos(1) < 0.3 && ((row_idx == 1 && pos(2) > 0.7) || ...
                    (row_idx == 2 && pos(2) > 0.45 && pos(2) <= 0.7) || ...
                    (row_idx == 3 && pos(2) > 0.2 && pos(2) <= 0.45) || ...
                    (row_idx == 4 && pos(2) <= 0.2))
                orig_subplot = ax;
                found = true;
                break;
            end
        end

        if found
            % 复制内容和属性
            children = get(orig_subplot, 'Children');
            for i = 1:length(children)
                copyobj(children(i), new_subplot);
            end

            % 改为（反转复制顺序）
            children = get(orig_subplot, 'Children');
            for i = length(children):-1:1
                copyobj(children(i), new_subplot);
            end

            % 复制轴属性
            set(new_subplot, 'XLim', get(orig_subplot, 'XLim'));
            set(new_subplot, 'YLim', get(orig_subplot, 'YLim'));
            set(new_subplot, 'XTick', get(orig_subplot, 'XTick'));
            set(new_subplot, 'YTick', get(orig_subplot, 'YTick'));
            set(new_subplot, 'FontSize', 10);

            % 使用简化方法创建图例 - 不依赖于原图的对象
            hold(new_subplot, 'on');

            % 完全脱离子图中的真实对象，创建虚拟图例项
            h1 = plot(new_subplot, NaN, NaN, 'bs', 'MarkerFaceColor', 'b', 'DisplayName', 'Vessel');
            h2 = plot(new_subplot, NaN, NaN, 'mo', 'MarkerFaceColor', 'm', 'DisplayName', 'Hypoxic');
            h3 = plot(new_subplot, NaN, NaN, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Normorxic');

            % 只使用需要的图例项
            leg_handles = [h1];
            leg_labels = {'Vessel'};

            % 根据原图中的图例字符串添加其他项
            orig_legend = legend(orig_subplot);
            if ~isempty(orig_legend)
                try
                    orig_legend_str = get(orig_legend, 'String');
                    if any(strcmp(orig_legend_str, 'Hypoxic'))
                        leg_handles = [leg_handles, h2];
                        leg_labels = [leg_labels, {'Hypoxic'}];
                    end
                    if any(strcmp(orig_legend_str, 'Normorxic'))
                        leg_handles = [leg_handles, h3];
                        leg_labels = [leg_labels, {'Normorxic'}];
                    end
                catch
                    % 如果无法获取原始图例，使用所有三个项
                    leg_handles = [h1, h2, h3];
                    leg_labels = {'Vessel', 'Hypoxic', 'Normorxic'};
                end
            end

            % 创建图例
            lgd = legend(new_subplot, leg_handles, leg_labels);
            set(lgd, 'FontSize', 6, 'Location', 'best');

            % 添加标签
            if col_idx == 1
                ylabel(new_subplot, sprintf('Time: t = %.1f', 14 + row_idx * 10), 'FontSize', 12);
            end
            if row_idx == 1
                title(sprintf('$\\mu = %.1g$', mut_value), 'FontSize', 10, 'Interpreter', 'latex')
            end
        end

        % 如果当前不是第一列，则去除 ytick 和 yticklabel
        if col_idx > 1
            set(new_subplot, 'YTickLabel', []);
        end

        % 如果当前不是最后一行，则去除 xtick 和 xticklabel
        if row_idx < 4
            set(new_subplot, 'XTickLabel', []);
        else
            xlabel(ax, 'x', 'FontSize', 10);
        end
    end
end

% 添加总标题和保存
sgtitle('Tumor Population Comparison for Different p_r Values', 'FontSize', 12);
set(fig3, 'Color', 'w');
% Save the combined figure
savefolder = 'Figures/vascular/(epsilon = 1e-2)/mutation_spon';
saveas(fig3, fullfile(savefolder, 'mutspon3.fig'), 'fig');
saveas(fig3, fullfile(savefolder, 'mutspon3.svg'), 'svg');

%% Figure4: 按照氧消耗率着色的肿瘤细胞
fig4 = figure('Position', [100, 100, 950, 1100]);
tlo4 = tiledlayout(4,3, 'TileSpacing','none','Padding','none');

for col_idx = 1:length(mut_rates) - 1
    mut_value = mut_rates(col_idx);

    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/mutation_spon/mr_%.1g', mut_value);
    filePath = fullfile(folderName, 'info.mat');
    info_data = load(filePath, 'info');
    info = info_data.info;
    params = info.params;

    for row_idx = 1:4
        ax = nexttile(tlo4, (row_idx - 1) * 3 + col_idx);

        tumor_cells = info.tumor_cells{row_idx};
        valid_tumor = ~cellfun('isempty', tumor_cells);

        if any(valid_tumor)
            % 提取细胞位置和氧消耗率（字段名称：oxygen_consumption_nor）
            pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
            pos_tumor = cat(1, pos_tumor{:});
            ocr = cellfun(@(tc) tc.oxygen_consumption_nor, tumor_cells(valid_tumor), 'UniformOutput', true);
            edges = linspace(0.5 * params.rhoo, 4 * params.rhoo, 11);
            bin_idx = discretize(ocr, edges);
            cmap = parula(10);  % 10级渐变色，可根据需要更换其他colormap

            scatter(pos_tumor(:,1), pos_tumor(:,2), 10, cmap(bin_idx, :), 'filled');
            colormap(cmap);
            clim([edges(1) edges(end)]);
            if col_idx == 3
                colorbar;
            end
            axis equal;
            xlabel('x', 'FontSize', 8);
            ylabel('y', 'FontSize', 8);
            if row_idx == 1
                title(sprintf('$\\mu = %.1g$', mut_value), 'FontSize', 10, 'Interpreter', 'latex');
            end
        end

        % 仅第一列显示 y-axis 的所有信息
        if col_idx == 1
            ylabel(ax, sprintf('Time: t = %.1f', 14 + row_idx * 10), 'FontSize', 12);
        else
            set(ax, 'YTickLabel', []);
        end

        % 仅最后一行显示 x-axis 的所有信息
        if row_idx < 4
            set(ax, 'XTickLabel', []);
        else
            xlabel(ax, 'x', 'FontSize', 12);
        end

        % 如果你希望刻度是 0:0.2:1，则在这里设置：
        set(ax, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
    end
end

set(fig4, 'Color', 'w');
% Save the combined figure
savefolder = 'Figures/vascular/(epsilon = 1e-2)/mutation_spon';
saveas(fig4, fullfile(savefolder, 'mutspon4.fig'), 'fig');
saveas(fig4, fullfile(savefolder, 'mutspon4.svg'), 'svg');

%% Figure5: 按照增殖率着色的肿瘤细胞
% 基准值 x = log(2)/params.age（未突变时的增殖率），控制范围为 [0.5*x, 4*x]，10等分
fig5 = figure('Position', [100, 100, 950, 1100]);
tlo5 = tiledlayout(4,3, 'TileSpacing','none','Padding','none');

for col_idx = 1:length(mut_rates) - 1
    mut_value = mut_rates(col_idx);

    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/mutation_spon/mr_%.1g', mut_value);
    filePath = fullfile(folderName, 'info.mat');
    info_data = load(filePath, 'info');
    info = info_data.info;
    params = info.params;

    for row_idx = 1:4
        ax = nexttile(tlo5, (row_idx - 1) * 3 + col_idx);

        tumor_cells = info.tumor_cells{row_idx};
        valid_tumor = ~cellfun('isempty', tumor_cells);

        if any(valid_tumor)
            pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
            pos_tumor = cat(1, pos_tumor{:});
            prate = cellfun(@(tc) tc.proliferation_rate, tumor_cells(valid_tumor));
            % 构造区间 [0.5*log(2)/params.age, 4*log(2)/params.age]并10等分
            edges = linspace(0.5 * log(2)/params.age, 4 * log(2)/params.age, 11);
            bin_idx = discretize(prate, edges);
            cmap = parula(10);

            scatter(pos_tumor(:,1), pos_tumor(:,2), 10, cmap(bin_idx, :), 'filled');
            colormap(cmap);
            clim([edges(1) edges(end)]);
            if col_idx == 3
                colorbar;
            end
            axis equal;
            xlabel('x', 'FontSize', 8);
            ylabel('y', 'FontSize', 8);
            if row_idx == 1
                title(sprintf('$\\mu = %.1g$', mut_value), 'FontSize', 10, 'Interpreter', 'latex');
            end
        end

        % 仅第一列显示 y-axis 的所有信息
        if col_idx == 1
            ylabel(ax, sprintf('Time: t = %.1f', 14 + row_idx * 10), 'FontSize', 12);
        else
            set(ax, 'YTickLabel', []);
        end

        % 仅最后一行显示 x-axis 的所有信息
        if row_idx < 4
            set(ax, 'XTickLabel', []);
        else
            xlabel(ax, 'x', 'FontSize', 12);
        end

        % 如果你希望刻度是 0:0.2:1，则在这里设置：
        set(ax, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
    end
end
set(fig5, 'Color', 'w');
% Save the combined figure
savefolder = 'Figures/vascular/(epsilon = 1e-2)/mutation_spon';
saveas(fig5, fullfile(savefolder, 'mutspon5.fig'), 'fig');
saveas(fig5, fullfile(savefolder, 'mutspon5.svg'), 'svg');

%
% 创建新图像
fig6 = figure('Position', [100, 100, 1250, 1150]);
tiledlayout(4, 4, 'TileSpacing','none','Padding','none');

% 遍历每个pr值
for col_idx = 1:length(mut_rates)
    mut_value = mut_rates(col_idx);

    % 构建文件路径
    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/mutation_spon/mr_%.1g', mut_value);
    info_data = load(fullfile(folderName, 'info.mat'));
    info = info_data.info;
    params = info.params;
    % filePath = fullfile(folderName, 'vessel1.fig');

    % orig_fig = openfig(filePath, 'invisible');
    % all_axes = findobj(orig_fig, 'Type', 'axes');

    % 遍历每一行
    for row_idx = 1:4
        figure(fig6);
        new_subplot = nexttile((row_idx - 1) * 4 + col_idx);

        drug = flipud(info.chem_field{row_idx}.drug);

        % Display the snapshot
        imagesc(drug);
         % Add colorbar to all columns but make non-last ones invisible
        if col_idx == 4
            colorbar;
        else
            cb = colorbar;
            cb.Visible = 'off'; % Invisible but reserves the space
        end
        axis tight equal;
        colormap(gca, 'parula');

        xlim([0, params.grid_size(2)]);
        ylim([0, params.grid_size(1)]);

        % Set axis labels and ticks
        if col_idx == 1
            ylabel(gca, sprintf('Time: t = %.1f', 14 + row_idx * 10), 'FontSize', 12);
        else
            set(new_subplot, 'YTickLabel', []);
        end

        if row_idx == 1
            title(sprintf('$\\mu = %.1g$', mut_value), 'FontSize', 10, 'Interpreter', 'latex')
        end

        if row_idx < 4
            set(new_subplot, 'XTickLabel', []);
        else
            xlabel('x', 'FontSize', 10);
        end

        if row_idx == 4
            set(gca, 'XTick', 0:params.grid_size(2) / 5:params.grid_size(2), 'XTickLabel', 0:0.2:1, ...
                'FontSize', 10);
        end
        if col_idx == 1
            set(gca, 'YTick', 0:params.grid_size(1) / 5:params.grid_size(1), 'YTickLabel', 1:-0.2:0, ...
                'FontSize', 10);
        end
    end

    % 如果当前不是第一列，则去除 ytick 和 yticklabel
    if col_idx > 1
        set(gca, 'YTickLabel', []);
    end
end
%
set(fig6, 'Color', 'w');
% % Save the combined figure
savefolder = 'Figures/vascular/(epsilon = 1e-2)/mutation_spon';
saveas(fig6, fullfile(savefolder, 'mutspon6.fig'), 'fig');
saveas(fig6, fullfile(savefolder, 'mutspon6.svg'), 'svg');


%% Figure7: 按照氧消耗率着色的肿瘤细胞
fig7 = figure('Position', [100, 100, 950, 1100]);
tiledlayout(4,3, 'TileSpacing','none','Padding','none');

for col_idx = 1:length(mut_rates) - 1
    mut_value = mut_rates(col_idx);

    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/mutation_spon/mr_%.1g', mut_value);
    filePath = fullfile(folderName, 'info.mat');
    info_data = load(filePath, 'info');
    info = info_data.info;
    % params = info.params;

    for row_idx = 1:4
        ax = nexttile((row_idx - 1) * 3 + col_idx);

        tumor_cells = info.tumor_cells{row_idx};
        valid_tumor = ~cellfun('isempty', tumor_cells);

        if any(valid_tumor)
            % 提取细胞位置和氧消耗率（字段名称：oxygen_consumption_nor）
            pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
            pos_tumor = cat(1, pos_tumor{:});
            ocr = cellfun(@(tc) tc.death_threshold, tumor_cells(valid_tumor));
            edges = linspace(min(ocr), max(ocr), 11);
            bin_idx = discretize(ocr, edges);
            cmap = parula(10);  % 10级渐变色，可根据需要更换其他colormap

            scatter(pos_tumor(:,1), pos_tumor(:,2), 10, cmap(bin_idx, :), 'filled');
            colormap(cmap);
            clim([edges(1) edges(end)]);
            if col_idx == 3
                colorbar;
            end
            axis equal;
            xlabel('x', 'FontSize', 8);
            ylabel('y', 'FontSize', 8);
            if row_idx == 1
                title(sprintf('$\\mu = %.1g$', mut_value), 'FontSize', 10, 'Interpreter', 'latex');
            end
        end

        % 仅第一列显示 y-axis 的所有信息
        if col_idx == 1
            ylabel(ax, sprintf('Time: t = %.1f', 14 + row_idx * 10), 'FontSize', 12);
        else
            set(ax, 'YTickLabel', []);
        end

        % 仅最后一行显示 x-axis 的所有信息
        if row_idx < 4
            set(ax, 'XTickLabel', []);
        else
            xlabel(ax, 'x', 'FontSize', 12);
        end

        % 如果你希望刻度是 0:0.2:1，则在这里设置：
        set(ax, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
    end
end

set(fig7, 'Color', 'w');
% Save the combined figure
savefolder = 'Figures/vascular/(epsilon = 1e-2)/mutation_spon';
% saveas(fig7, fullfile(savefolder, 'mutspon7.fig'), 'fig');
saveas(fig7, fullfile(savefolder, 'mutspon7.svg'), 'svg');

end
=======
function fig_mutspon()

clear all;
close all;

% Define the mutation rates we want to include
mut_rates = [1e-1, 1e-2, 1e-3, 1e-4];

%
%
% % Create a new figure for the combined plots
% fig1 = figure('Position', [100, 100, 1000, 800]);
% tlo1 = tiledlayout(4, 2, 'TileSpacing','tight','Padding','none');
% 
% % For each mutation rate
% for i = 1:length(mut_rates)
%     % Define folder path
%     folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/mutation_spon/mr_%.1g', mut_rates(i));
%     numericsFile = fullfile(folderName, 'numerics.mat');
% 
%     % Load the necessary data files
%     numerics_data = load(numericsFile);
% 
%     % Extract relevant data
%     numerics = numerics_data.numerics;
% 
%     % Get time vector
%     filename = 'Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat';
%     info = load(filename, 'info').info;
%     time_init = info.time;
%     max_time = size(numerics.tumor_num, 2);
%     time_vector = (time_init + (1:max_time)) / 10;
% 
%     % Left column: Plot tumor cell numbers
%     ax1 = nexttile(2 * i - 1);
% 
% % Plot tumor cell number with its components
% plot(time_vector, numerics.hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
% hold on;
% plot(time_vector, numerics.normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normoxic');
% plot(time_vector, numerics.tumor_num, 'b-', 'DisplayName', 'Tumor');
% hold off;
% 
% % Add mutation rate label at top-right in LaTeX
% text(0.75, 0.85, sprintf('$\\mu = 10^{-%d}$', i), ...
%     'Units', 'normalized', ...
%     'HorizontalAlignment', 'right', ...
%     'VerticalAlignment', 'top', ...
%     'FontSize', 10, ...
%     'Interpreter', 'latex');
% 
% if i == 4
%     xlabel('t', 'FontSize', 10);
% end
% ylabel('Cell Number', 'FontSize', 10);
% 
% % Set axis limits based on mutation rate
% xlim([time_init / 10, time_vector(end)]);
% 
% % Set tick labels
% xticks_vals = linspace(time_init / 10, time_vector(end), 6);
% xticks_labels = arrayfun(@(x) sprintf('%.2g', x), xticks_vals, 'UniformOutput', false);
% set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, 'FontSize', 8);
% 
% if i == 1
%     title(sprintf('Tumor Cell Number'), 'FontSize', 10);
% end
% if i == 4
%     legend('FontSize', 8, 'Location', 'best');
% end
% 
%     % Right column: Plot damage accumulation
%     ax2 = nexttile(2 * i);
% 
%     % Plot damage accumulation with error bars
%     plot(time_vector, numerics.dam_accum(1, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     hold on;
%     plot(time_vector, numerics.dam_accum(1, :) - numerics.dam_accum(2, :), 'b--', ...
%         'DisplayName', 'mean - std');
%     plot(time_vector, numerics.dam_accum(1, :) + numerics.dam_accum(2, :), 'b-.', ...
%         'DisplayName', 'mean + std');
% 
%     plot(time_vector, numerics.death_thres(1, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     plot(time_vector, numerics.death_thres(1, :) - numerics.death_thres(2, :), 'r--', ...
%         'DisplayName', 'mean - std');
%     plot(time_vector, numerics.death_thres(1, :) + numerics.death_thres(2, :), 'r-.', ...
%         'DisplayName', 'mean + std');
%     hold off;
% 
%     if i == 4
%         xlabel('t', 'FontSize', 10);
%     end
%     ylabel('Averaged Value', 'FontSize', 10);
% 
%     xlim([time_init / 10, time_vector(end)]);
%     % Set tick labels
%     xticks_vals = linspace(time_init / 10, time_vector(end), 6);
%     xticks_labels = arrayfun(@(x) sprintf('%.2g', x), xticks_vals, 'UniformOutput', false);
%     set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, 'FontSize', 8);
% 
%     if i == 1
%         title(sprintf('Average Damage Accumulation & \n Average Death Threshold'), 'FontSize', 10);
%     end
% 
%     if i == 4
%         legend('FontSize', 8, 'Location', 'northeast');
%     end
% 
%     % Add a text label for pr value (as a subfigure caption)
%     num_annotations = 4; % 或需要的总数
%     spacing = 0.95 / num_annotations;
%     annotation('textbox', [0, 0.945 - (i - 1) * spacing, 0.04, 0.04], ...
%         'String', sprintf('(%c)', 'a' + i - 1), ... % (a), (b), (c)
%         'FitBoxToText', 'on', ...
%         'EdgeColor', 'none', ...
%         'FontSize', 12);
% 
% end
% 
% % Set uniform background color
% set(fig1, 'Color', 'white');
% 
% % Save the combined figure
% savefolder = 'Figures/vascular/(epsilon = 1e-2)/mutation_spon';
% saveas(fig1, fullfile(savefolder, 'mutspon1.fig'), 'fig');
% saveas(fig1, fullfile(savefolder, 'mutspon1.svg'), 'svg');

%

% Create a new figure for the combined plots
fig2 = figure('Position', [100, 100, 1000, 900]);
tlo2 = tiledlayout(4, 2, 'TileSpacing','tight','Padding','none');

% For each mutation rate
for i = 1:length(mut_rates)
    % Define folder path
    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/mutation_spon/mr_%.1g', mut_rates(i));
    numericsFile = fullfile(folderName, 'numerics.mat');
    paramsFile = fullfile(folderName, 'params.mat');

    % Load the necessary data files
    numerics_data = load(numericsFile);
    params_data = load(paramsFile);

    % Extract relevant data
    numerics = numerics_data.numerics;
    params = params_data.params;

    % Get time vector
    filename = 'Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat';
    info = load(filename, 'info').info;
    time_init = info.time;
    max_time = size(numerics.tumor_num, 2);
    time_vector = (time_init + (1:max_time)) / 10;

    % Plot oxygen consumption histogram (left column)
    ax1 = nexttile(2 * i - 1);
    [X, Y] = meshgrid(linspace(0.5 * params.rhoo, 4 * params.rhoo, params.nbins), time_vector);
    s1 = surf(X, Y, numerics.oxygen_consumption_hist);

    % Improve appearance with better colormap and lighting
    colormap(gca, cool); % Use parula colormap (brighter than hot)
    s1.FaceAlpha = 0.8;    % Make surface more visible
    s1.EdgeColor = 'none'; % Remove grid lines for cleaner look

    xlabel('Oxygen Consumption', 'FontSize', 10);
    ylabel('Time', 'FontSize', 10);
    zlabel('Frequency', 'FontSize', 10);

    % Set X axis properties
    xlim([0.5 * params.rhoo, 4 * params.rhoo]);
    if i ~= 4
        ylim([time_init / 10, time_vector(end)]);
    else
        ylim([time_init / 10, 19.4]);
    end
    xticks_vals = 0.5 * params.rhoo : 0.35 * params.rhoo : 4 * params.rhoo;
    xticks_labels = arrayfun(@(x) sprintf('%.2f', x), xticks_vals, 'UniformOutput', false);
    if i ~= 4
        yticks_vals = linspace(time_init / 10, time_vector(end), 6);
        yticks_labels = arrayfun(@(x) sprintf('%.2g', x), yticks_vals, 'UniformOutput', false);
        set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, ...
            'YTick', yticks_vals, 'YTickLabel', yticks_labels, ...
            'FontSize', 8);
    else
        yticks_vals = linspace(time_init / 10, 19.4, 7);
        yticks_labels = arrayfun(@(x) sprintf('%.3g', x), yticks_vals, 'UniformOutput', false);
        set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, ...
            'YTick', yticks_vals, 'YTickLabel', yticks_labels, ...
            'FontSize', 8);
    end

    if i == 1
        title(sprintf('Oxygen Consumption'), 'FontSize', 12);
    end

    % Add mutation rate label at top-right in LaTeX
    text(0.75, 0.9, sprintf('$\\mu = 10^{-%d}$', i), ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 10, ...
    'Interpreter', 'latex');


    % Enhance lighting for better 3D appearance
    lighting phong;
    camlight('headlight');
    material shiny;
    view([-15, 50]);
    % colorbar;

    % Plot proliferation rate histogram (right column)
    ax2 = nexttile(2 * i);
    [X, Y] = meshgrid(linspace(0.5 * log(2) / params.age, 4 * log(2) / params.age, params.nbins), time_vector);
    s2 = surf(X, Y, numerics.proliferation_rate_hist);

    % Apply same visual improvements
    colormap(gca, cool);
    s2.FaceAlpha = 0.8;
    s2.EdgeColor = 'none';

    xlabel('Proliferation Rate', 'FontSize', 10);
    ylabel('Time', 'FontSize', 10);
    zlabel('Frequency', 'FontSize', 10);

    % Set X axis properties
    xlim([0.5 * log(2) / params.age, 4 * log(2) / params.age]);
    if i ~= 4
        ylim([time_init / 10, time_vector(end)]);
    else
        ylim([time_init / 10, 19.4]);
    end
    xticks_vals = 0.5 * log(2) / params.age : 0.35 * log(2) / params.age : 4 * log(2) / params.age;
    xticks_labels = arrayfun(@(x) sprintf('%.2f', x), xticks_vals, 'UniformOutput', false);
    if i ~= 4
        yticks_vals = linspace(time_init / 10, time_vector(end), 6);
        yticks_labels = arrayfun(@(x) sprintf('%.2g', x), yticks_vals, 'UniformOutput', false);
        set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, ...
            'YTick', yticks_vals, 'YTickLabel', yticks_labels, ...
            'FontSize', 8);
    else
        yticks_vals = linspace(time_init / 10, 19.4, 7);
        yticks_labels = arrayfun(@(x) sprintf('%.3g', x), yticks_vals, 'UniformOutput', false);
        set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, ...
            'YTick', yticks_vals, 'YTickLabel', yticks_labels, ...
            'FontSize', 8);
    end

    if i == 1
        title(sprintf('Proliferation Rate'), 'FontSize', 12);
    end

    lighting phong;
    camlight('headlight');
    material shiny;
    view([-15, 50]);
    % colorbar;

    % Add a text label for pr value (as a subfigure caption)
    num_annotations = 4; % 或需要的总数
    spacing = 0.95 / num_annotations;
    annotation('textbox', [0.01, 0.99 - i * spacing, 0.05, 0.05], ...
        'String', sprintf('(%c)', 'a' + i - 1), ... % (a), (b), (c)
        'FitBoxToText', 'on', ...
        'EdgeColor', 'none', ...
        'FontSize', 12);

end

% Set uniform background color
set(fig2, 'Color', 'white');

% Save the combined figure
savefolder = 'Figures/vascular/(epsilon = 1e-2)/mutation_spon';
saveas(fig2, fullfile(savefolder, 'mutspon2.fig'), 'fig');
saveas(fig2, fullfile(savefolder, 'mutspon2.svg'), 'svg');

%
% % 创建新图像
fig3 = figure('Position', [100, 100, 900, 1200]);
tlo3 = tiledlayout(4,3, 'TileSpacing','none','Padding','none');

% 遍历每个pr值
for col_idx = 1:3
    mut_value = mut_rates(col_idx);

    % 构建文件路径
    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/mutation_spon/mr_%.1g', mut_value);
    filePath = fullfile(folderName, 'vessel1.fig');

    orig_fig = openfig(filePath, 'invisible');
    all_axes = findobj(orig_fig, 'Type', 'axes');

    % 遍历每一行
    for row_idx = 1:4
        figure(fig3);
        new_subplot = nexttile(tlo3, (row_idx - 1) * 3 + col_idx);

        % 查找对应的原始子图
        found = false;
        for ax_idx = 1:length(all_axes)
            ax = all_axes(ax_idx);
            pos = get(ax, 'Position');

            if pos(1) < 0.3 && ((row_idx == 1 && pos(2) > 0.7) || ...
                    (row_idx == 2 && pos(2) > 0.45 && pos(2) <= 0.7) || ...
                    (row_idx == 3 && pos(2) > 0.2 && pos(2) <= 0.45) || ...
                    (row_idx == 4 && pos(2) <= 0.2))
                orig_subplot = ax;
                found = true;
                break;
            end
        end

        if found
            % 复制内容和属性
            children = get(orig_subplot, 'Children');
            for i = 1:length(children)
                copyobj(children(i), new_subplot);
            end

            % 改为（反转复制顺序）
            children = get(orig_subplot, 'Children');
            for i = length(children):-1:1
                copyobj(children(i), new_subplot);
            end

            % 复制轴属性
            set(new_subplot, 'XLim', get(orig_subplot, 'XLim'));
            set(new_subplot, 'YLim', get(orig_subplot, 'YLim'));
            set(new_subplot, 'XTick', get(orig_subplot, 'XTick'));
            set(new_subplot, 'YTick', get(orig_subplot, 'YTick'));
            set(new_subplot, 'FontSize', 10);

            % 使用简化方法创建图例 - 不依赖于原图的对象
            hold(new_subplot, 'on');

            % 完全脱离子图中的真实对象，创建虚拟图例项
            h1 = plot(new_subplot, NaN, NaN, 'bs', 'MarkerFaceColor', 'b', 'DisplayName', 'Vessel');
            h2 = plot(new_subplot, NaN, NaN, 'mo', 'MarkerFaceColor', 'm', 'DisplayName', 'Hypoxic');
            h3 = plot(new_subplot, NaN, NaN, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Normorxic');

            % 只使用需要的图例项
            leg_handles = [h1];
            leg_labels = {'Vessel'};

            % 根据原图中的图例字符串添加其他项
            orig_legend = legend(orig_subplot);
            if ~isempty(orig_legend)
                try
                    orig_legend_str = get(orig_legend, 'String');
                    if any(strcmp(orig_legend_str, 'Hypoxic'))
                        leg_handles = [leg_handles, h2];
                        leg_labels = [leg_labels, {'Hypoxic'}];
                    end
                    if any(strcmp(orig_legend_str, 'Normorxic'))
                        leg_handles = [leg_handles, h3];
                        leg_labels = [leg_labels, {'Normorxic'}];
                    end
                catch
                    % 如果无法获取原始图例，使用所有三个项
                    leg_handles = [h1, h2, h3];
                    leg_labels = {'Vessel', 'Hypoxic', 'Normorxic'};
                end
            end

            % 创建图例
            lgd = legend(new_subplot, leg_handles, leg_labels);
            set(lgd, 'FontSize', 6, 'Location', 'best');

            % 添加标签
            if col_idx == 1
                ylabel(new_subplot, sprintf('Time: t = %.1f', 14 + row_idx * 10), 'FontSize', 12);
            end
            if row_idx == 1
                title(sprintf('$\\mu = %.1g$', mut_value), 'FontSize', 10, 'Interpreter', 'latex')
            end
        end

        % 如果当前不是第一列，则去除 ytick 和 yticklabel
        if col_idx > 1
            set(new_subplot, 'YTickLabel', []);
        end

        % 如果当前不是最后一行，则去除 xtick 和 xticklabel
        if row_idx < 4
            set(new_subplot, 'XTickLabel', []);
        else
            xlabel(ax, 'x', 'FontSize', 10);
        end
    end
end

% 添加总标题和保存
sgtitle('Tumor Population Comparison for Different p_r Values', 'FontSize', 12);
set(fig3, 'Color', 'w');
% Save the combined figure
savefolder = 'Figures/vascular/(epsilon = 1e-2)/mutation_spon';
saveas(fig3, fullfile(savefolder, 'mutspon3.fig'), 'fig');
saveas(fig3, fullfile(savefolder, 'mutspon3.svg'), 'svg');

%% Figure4: 按照氧消耗率着色的肿瘤细胞
fig4 = figure('Position', [100, 100, 950, 1100]);
tlo4 = tiledlayout(4,3, 'TileSpacing','none','Padding','none');

for col_idx = 1:length(mut_rates) - 1
    mut_value = mut_rates(col_idx);

    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/mutation_spon/mr_%.1g', mut_value);
    filePath = fullfile(folderName, 'info.mat');
    info_data = load(filePath, 'info');
    info = info_data.info;
    params = info.params;

    for row_idx = 1:4
        ax = nexttile(tlo4, (row_idx - 1) * 3 + col_idx);

        tumor_cells = info.tumor_cells{row_idx};
        valid_tumor = ~cellfun('isempty', tumor_cells);

        if any(valid_tumor)
            % 提取细胞位置和氧消耗率（字段名称：oxygen_consumption_nor）
            pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
            pos_tumor = cat(1, pos_tumor{:});
            ocr = cellfun(@(tc) tc.oxygen_consumption_nor, tumor_cells(valid_tumor), 'UniformOutput', true);
            edges = linspace(0.5 * params.rhoo, 4 * params.rhoo, 11);
            bin_idx = discretize(ocr, edges);
            cmap = parula(10);  % 10级渐变色，可根据需要更换其他colormap

            scatter(pos_tumor(:,1), pos_tumor(:,2), 10, cmap(bin_idx, :), 'filled');
            colormap(cmap);
            clim([edges(1) edges(end)]);
            if col_idx == 3
                colorbar;
            end
            axis equal;
            xlabel('x', 'FontSize', 8);
            ylabel('y', 'FontSize', 8);
            if row_idx == 1
                title(sprintf('$\\mu = %.1g$', mut_value), 'FontSize', 10, 'Interpreter', 'latex');
            end
        end

        % 仅第一列显示 y-axis 的所有信息
        if col_idx == 1
            ylabel(ax, sprintf('Time: t = %.1f', 14 + row_idx * 10), 'FontSize', 12);
        else
            set(ax, 'YTickLabel', []);
        end

        % 仅最后一行显示 x-axis 的所有信息
        if row_idx < 4
            set(ax, 'XTickLabel', []);
        else
            xlabel(ax, 'x', 'FontSize', 12);
        end

        % 如果你希望刻度是 0:0.2:1，则在这里设置：
        set(ax, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
    end
end

set(fig4, 'Color', 'w');
% Save the combined figure
savefolder = 'Figures/vascular/(epsilon = 1e-2)/mutation_spon';
saveas(fig4, fullfile(savefolder, 'mutspon4.fig'), 'fig');
saveas(fig4, fullfile(savefolder, 'mutspon4.svg'), 'svg');

%% Figure5: 按照增殖率着色的肿瘤细胞
% 基准值 x = log(2)/params.age（未突变时的增殖率），控制范围为 [0.5*x, 4*x]，10等分
fig5 = figure('Position', [100, 100, 950, 1100]);
tlo5 = tiledlayout(4,3, 'TileSpacing','none','Padding','none');

for col_idx = 1:length(mut_rates) - 1
    mut_value = mut_rates(col_idx);

    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/mutation_spon/mr_%.1g', mut_value);
    filePath = fullfile(folderName, 'info.mat');
    info_data = load(filePath, 'info');
    info = info_data.info;
    params = info.params;

    for row_idx = 1:4
        ax = nexttile(tlo5, (row_idx - 1) * 3 + col_idx);

        tumor_cells = info.tumor_cells{row_idx};
        valid_tumor = ~cellfun('isempty', tumor_cells);

        if any(valid_tumor)
            pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
            pos_tumor = cat(1, pos_tumor{:});
            prate = cellfun(@(tc) tc.proliferation_rate, tumor_cells(valid_tumor));
            % 构造区间 [0.5*log(2)/params.age, 4*log(2)/params.age]并10等分
            edges = linspace(0.5 * log(2)/params.age, 4 * log(2)/params.age, 11);
            bin_idx = discretize(prate, edges);
            cmap = parula(10);

            scatter(pos_tumor(:,1), pos_tumor(:,2), 10, cmap(bin_idx, :), 'filled');
            colormap(cmap);
            clim([edges(1) edges(end)]);
            if col_idx == 3
                colorbar;
            end
            axis equal;
            xlabel('x', 'FontSize', 8);
            ylabel('y', 'FontSize', 8);
            if row_idx == 1
                title(sprintf('$\\mu = %.1g$', mut_value), 'FontSize', 10, 'Interpreter', 'latex');
            end
        end

        % 仅第一列显示 y-axis 的所有信息
        if col_idx == 1
            ylabel(ax, sprintf('Time: t = %.1f', 14 + row_idx * 10), 'FontSize', 12);
        else
            set(ax, 'YTickLabel', []);
        end

        % 仅最后一行显示 x-axis 的所有信息
        if row_idx < 4
            set(ax, 'XTickLabel', []);
        else
            xlabel(ax, 'x', 'FontSize', 12);
        end

        % 如果你希望刻度是 0:0.2:1，则在这里设置：
        set(ax, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
    end
end
set(fig5, 'Color', 'w');
% Save the combined figure
savefolder = 'Figures/vascular/(epsilon = 1e-2)/mutation_spon';
saveas(fig5, fullfile(savefolder, 'mutspon5.fig'), 'fig');
saveas(fig5, fullfile(savefolder, 'mutspon5.svg'), 'svg');

%
% 创建新图像
fig6 = figure('Position', [100, 100, 1250, 1150]);
tiledlayout(4, 4, 'TileSpacing','none','Padding','none');

% 遍历每个pr值
for col_idx = 1:length(mut_rates)
    mut_value = mut_rates(col_idx);

    % 构建文件路径
    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/mutation_spon/mr_%.1g', mut_value);
    info_data = load(fullfile(folderName, 'info.mat'));
    info = info_data.info;
    params = info.params;
    % filePath = fullfile(folderName, 'vessel1.fig');

    % orig_fig = openfig(filePath, 'invisible');
    % all_axes = findobj(orig_fig, 'Type', 'axes');

    % 遍历每一行
    for row_idx = 1:4
        figure(fig6);
        new_subplot = nexttile((row_idx - 1) * 4 + col_idx);

        drug = flipud(info.chem_field{row_idx}.drug);

        % Display the snapshot
        imagesc(drug);
         % Add colorbar to all columns but make non-last ones invisible
        if col_idx == 4
            colorbar;
        else
            cb = colorbar;
            cb.Visible = 'off'; % Invisible but reserves the space
        end
        axis tight equal;
        colormap(gca, 'parula');

        xlim([0, params.grid_size(2)]);
        ylim([0, params.grid_size(1)]);

        % Set axis labels and ticks
        if col_idx == 1
            ylabel(gca, sprintf('Time: t = %.1f', 14 + row_idx * 10), 'FontSize', 12);
        else
            set(new_subplot, 'YTickLabel', []);
        end

        if row_idx == 1
            title(sprintf('$\\mu = %.1g$', mut_value), 'FontSize', 10, 'Interpreter', 'latex')
        end

        if row_idx < 4
            set(new_subplot, 'XTickLabel', []);
        else
            xlabel('x', 'FontSize', 10);
        end

        if row_idx == 4
            set(gca, 'XTick', 0:params.grid_size(2) / 5:params.grid_size(2), 'XTickLabel', 0:0.2:1, ...
                'FontSize', 10);
        end
        if col_idx == 1
            set(gca, 'YTick', 0:params.grid_size(1) / 5:params.grid_size(1), 'YTickLabel', 1:-0.2:0, ...
                'FontSize', 10);
        end
    end

    % 如果当前不是第一列，则去除 ytick 和 yticklabel
    if col_idx > 1
        set(gca, 'YTickLabel', []);
    end
end
%
set(fig6, 'Color', 'w');
% % Save the combined figure
savefolder = 'Figures/vascular/(epsilon = 1e-2)/mutation_spon';
saveas(fig6, fullfile(savefolder, 'mutspon6.fig'), 'fig');
saveas(fig6, fullfile(savefolder, 'mutspon6.svg'), 'svg');


%% Figure7: 按照氧消耗率着色的肿瘤细胞
fig7 = figure('Position', [100, 100, 950, 1100]);
tiledlayout(4,3, 'TileSpacing','none','Padding','none');

for col_idx = 1:length(mut_rates) - 1
    mut_value = mut_rates(col_idx);

    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/mutation_spon/mr_%.1g', mut_value);
    filePath = fullfile(folderName, 'info.mat');
    info_data = load(filePath, 'info');
    info = info_data.info;
    % params = info.params;

    for row_idx = 1:4
        ax = nexttile((row_idx - 1) * 3 + col_idx);

        tumor_cells = info.tumor_cells{row_idx};
        valid_tumor = ~cellfun('isempty', tumor_cells);

        if any(valid_tumor)
            % 提取细胞位置和氧消耗率（字段名称：oxygen_consumption_nor）
            pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
            pos_tumor = cat(1, pos_tumor{:});
            ocr = cellfun(@(tc) tc.death_threshold, tumor_cells(valid_tumor));
            edges = linspace(min(ocr), max(ocr), 11);
            bin_idx = discretize(ocr, edges);
            cmap = parula(10);  % 10级渐变色，可根据需要更换其他colormap

            scatter(pos_tumor(:,1), pos_tumor(:,2), 10, cmap(bin_idx, :), 'filled');
            colormap(cmap);
            clim([edges(1) edges(end)]);
            if col_idx == 3
                colorbar;
            end
            axis equal;
            xlabel('x', 'FontSize', 8);
            ylabel('y', 'FontSize', 8);
            if row_idx == 1
                title(sprintf('$\\mu = %.1g$', mut_value), 'FontSize', 10, 'Interpreter', 'latex');
            end
        end

        % 仅第一列显示 y-axis 的所有信息
        if col_idx == 1
            ylabel(ax, sprintf('Time: t = %.1f', 14 + row_idx * 10), 'FontSize', 12);
        else
            set(ax, 'YTickLabel', []);
        end

        % 仅最后一行显示 x-axis 的所有信息
        if row_idx < 4
            set(ax, 'XTickLabel', []);
        else
            xlabel(ax, 'x', 'FontSize', 12);
        end

        % 如果你希望刻度是 0:0.2:1，则在这里设置：
        set(ax, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
    end
end

set(fig7, 'Color', 'w');
% Save the combined figure
savefolder = 'Figures/vascular/(epsilon = 1e-2)/mutation_spon';
% saveas(fig7, fullfile(savefolder, 'mutspon7.fig'), 'fig');
saveas(fig7, fullfile(savefolder, 'mutspon7.svg'), 'svg');

end
>>>>>>> d5569c1 (Initial commit)
