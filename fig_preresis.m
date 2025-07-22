<<<<<<< HEAD
function fig_preresis()

clear all;
close all;

% Define the pr values to extract and combine
pr_values = [0.2, 0.3, 1];

% Create a new figure with the size accommodating 3 rows
fig1 = figure('Position', [100, 100, 900, 600]);
tlo1 = tiledlayout(3, 2, 'TileSpacing','tight','Padding','none');

% Loop through each pr value
for i = 1:length(pr_values)
    % Define folder path
    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/vascular2/pr_%.1g', pr_values(i));
    numericsFile = fullfile(folderName, 'numerics.mat');

    % Load the necessary data files
    numerics_data = load(numericsFile);

    % Extract relevant data
    numerics = numerics_data.numerics;

    % Get time vector
    filename = 'Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat';
    info = load(filename, 'info').info;
    time_init = info.time;
    max_time = size(numerics.tumor_num, 2);
    time_vector = (time_init + (1:max_time)) / 10;

    % Left column: Plot tumor cell numbers
    nexttile(2 * i - 1);

    % Plot tumor cell number with its components
    plot(time_vector, numerics.hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, numerics.normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normoxic');
    plot(time_vector, numerics.tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;

    if i == 3
        xlabel('t', 'FontSize', 10);
    end

    text(0.55, 0.89, sprintf('$p_r = %.1f$', pr_values(i)), ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 10, ...
    'Interpreter', 'latex');

    ylabel('Cell Number', 'FontSize', 10);

    xlim([14, 34]);
    xticks_vals = 14:4:34;
    xticks_labels = arrayfun(@(x) sprintf('%.2g', x), xticks_vals, 'UniformOutput', false);

    set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, 'FontSize', 8);

    if i == 1
        title(sprintf('Tumor Cell Number'), 'FontSize', 10);
    end
    if i == length(pr_values)
        legend('FontSize', 8, 'Location', 'best');
    end

    % Right column: Plot damage accumulation
    nexttile(2 * i);

    % Plot damage accumulation with error bars
    errorbar(time_vector, numerics.dam_accum(1, :), numerics.dam_accum(2, :), 'b-', ...
        'DisplayName', 'Damage Accumulation');
    hold on;
    errorbar(time_vector, numerics.death_thres(1, :), numerics.death_thres(2, :), 'r-', ...
        'DisplayName', 'Death Threshold');
    hold off;

    if i == 3
        xlabel('t', 'FontSize', 10);
    end
    ylabel('Averaged Value', 'FontSize', 10);

    xlim([14, 34]);
    % Set tick labels
    xticks_vals = 14:4:34;
    xticks_labels = arrayfun(@(x) sprintf('%.2g', x), xticks_vals, 'UniformOutput', false);
    set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, 'FontSize', 8);

    if i == 1
        title(sprintf('Average Damage Accumulation & \n Average Death Threshold'), 'FontSize', 10);
    end

    if i == length(pr_values)
        legend('FontSize', 8, 'Location', 'best');
    end

    % Add a text label for pr value (as a subfigure caption)
    num_annotations = length(pr_values); % 或需要的总数
    spacing = 0.95 / num_annotations;
    annotation('textbox', [0, 0.95 - (i - 1) * spacing, 0.04, 0.04], ...
        'String', sprintf('(%c)', 'a' + i - 1), ... % (a), (b), (c)
        'FitBoxToText', 'on', ...
        'EdgeColor', 'none', ...
        'FontSize', 12);


end

% % Manually adjust the position of all subplots to avoid overlap
% all_axes = findall(fig1, 'Type', 'axes');
%
% % Adjust left column subplots (odd indices)
% for i = 1:2:length(all_axes)
%     pos = get(all_axes(i), 'Position');
%     % Reduce width and shift right to avoid being covered
%     set(all_axes(i), 'Position', [pos(1), pos(2), pos(3) * 0.88, pos(4)]);
% end
%
% % Adjust right column subplots (even indices)
% for i = 2:2:length(all_axes)
%     pos = get(all_axes(i), 'Position');
%     % Shift right to avoid covering left plots
%     set(all_axes(i), 'Position', [pos(1) * 1.1, pos(2), pos(3) * 0.88, pos(4)]);
% end

% % Add caption for the entire figure
% annotation('textbox', [0.25, 0.01, 0.5, 0.03], ...
%           'String', sprintf('Figure: p_{r} = 0, 0.3, 1 in (a), (b), and (c), respectively.'), ...
%           'FitBoxToText', 'on', ...
%           'EdgeColor', 'none', ...
%           'HorizontalAlignment', 'center', ...
%           'FontSize', 12);

% Save the combined figure
savefolder = 'Figures/vascular/(epsilon = 1e-2)/vascular2';
saveas(fig1, fullfile(savefolder, 'preresis1.fig'), 'fig');
saveas(fig1, fullfile(savefolder, 'preresis1.svg'), 'svg');


%%
% % 创建新图像
% fig2 = figure('Position', [100, 100, 900, 1200]);
% tlo2 = tiledlayout(4, 3, 'TileSpacing','none','Padding','none');
%
% % 遍历每个pr值
% for col_idx = 1:length(pr_values)
%     pr_value = pr_values(col_idx);
%
%     % 构建文件路径
%     folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/vascular2/pr_%.1g', pr_value);
%     filePath = fullfile(folderName, 'vessel1.fig');
%
%     % if col_idx == 4
%     %     folderName = 'Figures/vascular/(epsilon = 1e-2)/vascular2/pr_0';
%     %     filePath = fullfile(folderName, 'vessel1.fig');
%     % end
%
%     fprintf('正在处理文件: %s\n', filePath);
%
%     orig_fig = openfig(filePath, 'invisible');
%     all_axes = findobj(orig_fig, 'Type', 'axes');
%
%     % 遍历每一行
%     for row_idx = 1:4
%         figure(fig2);
%         new_subplot = nexttile(tlo2, (row_idx - 1) * 3 + col_idx);
%
%         % 查找对应的原始子图
%         found = false;
%         for ax_idx = 1:length(all_axes)
%             ax = all_axes(ax_idx);
%             pos = get(ax, 'Position');
%
%             if pos(1) < 0.3 && ((row_idx == 1 && pos(2) > 0.7) || ...
%                     (row_idx == 2 && pos(2) > 0.45 && pos(2) <= 0.7) || ...
%                     (row_idx == 3 && pos(2) > 0.2 && pos(2) <= 0.45) || ...
%                     (row_idx == 4 && pos(2) <= 0.2))
%                 orig_subplot = ax;
%                 found = true;
%                 break;
%             end
%         end
%
%         if found
%             % 复制内容和属性
%             children = get(orig_subplot, 'Children');
%             for i = 1:length(children)
%                 copyobj(children(i), new_subplot);
%             end
%
%             % 改为（反转复制顺序）
%             children = get(orig_subplot, 'Children');
%             for i = length(children):-1:1
%                 copyobj(children(i), new_subplot);
%             end
%
%             % 复制轴属性
%             set(new_subplot, 'XLim', get(orig_subplot, 'XLim'));
%             set(new_subplot, 'YLim', get(orig_subplot, 'YLim'));
%             set(new_subplot, 'XTick', get(orig_subplot, 'XTick'));
%             set(new_subplot, 'YTick', get(orig_subplot, 'YTick'));
%             set(new_subplot, 'FontSize', 10);
%
%             % 使用简化方法创建图例 - 不依赖于原图的对象
%             hold(new_subplot, 'on');
%
%             % 完全脱离子图中的真实对象，创建虚拟图例项
%             h1 = plot(new_subplot, NaN, NaN, 'bs', 'MarkerFaceColor', 'b', 'DisplayName', 'Vessel');
%             h2 = plot(new_subplot, NaN, NaN, 'mo', 'MarkerFaceColor', 'm', 'DisplayName', 'Hypoxic');
%             h3 = plot(new_subplot, NaN, NaN, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Normorxic');
%
%             % 只使用需要的图例项
%             leg_handles = [h1];
%             leg_labels = {'Vessel'};
%
%             % 根据原图中的图例字符串添加其他项
%             orig_legend = legend(orig_subplot);
%             if ~isempty(orig_legend)
%                 try
%                     orig_legend_str = get(orig_legend, 'String');
%                     if any(strcmp(orig_legend_str, 'Hypoxic'))
%                         leg_handles = [leg_handles, h2];
%                         leg_labels = [leg_labels, {'Hypoxic'}];
%                     end
%                     if any(strcmp(orig_legend_str, 'Normorxic'))
%                         leg_handles = [leg_handles, h3];
%                         leg_labels = [leg_labels, {'Normorxic'}];
%                     end
%                 catch
%                     % 如果无法获取原始图例，使用所有三个项
%                     leg_handles = [h1, h2, h3];
%                     leg_labels = {'Vessel', 'Hypoxic', 'Normorxic'};
%                 end
%             end
%
%             % 创建图例
%             lgd = legend(new_subplot, leg_handles, leg_labels);
%             set(lgd, 'FontSize', 6, 'Location', 'best');
%
%             % 添加标签
%             if col_idx == 1
%                 ylabel(new_subplot, sprintf('Time: t = %.1f', 14 + row_idx * 4), 'FontSize', 12);
%             else
%                 set(new_subplot, 'YTickLabel', []);
%             end
%             if row_idx == 1
%                 title(new_subplot, sprintf('p_r = %.1g', pr_value));
%             end
%
%             % 如果当前不是最后一行，则去除 xtick 和 xticklabel
%             if row_idx < 4
%                 set(new_subplot, 'XTickLabel', []);
%             else
%                 xlabel(ax, 'x', 'FontSize', 10);
%             end
%         end
%     end
% end
%
% % 添加总标题和保存
% % sgtitle('Tumor Population Comparison for Different p_r Values', 'FontSize', 12);
% set(fig2, 'Color', 'w');
% savefolder = 'Figures/vascular/(epsilon = 1e-2)/vascular2';
% saveas(fig2, fullfile(savefolder, 'preresis2.fig'), 'fig');
% saveas(fig2, fullfile(savefolder, 'preresis2.svg'), 'svg');
%
% fprintf('Combined figure created and saved in %s\n', savefolder);
=======
function fig_preresis()

clear all;
close all;

% Define the pr values to extract and combine
pr_values = [0.2, 0.3, 1];

% Create a new figure with the size accommodating 3 rows
fig1 = figure('Position', [100, 100, 900, 600]);
tlo1 = tiledlayout(3, 2, 'TileSpacing','tight','Padding','none');

% Loop through each pr value
for i = 1:length(pr_values)
    % Define folder path
    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/vascular2/pr_%.1g', pr_values(i));
    numericsFile = fullfile(folderName, 'numerics.mat');

    % Load the necessary data files
    numerics_data = load(numericsFile);

    % Extract relevant data
    numerics = numerics_data.numerics;

    % Get time vector
    filename = 'Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat';
    info = load(filename, 'info').info;
    time_init = info.time;
    max_time = size(numerics.tumor_num, 2);
    time_vector = (time_init + (1:max_time)) / 10;

    % Left column: Plot tumor cell numbers
    nexttile(2 * i - 1);

    % Plot tumor cell number with its components
    plot(time_vector, numerics.hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, numerics.normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normoxic');
    plot(time_vector, numerics.tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;

    if i == 3
        xlabel('t', 'FontSize', 10);
    end

    text(0.55, 0.89, sprintf('$p_r = %.1f$', pr_values(i)), ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 10, ...
    'Interpreter', 'latex');

    ylabel('Cell Number', 'FontSize', 10);

    xlim([14, 34]);
    xticks_vals = 14:4:34;
    xticks_labels = arrayfun(@(x) sprintf('%.2g', x), xticks_vals, 'UniformOutput', false);

    set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, 'FontSize', 8);

    if i == 1
        title(sprintf('Tumor Cell Number'), 'FontSize', 10);
    end
    if i == length(pr_values)
        legend('FontSize', 8, 'Location', 'best');
    end

    % Right column: Plot damage accumulation
    nexttile(2 * i);

    % Plot damage accumulation with error bars
    errorbar(time_vector, numerics.dam_accum(1, :), numerics.dam_accum(2, :), 'b-', ...
        'DisplayName', 'Damage Accumulation');
    hold on;
    errorbar(time_vector, numerics.death_thres(1, :), numerics.death_thres(2, :), 'r-', ...
        'DisplayName', 'Death Threshold');
    hold off;

    if i == 3
        xlabel('t', 'FontSize', 10);
    end
    ylabel('Averaged Value', 'FontSize', 10);

    xlim([14, 34]);
    % Set tick labels
    xticks_vals = 14:4:34;
    xticks_labels = arrayfun(@(x) sprintf('%.2g', x), xticks_vals, 'UniformOutput', false);
    set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, 'FontSize', 8);

    if i == 1
        title(sprintf('Average Damage Accumulation & \n Average Death Threshold'), 'FontSize', 10);
    end

    if i == length(pr_values)
        legend('FontSize', 8, 'Location', 'best');
    end

    % Add a text label for pr value (as a subfigure caption)
    num_annotations = length(pr_values); % 或需要的总数
    spacing = 0.95 / num_annotations;
    annotation('textbox', [0, 0.95 - (i - 1) * spacing, 0.04, 0.04], ...
        'String', sprintf('(%c)', 'a' + i - 1), ... % (a), (b), (c)
        'FitBoxToText', 'on', ...
        'EdgeColor', 'none', ...
        'FontSize', 12);


end

% % Manually adjust the position of all subplots to avoid overlap
% all_axes = findall(fig1, 'Type', 'axes');
%
% % Adjust left column subplots (odd indices)
% for i = 1:2:length(all_axes)
%     pos = get(all_axes(i), 'Position');
%     % Reduce width and shift right to avoid being covered
%     set(all_axes(i), 'Position', [pos(1), pos(2), pos(3) * 0.88, pos(4)]);
% end
%
% % Adjust right column subplots (even indices)
% for i = 2:2:length(all_axes)
%     pos = get(all_axes(i), 'Position');
%     % Shift right to avoid covering left plots
%     set(all_axes(i), 'Position', [pos(1) * 1.1, pos(2), pos(3) * 0.88, pos(4)]);
% end

% % Add caption for the entire figure
% annotation('textbox', [0.25, 0.01, 0.5, 0.03], ...
%           'String', sprintf('Figure: p_{r} = 0, 0.3, 1 in (a), (b), and (c), respectively.'), ...
%           'FitBoxToText', 'on', ...
%           'EdgeColor', 'none', ...
%           'HorizontalAlignment', 'center', ...
%           'FontSize', 12);

% Save the combined figure
savefolder = 'Figures/vascular/(epsilon = 1e-2)/vascular2';
saveas(fig1, fullfile(savefolder, 'preresis1.fig'), 'fig');
saveas(fig1, fullfile(savefolder, 'preresis1.svg'), 'svg');


%%
% % 创建新图像
% fig2 = figure('Position', [100, 100, 900, 1200]);
% tlo2 = tiledlayout(4, 3, 'TileSpacing','none','Padding','none');
%
% % 遍历每个pr值
% for col_idx = 1:length(pr_values)
%     pr_value = pr_values(col_idx);
%
%     % 构建文件路径
%     folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/vascular2/pr_%.1g', pr_value);
%     filePath = fullfile(folderName, 'vessel1.fig');
%
%     % if col_idx == 4
%     %     folderName = 'Figures/vascular/(epsilon = 1e-2)/vascular2/pr_0';
%     %     filePath = fullfile(folderName, 'vessel1.fig');
%     % end
%
%     fprintf('正在处理文件: %s\n', filePath);
%
%     orig_fig = openfig(filePath, 'invisible');
%     all_axes = findobj(orig_fig, 'Type', 'axes');
%
%     % 遍历每一行
%     for row_idx = 1:4
%         figure(fig2);
%         new_subplot = nexttile(tlo2, (row_idx - 1) * 3 + col_idx);
%
%         % 查找对应的原始子图
%         found = false;
%         for ax_idx = 1:length(all_axes)
%             ax = all_axes(ax_idx);
%             pos = get(ax, 'Position');
%
%             if pos(1) < 0.3 && ((row_idx == 1 && pos(2) > 0.7) || ...
%                     (row_idx == 2 && pos(2) > 0.45 && pos(2) <= 0.7) || ...
%                     (row_idx == 3 && pos(2) > 0.2 && pos(2) <= 0.45) || ...
%                     (row_idx == 4 && pos(2) <= 0.2))
%                 orig_subplot = ax;
%                 found = true;
%                 break;
%             end
%         end
%
%         if found
%             % 复制内容和属性
%             children = get(orig_subplot, 'Children');
%             for i = 1:length(children)
%                 copyobj(children(i), new_subplot);
%             end
%
%             % 改为（反转复制顺序）
%             children = get(orig_subplot, 'Children');
%             for i = length(children):-1:1
%                 copyobj(children(i), new_subplot);
%             end
%
%             % 复制轴属性
%             set(new_subplot, 'XLim', get(orig_subplot, 'XLim'));
%             set(new_subplot, 'YLim', get(orig_subplot, 'YLim'));
%             set(new_subplot, 'XTick', get(orig_subplot, 'XTick'));
%             set(new_subplot, 'YTick', get(orig_subplot, 'YTick'));
%             set(new_subplot, 'FontSize', 10);
%
%             % 使用简化方法创建图例 - 不依赖于原图的对象
%             hold(new_subplot, 'on');
%
%             % 完全脱离子图中的真实对象，创建虚拟图例项
%             h1 = plot(new_subplot, NaN, NaN, 'bs', 'MarkerFaceColor', 'b', 'DisplayName', 'Vessel');
%             h2 = plot(new_subplot, NaN, NaN, 'mo', 'MarkerFaceColor', 'm', 'DisplayName', 'Hypoxic');
%             h3 = plot(new_subplot, NaN, NaN, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Normorxic');
%
%             % 只使用需要的图例项
%             leg_handles = [h1];
%             leg_labels = {'Vessel'};
%
%             % 根据原图中的图例字符串添加其他项
%             orig_legend = legend(orig_subplot);
%             if ~isempty(orig_legend)
%                 try
%                     orig_legend_str = get(orig_legend, 'String');
%                     if any(strcmp(orig_legend_str, 'Hypoxic'))
%                         leg_handles = [leg_handles, h2];
%                         leg_labels = [leg_labels, {'Hypoxic'}];
%                     end
%                     if any(strcmp(orig_legend_str, 'Normorxic'))
%                         leg_handles = [leg_handles, h3];
%                         leg_labels = [leg_labels, {'Normorxic'}];
%                     end
%                 catch
%                     % 如果无法获取原始图例，使用所有三个项
%                     leg_handles = [h1, h2, h3];
%                     leg_labels = {'Vessel', 'Hypoxic', 'Normorxic'};
%                 end
%             end
%
%             % 创建图例
%             lgd = legend(new_subplot, leg_handles, leg_labels);
%             set(lgd, 'FontSize', 6, 'Location', 'best');
%
%             % 添加标签
%             if col_idx == 1
%                 ylabel(new_subplot, sprintf('Time: t = %.1f', 14 + row_idx * 4), 'FontSize', 12);
%             else
%                 set(new_subplot, 'YTickLabel', []);
%             end
%             if row_idx == 1
%                 title(new_subplot, sprintf('p_r = %.1g', pr_value));
%             end
%
%             % 如果当前不是最后一行，则去除 xtick 和 xticklabel
%             if row_idx < 4
%                 set(new_subplot, 'XTickLabel', []);
%             else
%                 xlabel(ax, 'x', 'FontSize', 10);
%             end
%         end
%     end
% end
%
% % 添加总标题和保存
% % sgtitle('Tumor Population Comparison for Different p_r Values', 'FontSize', 12);
% set(fig2, 'Color', 'w');
% savefolder = 'Figures/vascular/(epsilon = 1e-2)/vascular2';
% saveas(fig2, fullfile(savefolder, 'preresis2.fig'), 'fig');
% saveas(fig2, fullfile(savefolder, 'preresis2.svg'), 'svg');
%
% fprintf('Combined figure created and saved in %s\n', savefolder);
>>>>>>> d5569c1 (Initial commit)
end