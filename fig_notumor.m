<<<<<<< HEAD
function fig_notumor()

clear all;
close all;

% 创建新图，尺寸设为[1000,1200]
fig = figure('Position', [100, 100, 1000, 1200]);
set(fig, 'PaperPositionMode', 'auto');

%% 左侧处理：加载 vessel1.fig 的第一列子图
vessel_filepath = 'Figures/vascular/notumor/vessel1.fig';
vessel_fig = openfig(vessel_filepath, 'visible');
drawnow;

% 获取所有 axes 对象
vessel_axes = findall(vessel_fig, 'Type', 'axes');
fprintf('Found %d axes in vessel1.fig\n', length(vessel_axes));

% 筛选 x 坐标在 [0.12, 0.14] 内的 axes（左侧子图）
first_column_axes = [];
for i = 1:length(vessel_axes)
    pos = get(vessel_axes(i), 'Position');
    if pos(1) >= 0.12 && pos(1) <= 0.14
        first_column_axes = [first_column_axes, vessel_axes(i)];
        fprintf('Left col: Added axes %d with pos: [%.3f, %.3f, %.3f, %.3f]\n', i, pos);
    end
end

% 按 y 方向降序排序（从上到下）
pos_vals = cellfun(@(ax)get(ax, 'Position'), num2cell(first_column_axes), 'UniformOutput', false);
pos_vals = cellfun(@(p)p(2), pos_vals);
[~, idx] = sort(pos_vals, 'descend');
first_column_axes = first_column_axes(idx);

% 设置新的 ylabel 文本
left_labels = {'t = 5.76', 't = 11.52', 't = 17.28', 't = 23.04'};

% 完整复制左侧各子图的axes属性
for i = 1:min(4, length(first_column_axes))
    src_ax = first_column_axes(i);

    % 在新图中创建子图，位置与原图完全相同
    figure(fig);
    ax_new = axes('Position', get(src_ax, 'Position'));

    % 复制所有子对象（线条、文本等）
    copyobj(allchild(src_ax), ax_new);

    % 复制坐标轴范围和刻度
    set(ax_new, 'XLim', get(src_ax, 'XLim'));
    set(ax_new, 'YLim', get(src_ax, 'YLim'));
    set(ax_new, 'XTick', get(src_ax, 'XTick'));
    set(ax_new, 'YTick', get(src_ax, 'YTick'));
    set(ax_new, 'XTickLabel', get(src_ax, 'XTickLabel'));
    set(ax_new, 'YTickLabel', get(src_ax, 'YTickLabel'));
    set(ax_new, 'Box', get(src_ax, 'Box'));
    set(ax_new, 'XGrid', get(src_ax, 'XGrid'));
    set(ax_new, 'YGrid', get(src_ax, 'YGrid'));

    % 保留原始的 xlabel，但自定义 ylabel
    xlab = get(get(src_ax, 'XLabel'), 'String');
    if ~isempty(xlab)
        xlabel(ax_new, xlab, 'FontSize', 10);
    end
    ylabel(ax_new, left_labels{i}, 'FontSize', 10);

    % 仅第一个子图保留原始的 title
    if i == 1
        ttl = get(get(src_ax, 'Title'), 'String');
        if ~isempty(ttl)
            title(ax_new, ttl, 'FontSize', 12);
        end
    end
end

%% 右侧处理：使用直接绘制方法而不是图像捕获
num_filepath = 'Figures/vascular/notumor/num2.fig';
num_fig = openfig(num_filepath, 'visible');
drawnow;

% 获取所有 axes 对象
num_axes = findall(num_fig, 'Type', 'axes');
fprintf('Found %d axes in num2.fig\n', length(num_axes));

% 按 x 坐标升序排序（从左到右）
pos_vals_num = zeros(length(num_axes), 1);
for i = 1:length(num_axes)
    pos = get(num_axes(i), 'Position');
    pos_vals_num(i) = pos(1);
    fprintf('num2.fig axes %d: pos = [%.3f %.3f %.3f %.3f]\n', i, pos);
end
[~, idx_num] = sort(pos_vals_num, 'ascend');
num_axes = num_axes(idx_num);

% 定义新位置 - 确保右侧有足够空间且不会与左侧重叠
right_x = 0.5;
right_width = 0.4;
top_right_bottom = 0.548;
top_right_height = 0.374;
bottom_right_bottom = 0.110;
bottom_right_height = 0.374;

subplot_pos = {
    [right_x, top_right_bottom, right_width, top_right_height],
    [right_x, bottom_right_bottom, right_width, bottom_right_height]
    };

% 处理两个右侧子图 - 通过复制数据而不是捕获图像
for i = 1:min(2, length(num_axes))
    src_ax = num_axes(i);
    pos = subplot_pos{i};

    % 在新图中创建子图
    figure(fig);
    ax_new = axes('Position', pos);

    % 提取并复制源坐标轴上的所有线条对象
    lines = findobj(src_ax, 'Type', 'line');
    for j = 1:length(lines)
        % 获取线条数据
        x = get(lines(j), 'XData');
        y = get(lines(j), 'YData');
        linestyle = get(lines(j), 'LineStyle');
        linewidth = get(lines(j), 'LineWidth');
        color = get(lines(j), 'Color');

        % 在新坐标轴上重新绘制线条
        plot(ax_new, x, y, 'LineStyle', linestyle, 'LineWidth', linewidth, 'Color', color);
        hold(ax_new, 'on');
    end

    % 复制坐标轴设置
    xlim(ax_new, get(src_ax, 'XLim'));
    ylim(ax_new, get(src_ax, 'YLim'));

    % 添加标签和标题
    xlabel(ax_new, 'Position', 'FontSize', 10);
    if i == 1
        ylabel(ax_new, 'Horizontal Tip Number', 'FontSize', 10);
    else
        ylabel(ax_new, 'Horizontal Vessel Number', 'FontSize', 10);
    end

    % 只为第二个图添加新图例
    if i == 2
        % 删除原来的图例
        delete(findobj(ax_new, 'Type', 'legend'));

        % 创建新的自定义图例线条（蓝色，按要求的顺序和样式）
        h1 = plot(NaN, NaN, 'k-', 'LineWidth', 1.5, 'DisplayName', sprintf('t = 5.76'));
        h2 = plot(NaN, NaN, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('t = 11.52'));
        h3 = plot(NaN, NaN, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('t = 17.28'));
        h4 = plot(NaN, NaN, 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('t = 23.04'));

        % 创建图例，确保按正确顺序（从上到下）排列句柄
        legend([h1, h2, h3, h4], 'FontSize', 16, 'Location', 'northwest');

        % % 调整图例框属性来防止文本溢出
        % lgd.Box = 'on';                   % 显示图例框
        % % lgd.Interpreter = 'none';         % 防止特殊字符解释
        % % lgd.ItemTokenSize = [8, 18];      % 减小标记尺寸
        % lgd.Position(3) = lgd.Position(3) * 1.2; % 增加宽度
    end

    hold(ax_new, 'off');
end

% 关闭原始图形
close(vessel_fig);
close(num_fig);

% 刷新显示
drawnow;

% 保存合成图到指定目录
savefolder = 'Figures/vascular/notumor';
if ~exist(savefolder, 'dir')
    mkdir(savefolder);
end

% 保存为各种格式
saveas(fig, fullfile(savefolder, 'combined_vessel_tip_figure.fig'), 'fig');
saveas(fig, fullfile(savefolder, 'combined_vessel_tip_figure.svg'), 'svg');

fprintf('Combined figure saved successfully.\n');
=======
function fig_notumor()

clear all;
close all;

% 创建新图，尺寸设为[1000,1200]
fig = figure('Position', [100, 100, 1000, 1200]);
set(fig, 'PaperPositionMode', 'auto');

%% 左侧处理：加载 vessel1.fig 的第一列子图
vessel_filepath = 'Figures/vascular/notumor/vessel1.fig';
vessel_fig = openfig(vessel_filepath, 'visible');
drawnow;

% 获取所有 axes 对象
vessel_axes = findall(vessel_fig, 'Type', 'axes');
fprintf('Found %d axes in vessel1.fig\n', length(vessel_axes));

% 筛选 x 坐标在 [0.12, 0.14] 内的 axes（左侧子图）
first_column_axes = [];
for i = 1:length(vessel_axes)
    pos = get(vessel_axes(i), 'Position');
    if pos(1) >= 0.12 && pos(1) <= 0.14
        first_column_axes = [first_column_axes, vessel_axes(i)];
        fprintf('Left col: Added axes %d with pos: [%.3f, %.3f, %.3f, %.3f]\n', i, pos);
    end
end

% 按 y 方向降序排序（从上到下）
pos_vals = cellfun(@(ax)get(ax, 'Position'), num2cell(first_column_axes), 'UniformOutput', false);
pos_vals = cellfun(@(p)p(2), pos_vals);
[~, idx] = sort(pos_vals, 'descend');
first_column_axes = first_column_axes(idx);

% 设置新的 ylabel 文本
left_labels = {'t = 5.76', 't = 11.52', 't = 17.28', 't = 23.04'};

% 完整复制左侧各子图的axes属性
for i = 1:min(4, length(first_column_axes))
    src_ax = first_column_axes(i);

    % 在新图中创建子图，位置与原图完全相同
    figure(fig);
    ax_new = axes('Position', get(src_ax, 'Position'));

    % 复制所有子对象（线条、文本等）
    copyobj(allchild(src_ax), ax_new);

    % 复制坐标轴范围和刻度
    set(ax_new, 'XLim', get(src_ax, 'XLim'));
    set(ax_new, 'YLim', get(src_ax, 'YLim'));
    set(ax_new, 'XTick', get(src_ax, 'XTick'));
    set(ax_new, 'YTick', get(src_ax, 'YTick'));
    set(ax_new, 'XTickLabel', get(src_ax, 'XTickLabel'));
    set(ax_new, 'YTickLabel', get(src_ax, 'YTickLabel'));
    set(ax_new, 'Box', get(src_ax, 'Box'));
    set(ax_new, 'XGrid', get(src_ax, 'XGrid'));
    set(ax_new, 'YGrid', get(src_ax, 'YGrid'));

    % 保留原始的 xlabel，但自定义 ylabel
    xlab = get(get(src_ax, 'XLabel'), 'String');
    if ~isempty(xlab)
        xlabel(ax_new, xlab, 'FontSize', 10);
    end
    ylabel(ax_new, left_labels{i}, 'FontSize', 10);

    % 仅第一个子图保留原始的 title
    if i == 1
        ttl = get(get(src_ax, 'Title'), 'String');
        if ~isempty(ttl)
            title(ax_new, ttl, 'FontSize', 12);
        end
    end
end

%% 右侧处理：使用直接绘制方法而不是图像捕获
num_filepath = 'Figures/vascular/notumor/num2.fig';
num_fig = openfig(num_filepath, 'visible');
drawnow;

% 获取所有 axes 对象
num_axes = findall(num_fig, 'Type', 'axes');
fprintf('Found %d axes in num2.fig\n', length(num_axes));

% 按 x 坐标升序排序（从左到右）
pos_vals_num = zeros(length(num_axes), 1);
for i = 1:length(num_axes)
    pos = get(num_axes(i), 'Position');
    pos_vals_num(i) = pos(1);
    fprintf('num2.fig axes %d: pos = [%.3f %.3f %.3f %.3f]\n', i, pos);
end
[~, idx_num] = sort(pos_vals_num, 'ascend');
num_axes = num_axes(idx_num);

% 定义新位置 - 确保右侧有足够空间且不会与左侧重叠
right_x = 0.5;
right_width = 0.4;
top_right_bottom = 0.548;
top_right_height = 0.374;
bottom_right_bottom = 0.110;
bottom_right_height = 0.374;

subplot_pos = {
    [right_x, top_right_bottom, right_width, top_right_height],
    [right_x, bottom_right_bottom, right_width, bottom_right_height]
    };

% 处理两个右侧子图 - 通过复制数据而不是捕获图像
for i = 1:min(2, length(num_axes))
    src_ax = num_axes(i);
    pos = subplot_pos{i};

    % 在新图中创建子图
    figure(fig);
    ax_new = axes('Position', pos);

    % 提取并复制源坐标轴上的所有线条对象
    lines = findobj(src_ax, 'Type', 'line');
    for j = 1:length(lines)
        % 获取线条数据
        x = get(lines(j), 'XData');
        y = get(lines(j), 'YData');
        linestyle = get(lines(j), 'LineStyle');
        linewidth = get(lines(j), 'LineWidth');
        color = get(lines(j), 'Color');

        % 在新坐标轴上重新绘制线条
        plot(ax_new, x, y, 'LineStyle', linestyle, 'LineWidth', linewidth, 'Color', color);
        hold(ax_new, 'on');
    end

    % 复制坐标轴设置
    xlim(ax_new, get(src_ax, 'XLim'));
    ylim(ax_new, get(src_ax, 'YLim'));

    % 添加标签和标题
    xlabel(ax_new, 'Position', 'FontSize', 10);
    if i == 1
        ylabel(ax_new, 'Horizontal Tip Number', 'FontSize', 10);
    else
        ylabel(ax_new, 'Horizontal Vessel Number', 'FontSize', 10);
    end

    % 只为第二个图添加新图例
    if i == 2
        % 删除原来的图例
        delete(findobj(ax_new, 'Type', 'legend'));

        % 创建新的自定义图例线条（蓝色，按要求的顺序和样式）
        h1 = plot(NaN, NaN, 'k-', 'LineWidth', 1.5, 'DisplayName', sprintf('t = 5.76'));
        h2 = plot(NaN, NaN, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('t = 11.52'));
        h3 = plot(NaN, NaN, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('t = 17.28'));
        h4 = plot(NaN, NaN, 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('t = 23.04'));

        % 创建图例，确保按正确顺序（从上到下）排列句柄
        legend([h1, h2, h3, h4], 'FontSize', 16, 'Location', 'northwest');

        % % 调整图例框属性来防止文本溢出
        % lgd.Box = 'on';                   % 显示图例框
        % % lgd.Interpreter = 'none';         % 防止特殊字符解释
        % % lgd.ItemTokenSize = [8, 18];      % 减小标记尺寸
        % lgd.Position(3) = lgd.Position(3) * 1.2; % 增加宽度
    end

    hold(ax_new, 'off');
end

% 关闭原始图形
close(vessel_fig);
close(num_fig);

% 刷新显示
drawnow;

% 保存合成图到指定目录
savefolder = 'Figures/vascular/notumor';
if ~exist(savefolder, 'dir')
    mkdir(savefolder);
end

% 保存为各种格式
saveas(fig, fullfile(savefolder, 'combined_vessel_tip_figure.fig'), 'fig');
saveas(fig, fullfile(savefolder, 'combined_vessel_tip_figure.svg'), 'svg');

fprintf('Combined figure saved successfully.\n');
>>>>>>> d5569c1 (Initial commit)
end