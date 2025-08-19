function fig_notumor()

clear all;
close all;

% 打开 vessel1.fig
vessel_filepath = 'Figures/vascular/notumor/vessel1.fig';
vessel_fig = openfig(vessel_filepath, 'visible');
drawnow;

% 找到所有axes，筛选第一列子图 (x在[0.12,0.14])
vessel_axes = findall(vessel_fig, 'Type', 'axes');
first_column_axes = [];
for i = 1:length(vessel_axes)
    pos = get(vessel_axes(i), 'Position');
    if pos(1) >= 0.12 && pos(1) <= 0.14
        first_column_axes = [first_column_axes, vessel_axes(i)];
    end
end

% 按y轴位置降序排列（从上到下）
pos_vals = arrayfun(@(ax) get(ax, 'Position'), first_column_axes, 'UniformOutput', false);
pos_y = cellfun(@(p) p(2), pos_vals);
[~, idx] = sort(pos_y, 'descend');
first_column_axes = first_column_axes(idx);

% 新建2x2大图
fig = figure('Position', [100, 100, 600, 600]);
set(fig, 'PaperPositionMode', 'auto');

% 新图中4个位置，手动设定，保证2x2布局
% 每个子图大小
w = 0.33; h = 0.33;
new_positions = [
    0.1 0.57 w h;
    0.55 0.57 w h;
    0.1 0.12  w h;
    0.55 0.12  w h;
];
% 按你的映射顺序，原第一列4个对应新图子图索引：
% 原1 → 新1
% 原2 → 新2
% 原3 → 新3
% 原4 → 新4

% 标签和字体大小
fontsize = 14;
title_fontsize = 18;

for i = 1:4
    src_ax = first_column_axes(i);

    ax_new = axes('Position', new_positions(i,:));
    
    % 复制子图内容
    children = allchild(src_ax);
    for j = 1:length(children)
        if isa(children(j), 'matlab.graphics.chart.primitive.Line')
            new_line = copyobj(children(j), ax_new);
            set(new_line, 'LineWidth', max(get(children(j), 'LineWidth')*2, 1.5)); % 加粗线条
        else
            copyobj(children(j), ax_new);
        end
    end
    
    % 复制轴属性
    set(ax_new, ...
        'XLim', get(src_ax, 'XLim'), ...
        'YLim', get(src_ax, 'YLim'), ...
        'XTick', get(src_ax, 'XTick'), ...
        'YTick', get(src_ax, 'YTick'), ...
        'Box', get(src_ax, 'Box'), ...
        'XGrid', get(src_ax, 'XGrid'), ...
        'YGrid', get(src_ax, 'YGrid'), ...
        'FontSize', fontsize, ...
        'FontWeight', 'bold', ...
        'LineWidth', 1.5, ...
        'TickLength', get(src_ax, 'TickLength')*1.5);
    
    % 加大xlabel和ylabel字体
    xlab = get(get(src_ax, 'XLabel'), 'String');
    if ~isempty(xlab)
        xlabel(ax_new, xlab, 'FontSize', fontsize, 'FontWeight', 'bold');
    end
    ylab = get(get(src_ax, 'YLabel'), 'String');
    if ~isempty(ylab)
        ylabel(ax_new, ylab, 'FontSize', fontsize, 'FontWeight', 'bold');
    end

    % 第一行两个子图上方添加大标题“vessel”
    if i == 1 || i == 2
        title(ax_new, 'vessel', 'FontSize', title_fontsize, 'FontWeight', 'bold');
    end
    
end

% 关闭原fig，避免占用
close(vessel_fig);

% 保存新图
savefolder = 'Figures/vascular/notumor';
if ~exist(savefolder, 'dir')
    mkdir(savefolder);
end

saveas(fig, fullfile(savefolder, 'vessel_tip_2x2_figure.fig'), 'fig');
saveas(fig, fullfile(savefolder, 'vessel_tip_2x2_figure.svg'), 'svg');

end