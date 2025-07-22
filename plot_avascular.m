<<<<<<< HEAD
% This file is used to plot tumor development
% tumor shift from avascular to vascular stage
function plot_avascular()

clear all;
close all;

max_time = 4;
div_num = 4;
plot_int = max_time / div_num;
time_points = zeros(1, div_num);

for i = 1:length(time_points)
    time_points(i) = i * plot_int * 0.1;
end

% filepath to save data
folderName = 'Figures/2025_Jul/vascular0/grid_refinement';
time_init = 0;
save_path = folderName;

params = initializeParameters();
% This is where we make modifications
% in this simulation, we apply no drug 
params.Sd = 0;

initsetting  = struct();
initsetting.info = [];
initsetting.params = params;

[snapshots, numerics, info] = main(max_time, folderName, [], initsetting);
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;
vessel_onset_time = numerics.vessel_onset_time;
params = info.params;
tumor_cells = info.tumor_cells;
chem_field = info.chem_field;
vessel_agents = info.vessel_agents;
Nbins = params.nbins;
%%
% Create the main figure
fig1 = figure(1);
fig1.Visible = 'off';
fig1.Position = [100, 100, 1200, 1200];

% Create a 4\time 4 grid of subplots
for k = 1:4
    subplot(4, 4, (k - 1) * 4 + 1);
    plotTumorAndVessels(snapshots{1}{k}.tumor_cells, ...
        snapshots{1}{k}.vessel_agents, ...
        snapshots{1}{k}.params);
    if (snapshots{1}{k}.time - time_init) / plot_int == 1
        title('vessel', 'FontSize', 12);
    end
    % Set axis labels and ticks
    axis equal;
    xlabel('x', 'FontSize', 12);
    % Set the row title
    ylabel(sprintf('Time: t = %.1f', snapshots{1}{k}.time * 0.1), 'FontSize', 12);
    xlim([0, 1]);
    ylim([0, 1]);
    set(gca, 'XTick', 0:0.2:1, ...
        'YTick', 0:0.2:1, ...
        'FontSize', 10);
end

for i = 1:4
    for j = 2:4         % Loop over variables vessel, TAF, oxygen, drug (columns)
        subplot_idx = (i - 1) * 4 + j;
        subplot(4, 4, subplot_idx);

        % Get the corresponding snapshot
        snapshot = snapshots{j}{i};

        % Display the snapshot
        imagesc(snapshot);
        colorbar;
        axis tight equal;
        colormap(gca, 'parula');

        % Set the title
        if i == 1
            title(snapshots{5}{j}, 'FontSize', 12);
        end

        % Set axis labels and ticks
        xlabel('x', 'FontSize', 12);
        ylabel('y', 'FontSize', 12);
        xlim([0, params.grid_size(2)]);
        ylim([0, params.grid_size(1)]);
        set(gca, 'XTick', 0:params.grid_size(2) / 5:params.grid_size(2), 'XTickLabel', 0:0.2:1, ...
            'YTick', 0:params.grid_size(1) / 5:params.grid_size(1), 'YTickLabel', 1:-0.2:0, ...
            'FontSize', 10);
    end
end

filePath = fullfile(folderName, 'vessel1');
set(fig1, 'toolbar', 'none');
saveas(fig1, [filePath, '.fig'], 'fig');
% saveas(fig1, [filePath, '.eps'], 'eps');
saveas(fig1, [filePath, '.svg'], 'svg');

%%
% This plot is for tracking the cell number in our simulations
fig2 = figure(2);
fig2.Visible = 'off';
fig2.Position = [100, 100, 1200, 800];
time_vector = (time_init + (1:max_time)) / 10;

% subplot(1, 2, 1);
plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
hold on;
plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
hold off;

% Add vertical line for vessel onset
if ~isempty(vessel_onset_time) && ~isnan(vessel_onset_time)
    xline(vessel_onset_time * params.dt, 'm--', 'LineWidth', 2, ...
        'Label', 'Angiogenesis', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
end

legend('FontSize', 8, 'Location', 'best');
xlabel('t', 'FontSize', 12);
ylabel('Cell Number', 'FontSize', 12);
title('Tumor Cell Number', 'FontSize', 12);

filePath = fullfile(folderName, 'num1');
set(fig2, 'toolbar', 'none');
saveas(fig2, [filePath, '.fig'], 'fig');
saveas(fig2, [filePath, '.svg'], 'svg');

% runTumorModelAnalysis(numerics, tumor_cells, vessel_agents, chem_field, params, save_path, Nbins)
%%
% Save the 'params' variable to a .mat file
infoFilename = fullfile(folderName, 'info.mat');
save(infoFilename, 'info');
paramsFilename = fullfile(folderName, 'params.mat');
save(paramsFilename, 'params');
numericsFilename = fullfile(folderName, 'numerics.mat');
save(numericsFilename, 'numerics');

=======
% This file is used to plot tumor development
% tumor shift from avascular to vascular stage
function plot_avascular()

clear all;
close all;

max_time = 4;
div_num = 4;
plot_int = max_time / div_num;
time_points = zeros(1, div_num);

for i = 1:length(time_points)
    time_points(i) = i * plot_int * 0.1;
end

% filepath to save data
folderName = 'Figures/2025_Jul/vascular0/grid_refinement';
time_init = 0;
save_path = folderName;

params = initializeParameters();
% This is where we make modifications
% in this simulation, we apply no drug 
params.Sd = 0;

initsetting  = struct();
initsetting.info = [];
initsetting.params = params;

[snapshots, numerics, info] = main(max_time, folderName, [], initsetting);
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;
vessel_onset_time = numerics.vessel_onset_time;
params = info.params;
tumor_cells = info.tumor_cells;
chem_field = info.chem_field;
vessel_agents = info.vessel_agents;
Nbins = params.nbins;
%%
% Create the main figure
fig1 = figure(1);
fig1.Visible = 'off';
fig1.Position = [100, 100, 1200, 1200];

% Create a 4\time 4 grid of subplots
for k = 1:4
    subplot(4, 4, (k - 1) * 4 + 1);
    plotTumorAndVessels(snapshots{1}{k}.tumor_cells, ...
        snapshots{1}{k}.vessel_agents, ...
        snapshots{1}{k}.params);
    if (snapshots{1}{k}.time - time_init) / plot_int == 1
        title('vessel', 'FontSize', 12);
    end
    % Set axis labels and ticks
    axis equal;
    xlabel('x', 'FontSize', 12);
    % Set the row title
    ylabel(sprintf('Time: t = %.1f', snapshots{1}{k}.time * 0.1), 'FontSize', 12);
    xlim([0, 1]);
    ylim([0, 1]);
    set(gca, 'XTick', 0:0.2:1, ...
        'YTick', 0:0.2:1, ...
        'FontSize', 10);
end

for i = 1:4
    for j = 2:4         % Loop over variables vessel, TAF, oxygen, drug (columns)
        subplot_idx = (i - 1) * 4 + j;
        subplot(4, 4, subplot_idx);

        % Get the corresponding snapshot
        snapshot = snapshots{j}{i};

        % Display the snapshot
        imagesc(snapshot);
        colorbar;
        axis tight equal;
        colormap(gca, 'parula');

        % Set the title
        if i == 1
            title(snapshots{5}{j}, 'FontSize', 12);
        end

        % Set axis labels and ticks
        xlabel('x', 'FontSize', 12);
        ylabel('y', 'FontSize', 12);
        xlim([0, params.grid_size(2)]);
        ylim([0, params.grid_size(1)]);
        set(gca, 'XTick', 0:params.grid_size(2) / 5:params.grid_size(2), 'XTickLabel', 0:0.2:1, ...
            'YTick', 0:params.grid_size(1) / 5:params.grid_size(1), 'YTickLabel', 1:-0.2:0, ...
            'FontSize', 10);
    end
end

filePath = fullfile(folderName, 'vessel1');
set(fig1, 'toolbar', 'none');
saveas(fig1, [filePath, '.fig'], 'fig');
% saveas(fig1, [filePath, '.eps'], 'eps');
saveas(fig1, [filePath, '.svg'], 'svg');

%%
% This plot is for tracking the cell number in our simulations
fig2 = figure(2);
fig2.Visible = 'off';
fig2.Position = [100, 100, 1200, 800];
time_vector = (time_init + (1:max_time)) / 10;

% subplot(1, 2, 1);
plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
hold on;
plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
hold off;

% Add vertical line for vessel onset
if ~isempty(vessel_onset_time) && ~isnan(vessel_onset_time)
    xline(vessel_onset_time * params.dt, 'm--', 'LineWidth', 2, ...
        'Label', 'Angiogenesis', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
end

legend('FontSize', 8, 'Location', 'best');
xlabel('t', 'FontSize', 12);
ylabel('Cell Number', 'FontSize', 12);
title('Tumor Cell Number', 'FontSize', 12);

filePath = fullfile(folderName, 'num1');
set(fig2, 'toolbar', 'none');
saveas(fig2, [filePath, '.fig'], 'fig');
saveas(fig2, [filePath, '.svg'], 'svg');

% runTumorModelAnalysis(numerics, tumor_cells, vessel_agents, chem_field, params, save_path, Nbins)
%%
% Save the 'params' variable to a .mat file
infoFilename = fullfile(folderName, 'info.mat');
save(infoFilename, 'info');
paramsFilename = fullfile(folderName, 'params.mat');
save(paramsFilename, 'params');
numericsFilename = fullfile(folderName, 'numerics.mat');
save(numericsFilename, 'numerics');

>>>>>>> d5569c1 (Initial commit)
end