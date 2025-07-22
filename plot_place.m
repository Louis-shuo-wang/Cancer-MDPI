<<<<<<< HEAD
function plot_place()

clear all;
close all;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
time_points = zeros(1, div_num);

for i = 1:length(time_points)
    time_points(i) = i * plot_int * 0.1;
end

%%
for ijk = 1:3

    filename = 'Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat';
    info = load(filename, 'info').info;
    time_init = info.time;

    params = initializeParameters();
    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/place/place_%d', ijk);

    % Set which placement setting to test: 1, 2, or 3.
    resis_setting = ijk;   % Change this value to 1, 2, or 3

    if resis_setting == 1
        designated_positions = [0.5, 0.5];
    elseif resis_setting == 2
        designated_positions = [0.25, 0.25; 0.25, 0.75; 0.75, 0.25; 0.75, 0.75];
    elseif resis_setting == 3
        designated_positions = [0.25, 0.5; 0.75, 0.5; 0.5, 0.25; 0.5, 0.75];
    else
        error('Invalid resis_setting. Use 1, 2, or 3.');
    end

    % Get the pre-existing tumor cell indices from info.tumor_cells.
    tumor_cells = info.tumor_cells;

    % Find indices that are not empty.
    nonempty_idx = find(~cellfun('isempty', tumor_cells));

    % Total number of resistant cells: 1% of non-empty cells (at least 1).
    num_resis = max(1, ceil(params.resis_fraction * length(nonempty_idx)));

    % Decide how many cells to assign per designated region.
    num_regions = size(designated_positions, 1);
    num_per_region = ceil(num_resis / num_regions);

    % Tolerance for "closeness" (you may adjust this value)
    tol = 0.005;

    % Initialize an empty vector for selected indices.
    resis_idx = [];

    % Loop over each designated position and select cells that are sufficiently near.
    for r = 1:num_regions
        pos_target = designated_positions(r, :);
        % For every tumor cell in nonempty_idx, compute distance:
        distances = cellfun(@(tc) norm(tc.position - pos_target), tumor_cells(nonempty_idx));
        % Find indices within tolerance.
        idx_here = nonempty_idx(distances < tol);
        % If there are more available than needed, randomly sample num_per_region.
        while length(idx_here) < num_per_region
            tol = tol + 0.003;
            idx_here = nonempty_idx(distances < tol);
        end
        sel = randsample(idx_here, num_per_region, false);
        resis_idx = [resis_idx; sel(:)];
    end

    % If, after going through all regions, the total selected is different from num_resis,
    % then randomly sample (or pad) to make exactly num_resis.
    % if length(resis_idx) > num_resis
    %     resis_idx = randsample(resis_idx, num_resis, false);
    % elseif length(resis_idx) < num_resis
    %     % If insufficient cells were selected, add additional cells randomly from nonempty_idx not already chosen.
    %     available = setdiff(nonempty_idx, resis_idx);
    %     extra_needed = num_resis - length(resis_idx);
    %     if ~isempty(available)
    %         extra = randsample(available, min(extra_needed, length(available)), false);
    %         resis_idx = [resis_idx; extra(:)];
    %     end
    % end

    % Modify the selected tumor cells to have a higher death threshold (resistant cells)
    for i = 1:length(resis_idx)
        idx = resis_idx(i);
        tumor_cells{idx}.death_threshold = params.base_death_threshold * params.resistance_factor;
    end

    % Save the modified tumor cells back into info
    info.tumor_cells = tumor_cells;

    initsetting = struct();
    initsetting.info = info;
    initsetting.params = params;

    % main step
    [snapshots, numerics, info] = main(max_time, folderName, initsetting);
    hypoxic_num = numerics.hypoxic_num;
    normorxic_num = numerics.normorxic_num;
    tumor_num = numerics.tumor_num;
    dam_accum = numerics.dam_accum;
    death_thres = numerics.death_thres;
    params = info.params;

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
    fig2.Position = [100, 100, 1200, 400];
    time_vector = (time_init + (1:max_time)) / 10;

    subplot(1, 2, 1);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(1, 2, 2);
    hold on;
    % plot averaged damage accumulation with error bar
    errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
        'DisplayName', 'Damage Accumulation');
    errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
        'DisplayName', 'Death Threshold');
    xlabel('t', 'FontSize', 12);
    ylabel('Averaged Value', 'FontSize', 12);
    title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
    legend('Location', 'best', 'FontSize', 8);
    hold off;

    filePath = fullfile(folderName, 'num1');
    set(fig2, 'toolbar', 'none');
    saveas(fig2, [filePath, '.fig'], 'fig');
    % saveas(fig2, [filePath, '.eps'], 'eps');
    saveas(fig2, [filePath, '.svg'], 'svg');

    %%
    % Save the 'params' variable to a .mat file
    infoFilename = fullfile(folderName, 'info.mat');
    save(infoFilename, 'info');
    paramsFilename = fullfile(folderName, 'params.mat');
    save(paramsFilename, 'params');
    numericsFilename = fullfile(folderName, 'numerics.mat');
    save(numericsFilename, 'numerics');

end
=======
function plot_place()

clear all;
close all;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
time_points = zeros(1, div_num);

for i = 1:length(time_points)
    time_points(i) = i * plot_int * 0.1;
end

%%
for ijk = 1:3

    filename = 'Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat';
    info = load(filename, 'info').info;
    time_init = info.time;

    params = initializeParameters();
    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/place/place_%d', ijk);

    % Set which placement setting to test: 1, 2, or 3.
    resis_setting = ijk;   % Change this value to 1, 2, or 3

    if resis_setting == 1
        designated_positions = [0.5, 0.5];
    elseif resis_setting == 2
        designated_positions = [0.25, 0.25; 0.25, 0.75; 0.75, 0.25; 0.75, 0.75];
    elseif resis_setting == 3
        designated_positions = [0.25, 0.5; 0.75, 0.5; 0.5, 0.25; 0.5, 0.75];
    else
        error('Invalid resis_setting. Use 1, 2, or 3.');
    end

    % Get the pre-existing tumor cell indices from info.tumor_cells.
    tumor_cells = info.tumor_cells;

    % Find indices that are not empty.
    nonempty_idx = find(~cellfun('isempty', tumor_cells));

    % Total number of resistant cells: 1% of non-empty cells (at least 1).
    num_resis = max(1, ceil(params.resis_fraction * length(nonempty_idx)));

    % Decide how many cells to assign per designated region.
    num_regions = size(designated_positions, 1);
    num_per_region = ceil(num_resis / num_regions);

    % Tolerance for "closeness" (you may adjust this value)
    tol = 0.005;

    % Initialize an empty vector for selected indices.
    resis_idx = [];

    % Loop over each designated position and select cells that are sufficiently near.
    for r = 1:num_regions
        pos_target = designated_positions(r, :);
        % For every tumor cell in nonempty_idx, compute distance:
        distances = cellfun(@(tc) norm(tc.position - pos_target), tumor_cells(nonempty_idx));
        % Find indices within tolerance.
        idx_here = nonempty_idx(distances < tol);
        % If there are more available than needed, randomly sample num_per_region.
        while length(idx_here) < num_per_region
            tol = tol + 0.003;
            idx_here = nonempty_idx(distances < tol);
        end
        sel = randsample(idx_here, num_per_region, false);
        resis_idx = [resis_idx; sel(:)];
    end

    % If, after going through all regions, the total selected is different from num_resis,
    % then randomly sample (or pad) to make exactly num_resis.
    % if length(resis_idx) > num_resis
    %     resis_idx = randsample(resis_idx, num_resis, false);
    % elseif length(resis_idx) < num_resis
    %     % If insufficient cells were selected, add additional cells randomly from nonempty_idx not already chosen.
    %     available = setdiff(nonempty_idx, resis_idx);
    %     extra_needed = num_resis - length(resis_idx);
    %     if ~isempty(available)
    %         extra = randsample(available, min(extra_needed, length(available)), false);
    %         resis_idx = [resis_idx; extra(:)];
    %     end
    % end

    % Modify the selected tumor cells to have a higher death threshold (resistant cells)
    for i = 1:length(resis_idx)
        idx = resis_idx(i);
        tumor_cells{idx}.death_threshold = params.base_death_threshold * params.resistance_factor;
    end

    % Save the modified tumor cells back into info
    info.tumor_cells = tumor_cells;

    initsetting = struct();
    initsetting.info = info;
    initsetting.params = params;

    % main step
    [snapshots, numerics, info] = main(max_time, folderName, initsetting);
    hypoxic_num = numerics.hypoxic_num;
    normorxic_num = numerics.normorxic_num;
    tumor_num = numerics.tumor_num;
    dam_accum = numerics.dam_accum;
    death_thres = numerics.death_thres;
    params = info.params;

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
    fig2.Position = [100, 100, 1200, 400];
    time_vector = (time_init + (1:max_time)) / 10;

    subplot(1, 2, 1);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(1, 2, 2);
    hold on;
    % plot averaged damage accumulation with error bar
    errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
        'DisplayName', 'Damage Accumulation');
    errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
        'DisplayName', 'Death Threshold');
    xlabel('t', 'FontSize', 12);
    ylabel('Averaged Value', 'FontSize', 12);
    title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
    legend('Location', 'best', 'FontSize', 8);
    hold off;

    filePath = fullfile(folderName, 'num1');
    set(fig2, 'toolbar', 'none');
    saveas(fig2, [filePath, '.fig'], 'fig');
    % saveas(fig2, [filePath, '.eps'], 'eps');
    saveas(fig2, [filePath, '.svg'], 'svg');

    %%
    % Save the 'params' variable to a .mat file
    infoFilename = fullfile(folderName, 'info.mat');
    save(infoFilename, 'info');
    paramsFilename = fullfile(folderName, 'params.mat');
    save(paramsFilename, 'params');
    numericsFilename = fullfile(folderName, 'numerics.mat');
    save(numericsFilename, 'numerics');

end
>>>>>>> d5569c1 (Initial commit)
end