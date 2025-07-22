<<<<<<< HEAD
function plot_drug_mutspon()

clear all;
close all;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
time_points = zeros(1, div_num);
treatment = [10, 10; ...
             20, 5; ...
             30, 10/3; ...
             40, 5/2; ...
             50, 5; ...
             50, 10];

for i = 1:length(time_points)
    time_points(i) = i * plot_int * 0.1;
end

mut_rt = 1e-3;
%%
for ijk = 1:size(treatment, 1)
    filename = 'Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat';
    info = load(filename, 'info').info;
    time_init = info.time;

    params = initializeParameters();
    params.mutation_rate = mut_rt;
    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/drug_mutspon(1e-3)/treat_%d', ijk);

    initsetting = struct();
    initsetting.info = info;
    initsetting.params = params;

    % main step
    [snapshots, numerics, info] = main(max_time, folderName, treatment(ijk, :), initsetting);
    hypoxic_num = numerics.hypoxic_num;
    normorxic_num = numerics.normorxic_num;
    tumor_num = numerics.tumor_num;
    dam_accum = numerics.dam_accum;
    death_thres = numerics.death_thres;
    oxygen_consumption_hist = numerics.oxygen_consumption_hist;
    proliferation_rate_hist = numerics.proliferation_rate_hist;
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
    saveas(fig1, [filePath, '.svg'], 'svg');

    %%
    % This plot is for tracking the cell number in our simulations
    fig2 = figure(2);
    fig2.Visible = 'off';
    fig2.Position = [100, 100, 800, 400];
    time_vector = (time_init + (1:max_time)) / 10;

    subplot(2, 2, 1);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 6, 'Location', 'best');
    xlabel('t', 'FontSize', 10);
    ylabel('Cell Number', 'FontSize', 10);
    title('Tumor Cell Number', 'FontSize', 10);

    subplot(2, 2, 2);
    hold on;
    % plot averaged damage accumulation with error bar
    errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
        'DisplayName', 'Damage Accumulation');
    errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
        'DisplayName', 'Death Threshold');
    xlabel('t', 'FontSize', 10);
    ylabel('Averaged Value', 'FontSize', 10);
    title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 10);
    legend('Location', 'best', 'FontSize', 6);
    hold off;

    subplot(2, 2, 3);
    [X, Y] = meshgrid(linspace(0.5 * params.rhoo, 4 * params.rhoo, params.nbins), time_vector);
    surf(X, Y, oxygen_consumption_hist, 'FaceAlpha', 0.6);
    colormap(gca, hot);
    xlabel('oxygen consumption', 'FontSize', 10);
    ylabel('t', 'FontSize', 10);
    zlabel('frequency', 'FontSize', 10);
    xlim([0.5 * params.rhoo, 4 * params.rhoo]);
    xticks_vals = 0.5 * params.rhoo : 0.35 * params.rhoo : 4 * params.rhoo;
    xticks_labels = arrayfun(@(x) sprintf('%.2f', x), xticks_vals, 'UniformOutput', false);
    set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, 'FontSize', 8);
    title('Histogram Surface of Oxygen Consumption Rate', 'FontSize', 10);

    subplot(2, 2, 4);
    [X, Y] = meshgrid(linspace(0.5 * log(2) / params.age, 4 * log(2) / params.age, params.nbins), ...
        time_vector);
    surf(X, Y, proliferation_rate_hist, 'FaceAlpha', 0.6);
    colormap(gca, hot);
    xlabel('proliferation rate', 'FontSize', 10);
    ylabel('t', 'FontSize', 10);
    zlabel('frequency', 'FontSize', 10);
    xlim([0.5 * log(2) / params.age, 4 * log(2) / params.age]);
    xticks_vals = 0.5 * log(2) / params.age : 0.35 * log(2) / params.age : 4 * log(2) / params.age;
    xticks_labels = arrayfun(@(x) sprintf('%.2f', x), xticks_vals, 'UniformOutput', false);
    set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, 'FontSize', 8);
    title('Histogram Surface of Proliferation Rate', 'FontSize', 10);

    filePath = fullfile(folderName, 'num1');
    set(fig2, 'toolbar', 'none');
    saveas(fig2, [filePath, '.fig'], 'fig');
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
function plot_drug_mutspon()

clear all;
close all;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
time_points = zeros(1, div_num);
treatment = [10, 10; ...
             20, 5; ...
             30, 10/3; ...
             40, 5/2; ...
             50, 5; ...
             50, 10];

for i = 1:length(time_points)
    time_points(i) = i * plot_int * 0.1;
end

mut_rt = 1e-3;
%%
for ijk = 1:size(treatment, 1)
    filename = 'Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat';
    info = load(filename, 'info').info;
    time_init = info.time;

    params = initializeParameters();
    params.mutation_rate = mut_rt;
    folderName = sprintf('Figures/vascular/(epsilon = 1e-2)/drug_mutspon(1e-3)/treat_%d', ijk);

    initsetting = struct();
    initsetting.info = info;
    initsetting.params = params;

    % main step
    [snapshots, numerics, info] = main(max_time, folderName, treatment(ijk, :), initsetting);
    hypoxic_num = numerics.hypoxic_num;
    normorxic_num = numerics.normorxic_num;
    tumor_num = numerics.tumor_num;
    dam_accum = numerics.dam_accum;
    death_thres = numerics.death_thres;
    oxygen_consumption_hist = numerics.oxygen_consumption_hist;
    proliferation_rate_hist = numerics.proliferation_rate_hist;
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
    saveas(fig1, [filePath, '.svg'], 'svg');

    %%
    % This plot is for tracking the cell number in our simulations
    fig2 = figure(2);
    fig2.Visible = 'off';
    fig2.Position = [100, 100, 800, 400];
    time_vector = (time_init + (1:max_time)) / 10;

    subplot(2, 2, 1);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 6, 'Location', 'best');
    xlabel('t', 'FontSize', 10);
    ylabel('Cell Number', 'FontSize', 10);
    title('Tumor Cell Number', 'FontSize', 10);

    subplot(2, 2, 2);
    hold on;
    % plot averaged damage accumulation with error bar
    errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
        'DisplayName', 'Damage Accumulation');
    errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
        'DisplayName', 'Death Threshold');
    xlabel('t', 'FontSize', 10);
    ylabel('Averaged Value', 'FontSize', 10);
    title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 10);
    legend('Location', 'best', 'FontSize', 6);
    hold off;

    subplot(2, 2, 3);
    [X, Y] = meshgrid(linspace(0.5 * params.rhoo, 4 * params.rhoo, params.nbins), time_vector);
    surf(X, Y, oxygen_consumption_hist, 'FaceAlpha', 0.6);
    colormap(gca, hot);
    xlabel('oxygen consumption', 'FontSize', 10);
    ylabel('t', 'FontSize', 10);
    zlabel('frequency', 'FontSize', 10);
    xlim([0.5 * params.rhoo, 4 * params.rhoo]);
    xticks_vals = 0.5 * params.rhoo : 0.35 * params.rhoo : 4 * params.rhoo;
    xticks_labels = arrayfun(@(x) sprintf('%.2f', x), xticks_vals, 'UniformOutput', false);
    set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, 'FontSize', 8);
    title('Histogram Surface of Oxygen Consumption Rate', 'FontSize', 10);

    subplot(2, 2, 4);
    [X, Y] = meshgrid(linspace(0.5 * log(2) / params.age, 4 * log(2) / params.age, params.nbins), ...
        time_vector);
    surf(X, Y, proliferation_rate_hist, 'FaceAlpha', 0.6);
    colormap(gca, hot);
    xlabel('proliferation rate', 'FontSize', 10);
    ylabel('t', 'FontSize', 10);
    zlabel('frequency', 'FontSize', 10);
    xlim([0.5 * log(2) / params.age, 4 * log(2) / params.age]);
    xticks_vals = 0.5 * log(2) / params.age : 0.35 * log(2) / params.age : 4 * log(2) / params.age;
    xticks_labels = arrayfun(@(x) sprintf('%.2f', x), xticks_vals, 'UniformOutput', false);
    set(gca, 'XTick', xticks_vals, 'XTickLabel', xticks_labels, 'FontSize', 8);
    title('Histogram Surface of Proliferation Rate', 'FontSize', 10);

    filePath = fullfile(folderName, 'num1');
    set(fig2, 'toolbar', 'none');
    saveas(fig2, [filePath, '.fig'], 'fig');
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