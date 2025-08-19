function plot_parameter_sensitivity()
    clear all;
    close all;

mutation_rates = [1e-4, 1e-3, 1e-2];
prs = [0.1, 0.2, 0.3];
Sds = [1, 2, 4];
Dds = [0.25, 0.5, 1];
Sos = [2.8, 3.5, 4.2];

    
    % Get time vector (assuming it's the same for all cases)
    filename = 'Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat';
    info = load(filename, 'info').info;
    time_init = info.time;

%%
folderName = sprintf('parameter_sensitivity/mutation_rate/mutation_rate_%.1g', mutation_rates(1));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig1 = figure(1);
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

folderName = sprintf('parameter_sensitivity/mutation_rate/mutation_rate_%.1g', mutation_rates(2));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig2 = figure(2);
    fig2.Position = [100, 100, 1200, 1200];

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
    
folderName = sprintf('parameter_sensitivity/mutation_rate/mutation_rate_%.1g', mutation_rates(3));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig3 = figure(3);
    fig3.Position = [100, 100, 1200, 1200];

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
 
folderName = sprintf('parameter_sensitivity/mutation_rate/mutation_rate_%.1g', mutation_rates(1));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;

max_time = size(numerics.tumor_num, 2);
time_vector = (time_init + (1:max_time)) / 10;
fig4 = figure(4);

    fig4.Position = [100, 100, 1200, 400];
    title('\mu', 'Interpreter', 'latex');
    subplot(3, 2, 1);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 2);
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


    folderName = sprintf('parameter_sensitivity/mutation_rate/mutation_rate_%.1g', mutation_rates(2));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;

    subplot(3, 2, 3);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 4);
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

    folderName = sprintf('parameter_sensitivity/mutation_rate/mutation_rate_%.1g', mutation_rates(3));
    % folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(3));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;


    subplot(3, 2, 5);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 6);
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

    %%
    folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(1));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig5 = figure(5);
    fig5.Position = [100, 100, 1200, 1200];

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

folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(2));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig6 = figure(6);
    fig6.Position = [100, 100, 1200, 1200];

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
    
folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(3));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig7 = figure(7);
    fig7.Position = [100, 100, 1200, 1200];

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
 
folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(1));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;

max_time = size(numerics.tumor_num, 2);
time_vector = (time_init + (1:max_time)) / 10;
fig8 = figure(8);

    fig8.Position = [100, 100, 1200, 400];
    title('\mu', 'Interpreter', 'latex');
    subplot(3, 2, 1);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 2);
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


    folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(2));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;

    subplot(3, 2, 3);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 4);
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

    folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(3));
    % folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(3));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;


    subplot(3, 2, 5);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 6);
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

    %%
    folderName = sprintf('parameter_sensitivity/Sd/Sd_%.1g', Sds(1));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig9 = figure(9);
    fig9.Position = [100, 100, 1200, 1200];

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

folderName = sprintf('parameter_sensitivity/Sd/Sd_%.1g', Sds(2));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig10 = figure(10);
    fig10.Position = [100, 100, 1200, 1200];

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
    
folderName = sprintf('parameter_sensitivity/Sd/Sd_%.1g', Sds(3));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig11 = figure(11);
    fig11.Position = [100, 100, 1200, 1200];

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
 
folderName = sprintf('parameter_sensitivity/Sd/Sd_%.1g', Sds(1));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;

max_time = size(numerics.tumor_num, 2);
time_vector = (time_init + (1:max_time)) / 10;
fig12 = figure(12);

    fig12.Position = [100, 100, 1200, 400];
    title('\mu', 'Interpreter', 'latex');
    subplot(3, 2, 1);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 2);
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


    folderName = sprintf('parameter_sensitivity/Sd/Sd_%.1g', Sds(2));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;

    subplot(3, 2, 3);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 4);
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

    folderName = sprintf('parameter_sensitivity/Sd/Sd_%.1g', Sds(3));
    % folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(3));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;


    subplot(3, 2, 5);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 6);
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

        %%
    folderName = sprintf('parameter_sensitivity/Dd/Dd_%.2g', Dds(1));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig13 = figure(13);
    fig13.Position = [100, 100, 1200, 1200];

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

folderName = sprintf('parameter_sensitivity/Dd/Dd_%.1g', Dds(2));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig14 = figure(14);
    fig14.Position = [100, 100, 1200, 1200];

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
    
folderName = sprintf('parameter_sensitivity/Dd/Dd_%.1g', Dds(3));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig15 = figure(15);
    fig15.Position = [100, 100, 1200, 1200];

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
 
folderName = sprintf('parameter_sensitivity/Dd/Dd_%.2g', Dds(1));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;

max_time = size(numerics.tumor_num, 2);
time_vector = (time_init + (1:max_time)) / 10;
fig16 = figure(16);

    fig16.Position = [100, 100, 1200, 400];
    title('\mu', 'Interpreter', 'latex');
    subplot(3, 2, 1);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 2);
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


    folderName = sprintf('parameter_sensitivity/Dd/Dd_%.1g', Dds(2));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;

    subplot(3, 2, 3);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 4);
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

    folderName = sprintf('parameter_sensitivity/Dd/Dd_%.1g', Dds(3));
    % folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(3));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;


    subplot(3, 2, 5);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 6);
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

            %%
    folderName = sprintf('parameter_sensitivity/So/So_%.2g', Sos(1));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig17 = figure(17);
    fig17.Position = [100, 100, 1200, 1200];

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

folderName = sprintf('parameter_sensitivity/So/So_%.2g', Sos(2));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig18 = figure(18);
    fig18.Position = [100, 100, 1200, 1200];

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
    
folderName = sprintf('parameter_sensitivity/So/So_%.2g', Sos(3));
snapshotsFile = fullfile(folderName, 'snapshots.mat');
snapshots_data = load(snapshotsFile);
snapshots = snapshots_data.snapshots;

infoFile = fullfile(folderName, 'info.mat');
info_data = load(infoFile);
info = info_data.info;

max_time = 400;
div_num = 4;
plot_int = max_time / div_num;
params = info.params;

fig19 = figure(19);
    fig19.Position = [100, 100, 1200, 1200];

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
 
folderName = sprintf('parameter_sensitivity/So/So_%.2g', Sos(1));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;

max_time = size(numerics.tumor_num, 2);
time_vector = (time_init + (1:max_time)) / 10;
fig20 = figure(20);

    fig20.Position = [100, 100, 1200, 400];
    title('\mu', 'Interpreter', 'latex');
    subplot(3, 2, 1);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 2);
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


    folderName = sprintf('parameter_sensitivity/So/So_%.2g', Sos(2));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;

    subplot(3, 2, 3);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 4);
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

    folderName = sprintf('parameter_sensitivity/So/So_%.2g', Sos(3));
    % folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(3));
numericsFile = fullfile(folderName, 'numerics.mat');
numerics_data = load(numericsFile);
numerics = numerics_data.numerics;
dam_accum = numerics.dam_accum;
death_thres = numerics.death_thres;
hypoxic_num = numerics.hypoxic_num;
normorxic_num = numerics.normorxic_num;
tumor_num = numerics.tumor_num;


    subplot(3, 2, 5);
    plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
    hold on;
    plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
    plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
    hold off;
    legend('FontSize', 8, 'Location', 'best');
    xlabel('t', 'FontSize', 12);
    ylabel('Cell Number', 'FontSize', 12);
    title('Tumor Cell Number', 'FontSize', 12);

    subplot(3, 2, 6);
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