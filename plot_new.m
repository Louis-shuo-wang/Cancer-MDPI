function plot_parameter_sensitivity()
    clear all; close all;

    % === Parameter sets ===
    param_sets = struct( ...
        'mutation_rate', [1e-4, 1e-3, 1e-2], ...
        'p_r',            [0.1, 0.2, 0.3], ...
        'S_d',            [1, 2, 4], ...
        'D_d',            [0.25, 0.5, 1], ...
        'S_o',            [2.8, 3.5, 4.2] ...
    );

    param_labels = fieldnames(param_sets);   % 5 rows
    colors = lines(numel(param_labels));     % distinguishable colors

    % === Load common time info ===
    info = load('Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat', 'info');
    time_init = info.info.time;

    % === Collect results ===
    results = struct();
    for p = 1:numel(param_labels)
        pname = param_labels{p};
        values = param_sets.(pname);

        all_tumor = [];
        all_resis = [];
        all_hyp   = [];

        for v = 1:numel(values)
            folderName = sprintf('parameter_sensitivity/%s/%s_%g', pname, pname, values(v));
            numericsFile = fullfile(folderName, 'numerics.mat');

            if ~isfile(numericsFile)
                warning('Missing file: %s', numericsFile);
                continue;
            end

            numerics = load(numericsFile).numerics;
            max_time = size(numerics.tumor_num, 2);

            all_tumor = [all_tumor; numerics.tumor_num];
            all_resis = [all_resis; numerics.resis_frac];
            all_hyp   = [all_hyp;   numerics.r_hyp];
        end

        time_vector = (time_init + (1:max_time)) / 10;

        results.(pname).time = time_vector;
        results.(pname).tumor.mean = mean(all_tumor, 1);
        results.(pname).tumor.std  = std(all_tumor, 0, 1);
        results.(pname).resis.mean = mean(all_resis, 1);
        results.(pname).resis.std  = std(all_resis, 0, 1);
        results.(pname).hyp.mean   = mean(all_hyp, 1);
        results.(pname).hyp.std    = std(all_hyp, 0, 1);
    end

    % === Figure with 5x3 layout ===
    figure('Position',[100 100 1500 1200]);

    % Collect global y-limits per column
    ylims_tumor = get_global_ylim(results, param_labels, 'tumor');
    ylims_resis = get_global_ylim(results, param_labels, 'resis');
    ylims_hyp   = get_global_ylim(results, param_labels, 'hyp');

    for r = 1:numel(param_labels)
        pname = param_labels{r};
        if strcmp(pname, 'mutation_rate')
            plabel = '\mu';
        else
            plabel = pname;
        end
        time = results.(pname).time;

        % --- Tumor (col 1) ---
        subplot(numel(param_labels), 3, (r-1)*3 + 1);
        shaded_plot(time, results.(pname).tumor.mean, results.(pname).tumor.std, colors(r,:));
        ylim(ylims_tumor);
        if r == 1, title('Tumor Count'); end
        ylabel(plabel, 'Interpreter','latex'); % only first column has ylabel
        if r < numel(param_labels), set(gca,'XTickLabel',[]); end

        % --- Resistance (col 2) ---
        subplot(numel(param_labels), 3, (r-1)*3 + 2);
        shaded_plot(time, results.(pname).resis.mean, results.(pname).resis.std, colors(r,:));
        ylim(ylims_resis);
        if r == 1, title('Resistance Fraction'); end
        if r < numel(param_labels), set(gca,'XTickLabel',[]); end
        set(gca,'YTickLabel',[]); % no ylabel

        % --- Hypoxic (col 3) ---
        subplot(numel(param_labels), 3, (r-1)*3 + 3);
        shaded_plot(time, results.(pname).hyp.mean, results.(pname).hyp.std, colors(r,:));
        ylim(ylims_hyp);
        if r == 1, title('Hypoxic Fraction'); end
        if r < numel(param_labels), set(gca,'XTickLabel',[]); end
        set(gca,'YTickLabel',[]); % no ylabel
    end
end

%% === Helper: plot mean ± std shaded ===
function shaded_plot(time, mu, sig, c)
    fill([time, fliplr(time)], [mu+sig, fliplr(mu-sig)], ...
        c, 'FaceAlpha',0.2, 'EdgeColor','none'); hold on;
    plot(time, mu, '-', 'Color', c, 'LineWidth', 1.5);
    grid on;
end

%% === Helper: compute global y-limits ===
function ylims = get_global_ylim(results, param_labels, field)
    all_vals = [];
    for r = 1:numel(param_labels)
        pname = param_labels{r};
        mu = results.(pname).(field).mean;
        sig = results.(pname).(field).std;
        all_vals = [all_vals, mu+sig, mu-sig];
    end
    ylims = [min(all_vals), max(all_vals)];
end

%%
% folderName = sprintf('parameter_sensitivity/mutation_rate/mutation_rate_%.1g', mutation_rates(1));
% snapshotsFile = fullfile(folderName, 'snapshots.mat');
% snapshots_data = load(snapshotsFile);
% snapshots = snapshots_data.snapshots;
% 
% infoFile = fullfile(folderName, 'info.mat');
% info_data = load(infoFile);
% info = info_data.info;
% 
% max_time = 400;
% div_num = 4;
% plot_int = max_time / div_num;
% params = info.params;
% 
% fig1 = figure(1);
%     fig1.Position = [100, 100, 1200, 1200];
% 
%     % Create a 4\time 4 grid of subplots
%     for k = 1:4
%         subplot(4, 4, (k - 1) * 4 + 1);
%         plotTumorAndVessels(snapshots{1}{k}.tumor_cells, ...
%             snapshots{1}{k}.vessel_agents, ...
%             snapshots{1}{k}.params);
%         if (snapshots{1}{k}.time - time_init) / plot_int == 1
%             title('vessel', 'FontSize', 12);
%         end
%         % Set axis labels and ticks
%         axis equal;
%         xlabel('x', 'FontSize', 12);
%         % Set the row title
%         ylabel(sprintf('Time: t = %.1f', snapshots{1}{k}.time * 0.1), 'FontSize', 12);
%         xlim([0, 1]);
%         ylim([0, 1]);
%         set(gca, 'XTick', 0:0.2:1, ...
%             'YTick', 0:0.2:1, ...
%             'FontSize', 10);
%     end
% 
%     for i = 1:4
%         for j = 2:4         % Loop over variables vessel, TAF, oxygen, drug (columns)
%             subplot_idx = (i - 1) * 4 + j;
%             subplot(4, 4, subplot_idx);
% 
%             % Get the corresponding snapshot
%             snapshot = snapshots{j}{i};
% 
%             % Display the snapshot
%             imagesc(snapshot);
%             colorbar;
%             axis tight equal;
%             colormap(gca, 'parula');
% 
%             % Set the title
%             if i == 1
%                 title(snapshots{5}{j}, 'FontSize', 12);
%             end
% 
%             % Set axis labels and ticks
%             xlabel('x', 'FontSize', 12);
%             ylabel('y', 'FontSize', 12);
%             xlim([0, params.grid_size(2)]);
%             ylim([0, params.grid_size(1)]);
%             set(gca, 'XTick', 0:params.grid_size(2) / 5:params.grid_size(2), 'XTickLabel', 0:0.2:1, ...
%                 'YTick', 0:params.grid_size(1) / 5:params.grid_size(1), 'YTickLabel', 1:-0.2:0, ...
%                 'FontSize', 10);
%         end
%     end
% 
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
% max_time = size(numerics.tumor_num, 2);
% time_vector = (time_init + (1:max_time)) / 10;
% fig2 = figure(2);
% 
%     fig2.Position = [100, 100, 1200, 400];
%     title('\mu', 'Interpreter', 'latex');
%     subplot(3, 2, 1);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 2);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
% 
% 
%     folderName = sprintf('parameter_sensitivity/mutation_rate/mutation_rate_%.1g', mutation_rates(2));
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
%     subplot(3, 2, 3);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 4);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
% 
%     folderName = sprintf('parameter_sensitivity/mutation_rate/mutation_rate_%.1g', mutation_rates(3));
%     % folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(3));
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
% 
%     subplot(3, 2, 5);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 6);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
% 
%     filePath = fullfile(folderName, 'num1');
%     set(fig2, 'toolbar', 'none');
%     saveas(fig2, [filePath, '.fig'], 'fig');
%     % saveas(fig2, [filePath, '.eps'], 'eps');
%     saveas(fig2, [filePath, '.svg'], 'svg');

% %%
% folderName = sprintf('parameter_sensitivity/Dd/Dd_%.2g', Dds(1));
% snapshotsFile = fullfile(folderName, 'snapshots.mat');
% snapshots_data = load(snapshotsFile);
% snapshots = snapshots_data.snapshots;
% 
% infoFile = fullfile(folderName, 'info.mat');
% info_data = load(infoFile);
% info = info_data.info;
% 
% max_time = 400;
% div_num = 4;
% plot_int = max_time / div_num;
% params = info.params;
% 
% fig3 = figure(3);
%     fig3.Position = [100, 100, 1200, 1200];
% 
%     % Create a 4\time 4 grid of subplots
%     for k = 1:4
%         subplot(4, 4, (k - 1) * 4 + 1);
%         plotTumorAndVessels(snapshots{1}{k}.tumor_cells, ...
%             snapshots{1}{k}.vessel_agents, ...
%             snapshots{1}{k}.params);
%         if (snapshots{1}{k}.time - time_init) / plot_int == 1
%             title('vessel', 'FontSize', 12);
%         end
%         % Set axis labels and ticks
%         axis equal;
%         xlabel('x', 'FontSize', 12);
%         % Set the row title
%         ylabel(sprintf('Time: t = %.1f', snapshots{1}{k}.time * 0.1), 'FontSize', 12);
%         xlim([0, 1]);
%         ylim([0, 1]);
%         set(gca, 'XTick', 0:0.2:1, ...
%             'YTick', 0:0.2:1, ...
%             'FontSize', 10);
%     end
% 
%     for i = 1:4
%         for j = 2:4         % Loop over variables vessel, TAF, oxygen, drug (columns)
%             subplot_idx = (i - 1) * 4 + j;
%             subplot(4, 4, subplot_idx);
% 
%             % Get the corresponding snapshot
%             snapshot = snapshots{j}{i};
% 
%             % Display the snapshot
%             imagesc(snapshot);
%             colorbar;
%             axis tight equal;
%             colormap(gca, 'parula');
% 
%             % Set the title
%             if i == 1
%                 title(snapshots{5}{j}, 'FontSize', 12);
%             end
% 
%             % Set axis labels and ticks
%             xlabel('x', 'FontSize', 12);
%             ylabel('y', 'FontSize', 12);
%             xlim([0, params.grid_size(2)]);
%             ylim([0, params.grid_size(1)]);
%             set(gca, 'XTick', 0:params.grid_size(2) / 5:params.grid_size(2), 'XTickLabel', 0:0.2:1, ...
%                 'YTick', 0:params.grid_size(1) / 5:params.grid_size(1), 'YTickLabel', 1:-0.2:0, ...
%                 'FontSize', 10);
%         end
%     end
% 
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
% max_time = size(numerics.tumor_num, 2);
% time_vector = (time_init + (1:max_time)) / 10;
% fig4 = figure(4);
% 
%     fig4.Position = [100, 100, 1200, 400];
%     title('D_d', 'Interpreter', 'latex');
%     subplot(3, 2, 1);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 2);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
% 
% 
%     folderName = sprintf('parameter_sensitivity/Dd/Dd_%.1g', Dds(2));
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
%     subplot(3, 2, 3);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 4);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
% 
%     folderName = sprintf('parameter_sensitivity/Dd/Dd_%.1g', Dds(3));
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
% 
%     subplot(3, 2, 5);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 6);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
% 
%     %%
%     folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(1));
% snapshotsFile = fullfile(folderName, 'snapshots.mat');
% snapshots_data = load(snapshotsFile);
% snapshots = snapshots_data.snapshots;
% 
% infoFile = fullfile(folderName, 'info.mat');
% info_data = load(infoFile);
% info = info_data.info;
% 
% max_time = 400;
% div_num = 4;
% plot_int = max_time / div_num;
% params = info.params;
% 
% fig5 = figure(5);
%     fig5.Position = [100, 100, 1200, 1200];
% 
%     % Create a 4\time 4 grid of subplots
%     for k = 1:4
%         subplot(4, 4, (k - 1) * 4 + 1);
%         plotTumorAndVessels(snapshots{1}{k}.tumor_cells, ...
%             snapshots{1}{k}.vessel_agents, ...
%             snapshots{1}{k}.params);
%         if (snapshots{1}{k}.time - time_init) / plot_int == 1
%             title('vessel', 'FontSize', 12);
%         end
%         % Set axis labels and ticks
%         axis equal;
%         xlabel('x', 'FontSize', 12);
%         % Set the row title
%         ylabel(sprintf('Time: t = %.1f', snapshots{1}{k}.time * 0.1), 'FontSize', 12);
%         xlim([0, 1]);
%         ylim([0, 1]);
%         set(gca, 'XTick', 0:0.2:1, ...
%             'YTick', 0:0.2:1, ...
%             'FontSize', 10);
%     end
% 
%     for i = 1:4
%         for j = 2:4         % Loop over variables vessel, TAF, oxygen, drug (columns)
%             subplot_idx = (i - 1) * 4 + j;
%             subplot(4, 4, subplot_idx);
% 
%             % Get the corresponding snapshot
%             snapshot = snapshots{j}{i};
% 
%             % Display the snapshot
%             imagesc(snapshot);
%             colorbar;
%             axis tight equal;
%             colormap(gca, 'parula');
% 
%             % Set the title
%             if i == 1
%                 title(snapshots{5}{j}, 'FontSize', 12);
%             end
% 
%             % Set axis labels and ticks
%             xlabel('x', 'FontSize', 12);
%             ylabel('y', 'FontSize', 12);
%             xlim([0, params.grid_size(2)]);
%             ylim([0, params.grid_size(1)]);
%             set(gca, 'XTick', 0:params.grid_size(2) / 5:params.grid_size(2), 'XTickLabel', 0:0.2:1, ...
%                 'YTick', 0:params.grid_size(1) / 5:params.grid_size(1), 'YTickLabel', 1:-0.2:0, ...
%                 'FontSize', 10);
%         end
%     end
% 
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
% max_time = size(numerics.tumor_num, 2);
% time_vector = (time_init + (1:max_time)) / 10;
% fig6 = figure(6);
% 
%     fig6.Position = [100, 100, 1200, 400];
%     title('p_r', 'Interpreter', 'latex');
%     subplot(3, 2, 1);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 2);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
% 
% 
%     folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(2));
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
%     subplot(3, 2, 3);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 4);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
% 
%     folderName = sprintf('parameter_sensitivity/pr/pr_%.1g', prs(3));
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
% 
%     subplot(3, 2, 5);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 6);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
% 
%     %%
% 
%     folderName = sprintf('parameter_sensitivity/Sd/Sd_%.2g', Sds(1));
% snapshotsFile = fullfile(folderName, 'snapshots.mat');
% snapshots_data = load(snapshotsFile);
% snapshots = snapshots_data.snapshots;
% 
% infoFile = fullfile(folderName, 'info.mat');
% info_data = load(infoFile);
% info = info_data.info;
% 
% max_time = 400;
% div_num = 4;
% plot_int = max_time / div_num;
% params = info.params;
% 
% fig7 = figure(7);
%     fig7.Position = [100, 100, 1200, 1200];
% 
%     % Create a 4\time 4 grid of subplots
%     for k = 1:4
%         subplot(4, 4, (k - 1) * 4 + 1);
%         plotTumorAndVessels(snapshots{1}{k}.tumor_cells, ...
%             snapshots{1}{k}.vessel_agents, ...
%             snapshots{1}{k}.params);
%         if (snapshots{1}{k}.time - time_init) / plot_int == 1
%             title('vessel', 'FontSize', 12);
%         end
%         % Set axis labels and ticks
%         axis equal;
%         xlabel('x', 'FontSize', 12);
%         % Set the row title
%         ylabel(sprintf('Time: t = %.1f', snapshots{1}{k}.time * 0.1), 'FontSize', 12);
%         xlim([0, 1]);
%         ylim([0, 1]);
%         set(gca, 'XTick', 0:0.2:1, ...
%             'YTick', 0:0.2:1, ...
%             'FontSize', 10);
%     end
% 
%     for i = 1:4
%         for j = 2:4         % Loop over variables vessel, TAF, oxygen, drug (columns)
%             subplot_idx = (i - 1) * 4 + j;
%             subplot(4, 4, subplot_idx);
% 
%             % Get the corresponding snapshot
%             snapshot = snapshots{j}{i};
% 
%             % Display the snapshot
%             imagesc(snapshot);
%             colorbar;
%             axis tight equal;
%             colormap(gca, 'parula');
% 
%             % Set the title
%             if i == 1
%                 title(snapshots{5}{j}, 'FontSize', 12);
%             end
% 
%             % Set axis labels and ticks
%             xlabel('x', 'FontSize', 12);
%             ylabel('y', 'FontSize', 12);
%             xlim([0, params.grid_size(2)]);
%             ylim([0, params.grid_size(1)]);
%             set(gca, 'XTick', 0:params.grid_size(2) / 5:params.grid_size(2), 'XTickLabel', 0:0.2:1, ...
%                 'YTick', 0:params.grid_size(1) / 5:params.grid_size(1), 'YTickLabel', 1:-0.2:0, ...
%                 'FontSize', 10);
%         end
%     end
% 
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
% max_time = size(numerics.tumor_num, 2);
% time_vector = (time_init + (1:max_time)) / 10;
% fig8 = figure(8);
%     fig8.Position = [100, 100, 1200, 400];
%     title('S_d', 'Interpreter', 'latex');
%     subplot(3, 2, 1);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 2);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
% 
% 
%     folderName = sprintf('parameter_sensitivity/Sd/Sd_%.1g', Sds(2));
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
%     subplot(3, 2, 3);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 4);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
% 
%     folderName = sprintf('parameter_sensitivity/Sd/Sd_%.1g', Sds(3));
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
% 
%     subplot(3, 2, 5);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 6);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
% 
%     %%
% 
% folderName = sprintf('parameter_sensitivity/So/So_%.2g', Sos(1));
% snapshotsFile = fullfile(folderName, 'snapshots.mat');
% snapshots_data = load(snapshotsFile);
% snapshots = snapshots_data.snapshots;
% 
% infoFile = fullfile(folderName, 'info.mat');
% info_data = load(infoFile);
% info = info_data.info;
% 
% max_time = 400;
% div_num = 4;
% plot_int = max_time / div_num;
% params = info.params;
% 
% fig9 = figure(9);
%     fig9.Position = [100, 100, 1200, 1200];
% 
%     % Create a 4\time 4 grid of subplots
%     for k = 1:4
%         subplot(4, 4, (k - 1) * 4 + 1);
%         plotTumorAndVessels(snapshots{1}{k}.tumor_cells, ...
%             snapshots{1}{k}.vessel_agents, ...
%             snapshots{1}{k}.params);
%         if (snapshots{1}{k}.time - time_init) / plot_int == 1
%             title('vessel', 'FontSize', 12);
%         end
%         % Set axis labels and ticks
%         axis equal;
%         xlabel('x', 'FontSize', 12);
%         % Set the row title
%         ylabel(sprintf('Time: t = %.1f', snapshots{1}{k}.time * 0.1), 'FontSize', 12);
%         xlim([0, 1]);
%         ylim([0, 1]);
%         set(gca, 'XTick', 0:0.2:1, ...
%             'YTick', 0:0.2:1, ...
%             'FontSize', 10);
%     end
% 
%     for i = 1:4
%         for j = 2:4         % Loop over variables vessel, TAF, oxygen, drug (columns)
%             subplot_idx = (i - 1) * 4 + j;
%             subplot(4, 4, subplot_idx);
% 
%             % Get the corresponding snapshot
%             snapshot = snapshots{j}{i};
% 
%             % Display the snapshot
%             imagesc(snapshot);
%             colorbar;
%             axis tight equal;
%             colormap(gca, 'parula');
% 
%             % Set the title
%             if i == 1
%                 title(snapshots{5}{j}, 'FontSize', 12);
%             end
% 
%             % Set axis labels and ticks
%             xlabel('x', 'FontSize', 12);
%             ylabel('y', 'FontSize', 12);
%             xlim([0, params.grid_size(2)]);
%             ylim([0, params.grid_size(1)]);
%             set(gca, 'XTick', 0:params.grid_size(2) / 5:params.grid_size(2), 'XTickLabel', 0:0.2:1, ...
%                 'YTick', 0:params.grid_size(1) / 5:params.grid_size(1), 'YTickLabel', 1:-0.2:0, ...
%                 'FontSize', 10);
%         end
%     end
% 
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
% max_time = size(numerics.tumor_num, 2);
% time_vector = (time_init + (1:max_time)) / 10;
% fig10 = figure(10);
%     fig5.Position = [100, 100, 1200, 400];
%     title('S_o', 'Interpreter', 'latex');
%     subplot(3, 2, 1);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 2);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
% 
% 
%     folderName = sprintf('parameter_sensitivity/So/So_%.2g', Sos(2));
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
%     subplot(3, 2, 3);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 4);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
% 
%     folderName = sprintf('parameter_sensitivity/So/SO_%.2g', Sos(3));
% numericsFile = fullfile(folderName, 'numerics.mat');
% numerics_data = load(numericsFile);
% numerics = numerics_data.numerics;
% dam_accum = numerics.dam_accum;
% death_thres = numerics.death_thres;
% hypoxic_num = numerics.hypoxic_num;
% normorxic_num = numerics.normorxic_num;
% tumor_num = numerics.tumor_num;
% 
% 
%     subplot(3, 2, 5);
%     plot(time_vector, hypoxic_num, 'ro-', 'MarkerSize', 2.5, 'DisplayName', 'Hypoxic');
%     hold on;
%     plot(time_vector, normorxic_num, 'kx--', 'MarkerSize', 2.5, 'DisplayName', 'Normorxic');
%     plot(time_vector, tumor_num, 'b-', 'DisplayName', 'Tumor');
%     hold off;
%     legend('FontSize', 8, 'Location', 'best');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Cell Number', 'FontSize', 12);
%     title('Tumor Cell Number', 'FontSize', 12);
% 
%     subplot(3, 2, 6);
%     hold on;
%     % plot averaged damage accumulation with error bar
%     errorbar(time_vector, dam_accum(1, :), dam_accum(2, :), 'b-', ...
%         'DisplayName', 'Damage Accumulation');
%     errorbar(time_vector, death_thres(1, :), death_thres(2, :), 'r-', ...
%         'DisplayName', 'Death Threshold');
%     xlabel('t', 'FontSize', 12);
%     ylabel('Averaged Value', 'FontSize', 12);
%     title(sprintf('Averaged Damage Accumulation & \n Averaged Death Threshold'), 'FontSize', 12);
%     legend('Location', 'best', 'FontSize', 8);
%     hold off;
