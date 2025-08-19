function plot_extended_therapy_paper()
    clear all; close all;

    
    % === Load common time info ===
    info = load('Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat', 'info');
    time_init = info.info.time;

    % === Collect results ===
    
    % results = struct();
    %     all_tumor = [];
    %     all_resis = [];
    %     all_hyp   = [];
max_time = 400;
        %% metronomic
            folderName = 'extended_therapy/metronomic';
             numericsFile = fullfile(folderName, 'numerics.mat');
          
            numerics = load(numericsFile).numerics;

            tumor_metro = numerics.tumor_num;
            resis_metro = numerics.resis_frac;
            hyp_metro   = numerics.r_hyp;

        %% adaptive
folderName = 'extended_therapy/adaptive';
             numericsFile = fullfile(folderName, 'numerics.mat');
          
            numerics = load(numericsFile).numerics;

            tumor_adap = numerics.tumor_num;
            resis_adap = numerics.resis_frac;
            hyp_adap   = numerics.r_hyp;

            %% antiangiogenic_preexisting
            folderName = 'extended_therapy/antiangiogenic_preexisting';
             numericsFile = fullfile(folderName, 'numerics.mat');
          
            numerics = load(numericsFile).numerics;

            tumor_pre = numerics.tumor_num;
            resis_pre = numerics.resis_frac;
            hyp_pre   = numerics.r_hyp;


            %% antiangiogenic_spontaneous
            folderName = 'extended_therapy/antiangiogenic_spontaneous';
             numericsFile = fullfile(folderName, 'numerics.mat');
          
            numerics = load(numericsFile).numerics;

            tumor_spon = numerics.tumor_num;
            resis_spon = numerics.resis_frac;
            hyp_spon   = numerics.r_hyp;


        time_vector = (time_init + (1:max_time)) / 10;

        % results.time = time_vector;
        % results.tumor.mean = mean(all_tumor, 1);
        % results.tumor.std  = std(all_tumor, 0, 1);
        % results.resis.mean = mean(all_resis, 1);
        % results.resis.std  = std(all_resis, 0, 1);
        % results.hyp.mean   = mean(all_hyp, 1);
        % results.hyp.std    = std(all_hyp, 0, 1);


% === Figure with 5x3 layout ===
fig = figure('Position',[100 100 1200 800]);
% Collect global y-limits per column
% ylims_tumor = get_global_ylim(results, param_labels, 'tumor');
% ylims_resis = get_global_ylim(results, param_labels, 'resis');
% ylims_hyp = get_global_ylim(results, param_labels, 'hyp');

tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact');


% --- Tumor ---
nexttile;
plot(time_vector, tumor_metro, 'r-', 'LineWidth', 1.5); hold on;
plot(time_vector, tumor_adap, 'b--', 'LineWidth', 1.5);
plot(time_vector, tumor_pre, 'g-.', 'LineWidth', 1.5);
plot(time_vector, tumor_spon, 'k:', 'LineWidth', 1.5);
set(gca,'FontSize',14);
xlabel('t','FontSize',18);
xlim([14,54]);
ylabel('number','FontSize',18);
title('Tumor Count','FontSize',18);

% --- Resistance ---
nexttile;
plot(time_vector, resis_metro, 'r-', 'LineWidth', 1.5); hold on;
plot(time_vector, resis_adap, 'b--', 'LineWidth', 1.5);
plot(time_vector, resis_pre, 'g-.', 'LineWidth', 1.5);
plot(time_vector, resis_spon, 'k:', 'LineWidth', 1.5);
set(gca,'FontSize',14);
xlabel('t','FontSize',18);
xlim([14,54]);
ylabel('fraction','FontSize',18);
title('Resistance Fraction','FontSize',18);

% % --- Hypoxic ---
% nexttile;
% plot(time_vector, hyp_metro, 'r-', 'LineWidth', 1.5); hold on;
% plot(time_vector, hyp_adap, 'b--', 'LineWidth', 1.5);
% plot(time_vector, hyp_pre, 'g-.', 'LineWidth', 1.5);
% plot(time_vector, hyp_spon, 'k:', 'LineWidth', 1.5);
% set(gca,'FontSize',14);
% xlabel('t','FontSize',18);
% ylabel('fraction','FontSize',18);
% title('Hypoxic Fraction','FontSize',18);

% Legend (only once)
lgd = legend({'metronomic','adaptive','com\_preexisting','com\_spontaneous'},...
    'FontSize',16,'Location','best');
end