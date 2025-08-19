function plot_parameter_sensitivity_paper()
    clear all; close all;

    % === Parameter sets ===
    param_sets = struct( ...
        'mutation_rate', [1e-4, 1e-3, 1e-2], ...
        'pr',            [0.1, 0.2, 0.3], ...
        'Sd',            [1, 2, 4], ...
        'Dd',            [0.25, 0.5, 1], ...
        'So',            [2.8, 3.5, 4.2] ...
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
            if p == 4 
                if v == 1
            folderName = sprintf('parameter_sensitivity/%s/%s_%.2g', pname, pname, values(v));
             numericsFile = fullfile(folderName, 'numerics.mat');
                else
            folderName = sprintf('parameter_sensitivity/%s/%s_%.1g', pname, pname, values(v));    
             numericsFile = fullfile(folderName, 'numerics.mat');
                end
            elseif p == 5
                folderName = sprintf('parameter_sensitivity/%s/%s_%.2g', pname, pname, values(v));
                 numericsFile = fullfile(folderName, 'numerics.mat');
            else 
                folderName = sprintf('parameter_sensitivity/%s/%s_%.1g', pname, pname, values(v));
            numericsFile = fullfile(folderName, 'numerics.mat');

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
ylims_hyp = get_global_ylim(results, param_labels, 'hyp');

tiledlayout(numel(param_labels), 3, 'TileSpacing','compact', 'Padding','compact');

for r = 1:numel(param_labels)
    pname = param_labels{r};
    if strcmp(pname, 'mutation_rate')
        plabel = '$\mu$';
    else
        plabel = pname;
    end
    time = results.(pname).time;

    % --- Tumor ---
    nexttile;
    shaded_plot(time, results.(pname).tumor.mean, results.(pname).tumor.std, colors(r,:));
    xlim([14,54]); ylim(ylims_tumor);
    set(gca,'FontSize',14);  % xtick, ytick fontsize
    if r < numel(param_labels)
        set(gca,'XTickLabel',[]);
    else
        xlabel('t','FontSize',18);
    end
    if r == 1, title('Tumor Count','FontSize',18); end
    ylabel(plabel,'Interpreter','latex','FontSize',16);
ylabel(plabel,'Interpreter','latex','FontSize',18);

    % --- Resistance ---
    nexttile;
    shaded_plot(time, results.(pname).resis.mean, results.(pname).resis.std, colors(r,:));
    xlim([14,54]); ylim(ylims_resis);
    set(gca,'FontSize',14);
    if r < numel(param_labels)
        set(gca,'XTickLabel',[]);
    else
        xlabel('t','FontSize',18);
    end
    if r == 1, title('Resistance Fraction','FontSize',18); end
ylabel(plabel,'Interpreter','latex','FontSize',18);

    % --- Hypoxic ---
    nexttile;
    shaded_plot(time, results.(pname).hyp.mean, results.(pname).hyp.std, colors(r,:));
    xlim([14,54]); ylim(ylims_hyp);
    set(gca,'FontSize',14);
    if r < numel(param_labels)
        set(gca,'XTickLabel',[]);
    else
        xlabel('t','FontSize',18);
    end
    if r == 1
        title('Hypoxic Fraction','FontSize',18);
    end
    lgd = legend({'Mean \pm Std','Mean'},'FontSize',18,'Location','best');
    ylabel(plabel,'Interpreter','latex','FontSize',18);
end
%% === Helper: plot mean ± std shaded ===
function shaded_plot(time, mu, sig, c)
time = time(:).';
mu = mu(:).';
sig = sig(:).';
fill([time, fliplr(time)], [mu+sig, fliplr(mu-sig)], ...
            c, 'FaceAlpha',0.2, 'EdgeColor','none'); hold on;
        plot(time, mu, '-', 'Color', c, 'LineWidth', 1.8);
grid on;
end
%% === Helper: compute global y-limits ===
function ylims = get_global_ylim(results, param_labels, field)
all_vals = [];
for r = 1:numel(param_labels)
pname = param_labels{r};
mu = results.(pname).(field).mean;
sig = results.(pname).(field).std;
mu = mu(:).';
sig = sig(:).';
all_vals = [all_vals, mu+sig, mu-sig];
end
ylims = [min(all_vals), max(all_vals)];
end

end