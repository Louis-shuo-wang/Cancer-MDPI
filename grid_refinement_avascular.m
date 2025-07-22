<<<<<<< HEAD
clear; close all;

folderName = 'Figures/2025_Jul/vascular0/grid_refinement';

time_init = 0;
save_path = folderName;
hList   = [0.02, 0.01, 0.005];

params = initializeParameters();
max_time = 300;
time = (1:max_time)*params.dt;

params.age = 1e5;   % Stop tumor proliferation
% This is where we make modifications
% in this simulation, we apply no drug
params.Sd = 0;
params.initial_tumor_number = 0;
params.initial_tip_number = 0;
diffSum = zeros(max_time, numel(hList));
numericS = cell(1,numel(hList));

for k = 1:numel(hList)
    h = hList(k);
    % grid and temporal spacings

    % % start parallel pool if not already running
    % if isempty(gcp('nocreate'))
    %     parpool;
    % end
    params.grid_size = [1/h, 1/h];        % Size of simulation domain
    params.dg = 1 ./ params.grid_size;
    params.dy = params.dg(1);
    params.dx = params.dg(2);
    initsetting  = struct();
    initsetting.info = [];
    initsetting.params = params;
    [snapshots, numerics, info] = main(max_time, folderName, [], initsetting);
    diffSum(:,k) = numerics.diffSum(:);
    numericS{1,k}.snapshots = snapshots;
    numericS{1,k}.info = info;
    numericS{1,k}.numerics = numerics;
    fprintf('h=%.4g → maxΔsum=%.3e\n', h, max(numerics.diffSum(:)));
end

% Ensure save directory exists
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

% plot convergence
figure; hold on;
plot(time, log10(diffSum(:,1)), '-b', 'LineWidth',1, 'DisplayName', sprintf('h = %.4g',hList(1)));
plot(time, log10(diffSum(:,2)), '--r', 'LineWidth',1, 'DisplayName', sprintf('h = %.4g',hList(2)));
plot(time, log10(diffSum(:,3)), ':k', 'LineWidth',1, 'DisplayName', sprintf('h = %.4g',hList(3)));
hold off;
grid on;
xlabel('t');
ylabel('$\log_{10}\left( \sum_{u} \max_{i,j} \left| u_{i,j}^{k+1} - u_{i,j}^{k} \right| \right)$', 'Interpreter', 'latex');
title('Step‐change sum vs. time for different h');
legend('Location','best');

% Save figures and data
svg_file = fullfile(save_path, 'figure_stepchange_spatial_noABM.svg');
fig_file = fullfile(save_path, 'figure_stepchange_spatial_noABM.fig');

print(gcf, svg_file, '-dsvg');
savefig(fig_file);

data_file = fullfile(save_path, 'figure_stepchange_spatial_noABM_data.mat');
=======
clear; close all;

folderName = 'Figures/2025_Jul/vascular0/grid_refinement';

time_init = 0;
save_path = folderName;
hList   = [0.02, 0.01, 0.005];

params = initializeParameters();
max_time = 300;
time = (1:max_time)*params.dt;

params.age = 1e5;   % Stop tumor proliferation
% This is where we make modifications
% in this simulation, we apply no drug
params.Sd = 0;
params.initial_tumor_number = 0;
params.initial_tip_number = 0;
diffSum = zeros(max_time, numel(hList));
numericS = cell(1,numel(hList));

for k = 1:numel(hList)
    h = hList(k);
    % grid and temporal spacings

    % % start parallel pool if not already running
    % if isempty(gcp('nocreate'))
    %     parpool;
    % end
    params.grid_size = [1/h, 1/h];        % Size of simulation domain
    params.dg = 1 ./ params.grid_size;
    params.dy = params.dg(1);
    params.dx = params.dg(2);
    initsetting  = struct();
    initsetting.info = [];
    initsetting.params = params;
    [snapshots, numerics, info] = main(max_time, folderName, [], initsetting);
    diffSum(:,k) = numerics.diffSum(:);
    numericS{1,k}.snapshots = snapshots;
    numericS{1,k}.info = info;
    numericS{1,k}.numerics = numerics;
    fprintf('h=%.4g → maxΔsum=%.3e\n', h, max(numerics.diffSum(:)));
end

% Ensure save directory exists
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

% plot convergence
figure; hold on;
plot(time, log10(diffSum(:,1)), '-b', 'LineWidth',1, 'DisplayName', sprintf('h = %.4g',hList(1)));
plot(time, log10(diffSum(:,2)), '--r', 'LineWidth',1, 'DisplayName', sprintf('h = %.4g',hList(2)));
plot(time, log10(diffSum(:,3)), ':k', 'LineWidth',1, 'DisplayName', sprintf('h = %.4g',hList(3)));
hold off;
grid on;
xlabel('t');
ylabel('$\log_{10}\left( \sum_{u} \max_{i,j} \left| u_{i,j}^{k+1} - u_{i,j}^{k} \right| \right)$', 'Interpreter', 'latex');
title('Step‐change sum vs. time for different h');
legend('Location','best');

% Save figures and data
svg_file = fullfile(save_path, 'figure_stepchange_spatial_noABM.svg');
fig_file = fullfile(save_path, 'figure_stepchange_spatial_noABM.fig');

print(gcf, svg_file, '-dsvg');
savefig(fig_file);

data_file = fullfile(save_path, 'figure_stepchange_spatial_noABM_data.mat');
>>>>>>> d5569c1 (Initial commit)
save(data_file, 'diffSum', 'numericS', 'time', '-v7.3');