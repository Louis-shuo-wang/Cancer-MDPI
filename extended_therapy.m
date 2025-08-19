clear all; close all;

filename = 'Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat';
info = load(filename, 'info').info;
time_init = info.time;
params = initializeParameters();

mainFolder = 'extended_therapy';
if ~exist(mainFolder, 'dir')
    mkdir(mainFolder);
end

max_time = 400;

for i = 3
%%  metronomic therapy
if i == 1
params.Sd = 2;  % Adjusted to 2 for consistency with low-dose baseline, avoiding failure at lower rates
params.mutation_rate = 0;
tumor_cells = info.tumor_cells;
    nonempty_idx = find(~cellfun('isempty', tumor_cells));
    num_resis = max(1, ceil(params.resis_fraction * length(nonempty_idx)));
    resis_idx = randsample(nonempty_idx, num_resis, false);

    for j = 1:length(resis_idx)
        idx = resis_idx(j);
        tumor_cells{idx}.death_threshold = params.base_death_threshold * params.resistance_factor;
    end
    info.tumor_cells = tumor_cells;

    initsetting = struct();
    initsetting.info = info;
    initsetting.params = params;

    folderName = sprintf('%s/%s', mainFolder, 'metronomic');
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end

    [snapshots, numerics, info] = main(max_time, folderName, [], initsetting);

        save(fullfile(folderName, 'snapshots.mat'), 'snapshots');
        save(fullfile(folderName, 'info.mat'), 'info');
        save(fullfile(folderName, 'params.mat'), 'params');
        save(fullfile(folderName, 'numerics.mat'), 'numerics');

%% adaptive therapy
elseif i == 2
params.mutation_rate = 1e-2;
tumor_cells = info.tumor_cells;
num_tumor = sum(~cellfun('isempty', tumor_cells));

    initsetting = struct();
    initsetting.info = info;
    initsetting.params = params;

   folderName = sprintf('%s/%s', mainFolder, 'adaptive');
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end

[snapshots, numerics, info] = main_extended_therapy(max_time, folderName, [], initsetting);

        save(fullfile(folderName, 'snapshots.mat'), 'snapshots');
        save(fullfile(folderName, 'info.mat'), 'info');
        save(fullfile(folderName, 'params.mat'), 'params');
        save(fullfile(folderName, 'numerics.mat'), 'numerics');

% Anti-angiogenic therapy (spontaneous, low-mu)
elseif i == 3
    params.Sd = 5;
    params.So = 2.45;  % Reduced for normalization
    params.cbr = 0.5;  % Conceptual reduction for any residual branching
    params.mutation_rate = 1e-4;  % Low-mu for spontaneous

    % Add vessel pruning for anti-angiogenic effect on static network
    % Assuming info.V_t is the cell array or list of vessel agents
    nonempty_v = find(~cellfun('isempty', info.vessel_agents));
    num_v = length(nonempty_v);
    prune_prob = 0.3;  % Prune 30% of vessels
    prune_idx = rand(num_v,1) < prune_prob;
    prune_selected = nonempty_v(prune_idx);
    info.vessel_agents(prune_selected) = {[]};  % Remove pruned vessels (set to empty)

    % Increase drug diffusion for improved penetration post-normalization
    params.Dd = params.Dd * 1.2;  % 20% increase

    initsetting = struct();
    initsetting.info = info;
    initsetting.params = params;

 folderName = sprintf('%s/%s_%s', mainFolder, 'antiangiogenic', 'spontaneous');
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end
  [snapshots, numerics, info] = main_extended_therapy(max_time, folderName, [], initsetting);

        save(fullfile(folderName, 'snapshots.mat'), 'snapshots');
        save(fullfile(folderName, 'info.mat'), 'info');
        save(fullfile(folderName, 'params.mat'), 'params');
        save(fullfile(folderName, 'numerics.mat'), 'numerics');

% Anti-angiogenic therapy (preexisting)
else

    params.Sd = 5;
    params.So = 2.45;  % Reduced for normalization
    params.cbr = 0.5;  % Conceptual reduction
    params.mutation_rate = 0;

    tumor_cells = info.tumor_cells;
    nonempty_idx = find(~cellfun('isempty', tumor_cells));
    num_resis = max(1, ceil(params.resis_fraction * length(nonempty_idx)));
    resis_idx = randsample(nonempty_idx, num_resis, false);

    for j = 1:length(resis_idx)
        idx = resis_idx(j);
        tumor_cells{idx}.death_threshold = params.base_death_threshold * params.resistance_factor;
    end
    info.tumor_cells = tumor_cells;

    % Add vessel pruning for anti-angiogenic effect on static network
    nonempty_v = find(~cellfun('isempty', info.vessel_agents));
    num_v = length(nonempty_v);
    prune_prob = 0.3;  % Prune 30% of vessels
    prune_idx = rand(num_v,1) < prune_prob;
    prune_selected = nonempty_v(prune_idx);
    info.vessel_agents(prune_selected) = {[]};  % Remove pruned vessels

    % Increase drug diffusion for improved penetration
    params.Dd = params.Dd * 1.2;  % 20% increase

    initsetting = struct();
    initsetting.info = info;
    initsetting.params = params;

  folderName = sprintf('%s/%s_%s', mainFolder, 'antiangiogenic', 'preexisting');
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end
  [snapshots, numerics, info] = main_extended_therapy(max_time, folderName, [], initsetting);

        save(fullfile(folderName, 'snapshots.mat'), 'snapshots');
        save(fullfile(folderName, 'info.mat'), 'info');
        save(fullfile(folderName, 'params.mat'), 'params');
        save(fullfile(folderName, 'numerics.mat'), 'numerics');

end
end