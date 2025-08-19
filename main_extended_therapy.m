function [all_snapshots, numerics, info] = main_extended_therapy(max_time, folderName, treatment, initsetting)

close all;
tic;

if nargin < 3
    treatment = [];
    params = initializeParameters();
    [tumor_cells, tip_cells, vessel_agents, chem_field] = initializeAgents(params);
    time = 0;
else
    if ~isempty(treatment)
        treat_on = treatment(1);
        dose = treatment(2);
    end

    info = initsetting.info;
    params = initsetting.params;

    if isempty(params)
        params = initializeParameters();
    end

    if ~isempty(info)
        tumor_cells = info.tumor_cells;
        tip_cells = info.tip_cells;
        vessel_agents = info.vessel_agents;
        chem_field = info.chem_field;
        time = info.time;
    else
        [tumor_cells, tip_cells, vessel_agents, chem_field] = initializeAgents(params);
        time = 0;
    end
end

if ~exist(folderName, 'dir')
    mkdir(folderName);
end

init_num_tumor = sum(~cellfun('isempty', tumor_cells));


div_num = 4;
plot_int = max_time / div_num;

variables = {'vessel', 'TAF', 'oxygen', 'drug'};

vessel_snapshots = cell(1, 4);
TAF_snapshots = cell(1, 4);
oxygen_snapshots = cell(1, 4);
drug_snapshots = cell(1, 4);
Tumor_cells = cell(1, 4);
Tip_cells = cell(1, 4);
Vessel_agents = cell(1, 4);
Chem_field = cell(1, 4);

dam_accum = zeros(2, max_time);
death_thres = zeros(2, max_time);
oxygen_consumption = cell(1, max_time);
oxygen_consumption_hist = zeros(max_time, params.nbins);
proliferation_rate = cell(1, max_time);
proliferation_rate_hist = zeros(max_time, params.nbins);

hypoxic_num = zeros(1, max_time);
normorxic_num = zeros(1, max_time);
tumor_num = zeros(1, max_time);
resis_frac = zeros(1, max_time);

peak_tumor_num = init_num_tumor;       % 历史峰值
dose_state = 'on';                     % 状态: 'on' = 给药, 'skip' = 暂停

for t_loc = 1:max_time
    fprintf('time step: %d\n', t_loc);
    t_glob = t_loc + time;

    if ~isempty(treatment)
        if mod(t_loc, 50) <= treat_on
            params.Sd = dose;
        else
            params.Sd = 0;
        end
    end

    %% Key update for PDE ABM variable
    [tumor_cells, tip_cells, vessel_agents] = updatePositions(tumor_cells, tip_cells, vessel_agents, chem_field, params);

    for n = 1:params.div_time
        chem_field = updateChemicalEnvironmentAngiogenesis(tumor_cells, vessel_agents, chem_field, params);
    end

    [tip_cells, vessel_agents] = updateAngiogenesis(tip_cells, vessel_agents, chem_field, t_glob, params);

    for n = 1:params.div_time
        chem_field = updateChemicalEnvironment(tumor_cells, vessel_agents, chem_field, params);
    end

    [tumor_cells] = updateTumorCells(tumor_cells, chem_field, params);


    % Record tumor cell metrics
    valid_tumor = ~cellfun('isempty', tumor_cells);
    if any(valid_tumor, 'all')
        oxygen_tumor = cellfun(@(tc) tc.oxygen, tumor_cells(valid_tumor));
        idx_hypoxic = oxygen_tumor <= params.ohyp;
        idx_normorxic = oxygen_tumor > params.ohyp;

        hypoxic_num(t_loc) = sum(idx_hypoxic);
        normorxic_num(t_loc) = sum(idx_normorxic);
        tumor_num(t_loc) = sum(valid_tumor);


% Adaptive dosing logic
if sum(valid_tumor) <= 0.5 * peak_tumor_num
    % Tumor burden < 50% of peak -> skip dose
    params.Sd = 0;
    dose_state = 'skip';
elseif strcmp(dose_state, 'skip') && sum(valid_tumor) > 0.8 * peak_tumor_num
    % Tumor rebounds above 50% of peak after skipping -> resume at half-dose
    params.Sd = 0.5 * 0.05;
    dose_state = 'on';
end


        dam_tumor = cellfun(@(tc) tc.damage, tumor_cells(valid_tumor));
        dam_accum(:, t_loc) = [mean(dam_tumor); std(dam_tumor)];

        thres_tumor = cellfun(@(tc) tc.death_threshold, tumor_cells(valid_tumor));
        idx_resis = thres_tumor >= 1.5 * params.base_death_threshold;
        resis_frac(t_loc) = sum(idx_resis) / sum(valid_tumor);
        death_thres(:, t_loc) = [mean(thres_tumor); std(thres_tumor)];

        oxygen_cons_tumor = cellfun(@(tc) tc.oxygen_consumption_nor, tumor_cells(valid_tumor));
        edges1 = linspace(0.5 * params.rhoo, 4 * params.rhoo, params.nbins + 1);
        oxygen_consumption{t_loc} = oxygen_cons_tumor;
        oxygen_consumption_hist(t_loc, :) = histcounts(oxygen_cons_tumor, edges1, 'Normalization', 'probability');

        prolif_rate_tumor = cellfun(@(tc) tc.proliferation_rate, tumor_cells(valid_tumor));
        edges2 = linspace(0.5 * log(2) / params.age, 4 * log(2) / params.age, params.nbins + 1);
        proliferation_rate{t_loc} = prolif_rate_tumor;
        proliferation_rate_hist(t_loc, :) = histcounts(prolif_rate_tumor, edges2, 'Normalization', 'probability');
    end

    % Save snapshots
    if mod(t_loc, plot_int) == 0
        idx = t_loc / plot_int;
        TAF_snapshots{idx} = flipud(chem_field.TAF);
        oxygen_snapshots{idx} = flipud(chem_field.oxygen);
        drug_snapshots{idx} = flipud(chem_field.drug);
        Tumor_cells{idx} = tumor_cells;
        Tip_cells{idx} = tip_cells;
        Vessel_agents{idx} = vessel_agents;
        Chem_field{idx} = chem_field;

        vessel = struct();
        vessel.tumor_cells = tumor_cells;
        vessel.tip_cells = tip_cells;
        vessel.vessel_agents = vessel_agents;
        vessel.time = t_glob;
        vessel.params = params;
        vessel_snapshots{idx} = vessel;
    end

    if t_loc == max_time
        info = struct();
        info.tumor_cells = Tumor_cells;
        info.tip_cells = Tip_cells;
        info.vessel_agents = Vessel_agents;
        info.chem_field = Chem_field;
        info.params = params;
        info.time = t_glob;
    end
end

%%
% Compute hypoxia/normoxia ratios
r_hyp = hypoxic_num ./ max(tumor_num, 1);      % avoid division by zero
r_norm = normorxic_num ./ max(tumor_num, 1);

%%
% Outputs
all_snapshots = {vessel_snapshots, TAF_snapshots, oxygen_snapshots, drug_snapshots, variables};

numerics = struct();
numerics.hypoxic_num = hypoxic_num;
numerics.normorxic_num = normorxic_num;
numerics.tumor_num = tumor_num;
numerics.dam_accum = dam_accum;
numerics.death_thres = death_thres;
numerics.oxygen_consumption = oxygen_consumption;
numerics.oxygen_consumption_hist = oxygen_consumption_hist;
numerics.proliferation_rate = proliferation_rate;
numerics.proliferation_rate_hist = proliferation_rate_hist;
numerics.resis_frac = resis_frac;
numerics.r_hyp = r_hyp;
numerics.r_norm = r_norm;

toc;
end