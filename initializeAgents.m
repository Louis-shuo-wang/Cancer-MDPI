<<<<<<< HEAD
% initializeAgents.m
function [tumor_cells, tip_cells, vessel_agents, chem_field] = initializeAgents(params)

% Initialize empty arrays for agents
tumor_cells = cell(1, 500 * params.initial_tumor_number);
tip_cells = cell(1, 50 * params.initial_tip_number);
vessel_agents = cell(1, params.grid_size(1) * params.grid_size(2));

% Initialize chemical fields
chem_field = struct();
% chem_field.oxygen = params.initial_oxygen_concentration * ones(params.grid_size);
% chem_field.drug = params.initial_drug_concentration * ones(params.grid_size);
[X, Y] = meshgrid(params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2, ...
                  params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2);
% chem_field.endo = exp(- Y ^ 2 / 0.001) * sin(12 * pi * X) .^ 2;
% % chem_field.TAF = 5 * Y;
chem_field.oxygen = 4*cos(2*pi*Y)*cos(2*pi*X);
chem_field.drug   = cos(4*pi*Y)*cos(4*pi*X);
chem_field.TAF    = cos(2*pi*Y)*cos(2*pi*X);
% chem_field.TAF = 1 - sqrt((X - params.tumor_center(1)) .^ 2 + (Y - params.tumor_center(2)) .^ 2);

new_cell = struct(...
    'oxygen', params.initial_oxygen_concentration, ...
    'drug', 0, ...
    'damage', 0, ...
    'death_threshold', params.base_death_threshold, ...
    'exposure_time', 0, ...
    'type', 'normoxic', ...
    'oxygen_consumption_nor', params.rhoo, ...
    'oxygen_consumption', params.rhoo);

% Add initial tumor cells
for i = 1:params.initial_tumor_number

    time_mat = params.age;
    rate_prol = log(2) / time_mat;
    tumor_cells{i} = new_cell;
    tumor_cells{i}.maturation_time = time_mat;
    tumor_cells{i}.age = tumor_cells{i}.maturation_time * rand();
    tumor_cells{i}.proliferation_rate = rate_prol;
    tumor_cells{i}.position = (params.tumor_center + 0.2 * rand(1, 2) .* [cos(2 * pi * rand()), sin(2 * pi * rand())]);
    tumor_cells{i}.label = i;
end

new_tip = struct(...
    'age', params.psi, ...
    'proliferation_time', 0);

pos = {[0.25, 0.25], ...
       [0.25, 0.75], ...
       [0.75, 0.25], ...
       [0.75, 0.75]};

% Add initial tips
for i = 1:params.initial_tip_number

    % Create new tip cell
    tip_cells{i} = new_tip;
    tip_cells{i}.position = pos{i};
    tip_cells{i}.order = i;

    % Add initial vessel agent
    pos_idx = floor(tip_cells{i}.position .* params.grid_size) + [1, 1];
    ind = sub2ind(params.grid_size, pos_idx(2), pos_idx(1));
    vessel_agents{ind} = struct('position', tip_cells{i}.position, ...
                                'order', i);
end
=======
% initializeAgents.m
function [tumor_cells, tip_cells, vessel_agents, chem_field] = initializeAgents(params)

% Initialize empty arrays for agents
tumor_cells = cell(1, 500 * params.initial_tumor_number);
tip_cells = cell(1, 50 * params.initial_tip_number);
vessel_agents = cell(1, params.grid_size(1) * params.grid_size(2));

% Initialize chemical fields
chem_field = struct();
% chem_field.oxygen = params.initial_oxygen_concentration * ones(params.grid_size);
% chem_field.drug = params.initial_drug_concentration * ones(params.grid_size);
[X, Y] = meshgrid(params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2, ...
                  params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2);
% chem_field.endo = exp(- Y ^ 2 / 0.001) * sin(12 * pi * X) .^ 2;
% % chem_field.TAF = 5 * Y;
chem_field.oxygen = 4*cos(2*pi*Y)*cos(2*pi*X);
chem_field.drug   = cos(4*pi*Y)*cos(4*pi*X);
chem_field.TAF    = cos(2*pi*Y)*cos(2*pi*X);
% chem_field.TAF = 1 - sqrt((X - params.tumor_center(1)) .^ 2 + (Y - params.tumor_center(2)) .^ 2);

new_cell = struct(...
    'oxygen', params.initial_oxygen_concentration, ...
    'drug', 0, ...
    'damage', 0, ...
    'death_threshold', params.base_death_threshold, ...
    'exposure_time', 0, ...
    'type', 'normoxic', ...
    'oxygen_consumption_nor', params.rhoo, ...
    'oxygen_consumption', params.rhoo);

% Add initial tumor cells
for i = 1:params.initial_tumor_number

    time_mat = params.age;
    rate_prol = log(2) / time_mat;
    tumor_cells{i} = new_cell;
    tumor_cells{i}.maturation_time = time_mat;
    tumor_cells{i}.age = tumor_cells{i}.maturation_time * rand();
    tumor_cells{i}.proliferation_rate = rate_prol;
    tumor_cells{i}.position = (params.tumor_center + 0.2 * rand(1, 2) .* [cos(2 * pi * rand()), sin(2 * pi * rand())]);
    tumor_cells{i}.label = i;
end

new_tip = struct(...
    'age', params.psi, ...
    'proliferation_time', 0);

pos = {[0.25, 0.25], ...
       [0.25, 0.75], ...
       [0.75, 0.25], ...
       [0.75, 0.75]};

% Add initial tips
for i = 1:params.initial_tip_number

    % Create new tip cell
    tip_cells{i} = new_tip;
    tip_cells{i}.position = pos{i};
    tip_cells{i}.order = i;

    % Add initial vessel agent
    pos_idx = floor(tip_cells{i}.position .* params.grid_size) + [1, 1];
    ind = sub2ind(params.grid_size, pos_idx(2), pos_idx(1));
    vessel_agents{ind} = struct('position', tip_cells{i}.position, ...
                                'order', i);
end
>>>>>>> d5569c1 (Initial commit)
end