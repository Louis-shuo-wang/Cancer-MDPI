<<<<<<< HEAD
function [tumor_cells, div_log] = updateTumorCells(tumor_cells, chem_field, params, div_log)

div_log_array = cell(1, 2* length(tumor_cells));  % 每个worker的局部日志

% processTumorCells.m
% Preallocate a cell array to store the updated tumor cells
updated_tumor_cells = cell(size(tumor_cells));

% Preallocate daughter_cells
daughter_cells = cell(length(tumor_cells), 2);

parfor i = 1:length(tumor_cells)
    div_log_temp = initDivLogEntry();

    temp_cell = tumor_cells{i};   % Create a temporary variable to store the current cell
    daughter_cells_temp = cell(1, 2);

    if isempty(temp_cell)
        updated_tumor_cells{i} = temp_cell;
        % Ensure consistent assignment for daughter_cells for this iteration
        daughter_cells_temp = {[], []};
    end

    % handle mutation
    if ~isempty(temp_cell)
        temp_cell = handleMutation(temp_cell,params)
    end

    if ~isempty(temp_cell)
        % Update cellular drug and oxygen uptake (pro2)
        temp_cell. oxygen = updateCellularOxygenUptake(temp_cell, chem_field, params)
        
        % Check oxygen levels and handle cell fate
        if temp_cell.oxygen <= params.oapop
            % Remove cell (apoptosis) (pro4a)
            temp_cell = [];
        end

        if ~isempty(temp_cell)
            % handle mutation
            temp_cell = handleMutation(temp_cell,params)

            if temp_cell.oxygen <= params.ohyp
                % Cell is hypoxic - no age increase
                temp_cell.type = 'hypoxic';
                temp_cell.oxygen_consumption = temp_cell.oxygen_consumption_nor / 2;
            else
                % Cell is normoxic - increase age (pro1a)
                temp_cell.type = 'normoxic';
                temp_cell.oxygen_consumption = temp_cell.oxygen_consumption_nor;
                temp_cell.age = temp_cell.age + params.dt;
            end

            % Check for cell maturity and proliferation
            if strcmp(temp_cell.type, 'normoxic') && ...
                    temp_cell.age >= temp_cell.maturation_time

                % Handle proliferation if space available
                if hasAvailableSpace(temp_cell, tumor_cells, params)
                    [daughter1, daughter2] = performProliferation(temp_cell, chem_field, params);
                    div_log_temp.pre_damage(end+1)     = temp_cell.damage;
                    div_log_temp.death_thresh(end+1)   = temp_cell.death_threshold;
                    div_log_temp.post_damage1(end+1)   = daughter1.damage;
                    div_log_temp.post_damage2(end+1)   = daughter2.damage;
                    temp_cell = [];   % Delete parent cell after proliferation
                    daughter_cells_temp = {daughter1, daughter2};
                end
            end

            % Update damage and resistance (pro5, pro6)
            if ~isempty(temp_cell)
                temp_cell = updateDamageAndResistance(temp_cell, chem_field, params);
                % Check cell survival
                if temp_cell.damage > temp_cell.death_threshold
                    temp_cell = [];
                end
            end
        end
        updated_tumor_cells{i} = temp_cell;      % Assign the updated cell to the preallocated cell array
    end
    % Assign the entire row for iteration i with a consistent indexing pattern
    daughter_cells(i, :) = daughter_cells_temp;
    div_log_array{i} = div_log_temp;
end

div_log.pre_damage   = [];
div_log.death_thresh = [];
div_log.post_damage1 = [];
div_log.post_damage2 = [];

for i = 1:length(div_log_array)
    if isempty(div_log_array{i})
        continue;
    end
    div_log.pre_damage   = [div_log.pre_damage,   div_log_array{i}.pre_damage];
    div_log.death_thresh = [div_log.death_thresh, div_log_array{i}.death_thresh];
    div_log.post_damage1 = [div_log.post_damage1, div_log_array{i}.post_damage1];
    div_log.post_damage2 = [div_log.post_damage2, div_log_array{i}.post_damage2];
end

% Combine the results from all workers
for j = 1:length(tumor_cells)
    tumor_cells{j} = updated_tumor_cells{j};
end

% Add new daughter cells to tumor_cells after the parfor loop
for k = 1:size(daughter_cells,1)
    for kk = 1:2
        if ~isempty(daughter_cells{k, kk})
            empty_tumor = cellfun('isempty', tumor_cells);
            idx_tumor = find(empty_tumor, 1);
            if ~isempty(idx_tumor)
                tumor_cells{idx_tumor} = daughter_cells{k, kk};
            else
                tumor_cells{end + 1} = daughter_cells{k, kk};
            end
        end
    end
end
end

function log = initDivLogEntry()
log.pre_damage = [];
log.death_thresh = [];
log.post_damage1 = [];
log.post_damage2 = [];
end

% handleMutation.m
function cell = handleMutation(cell, params)
if rand() <= params.mutation_rate * params.dt
    cell = applyMutation(cell, params);
else
    cell.oxygen_consumption_nor = cell.oxygen_consumption_nor;
    cell.death_threshold = cell.death_threshold;
    cell.proliferation_rate = cell.proliferation_rate;
end
end

% applyLinearMutation.m
function cell = applyMutation(cell, params)

oxygen_factor = 0.7 + 1 * rand();
threshold_factor = 0.7 + 1 * rand();
proliferation_factor = 0.7 + 1 * rand();

% Apply mutations with constraints
cell.oxygen_consumption_nor = constrainValue(cell.oxygen_consumption_nor * oxygen_factor, ...
    0.5 * params.rhoo, 4 * params.rhoo);
cell.death_threshold = constrainValue(cell.death_threshold * threshold_factor, ...
    0.5 * params.base_death_threshold, 4 * params.base_death_threshold);
cell.proliferation_rate = constrainValue(cell.proliferation_rate * proliferation_factor, ...
    0.5 * log(2) / params.age, 4 * log(2) / params.age);
end

function value = constrainValue(value, min_val, max_val)
value = min(max(value, min_val), max_val);
end

function available = hasAvailableSpace(cell, tumor_cells, params)
% Check neighboring positions for available space
% Count the number of neighbouring tumor cells
valid_tumor = ~cellfun('isempty', tumor_cells);
% valid_vessel = ~cellfun('isempty', vessel_agents);
pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
pos_tumor = cat(1, pos_tumor{:});
nbhd_number = sum(vecnorm(cell.position - pos_tumor, 2, 2) <= params.Rc);

if nbhd_number < params.max_proliferation_density + 1 % the calculation of nbhd number involving the current tumor cell
    available = true;
else
    available = false;
end
end

function tumor_cell = updateDamageAndResistance(tumor_cell, chem_field, params)
% Update drug-induced damage
position = tumor_cell.position;
drug_level = interp2(params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2, ...
    (params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2)', ...
    chem_field.drug, ...
    position(1), position(2), 'spline');
damage_change = drug_level;
tumor_cell.damage = tumor_cell.damage + (damage_change - params.pr * tumor_cell.damage) * params.dt;
tumor_cell.drug = tumor_cell.drug + damage_change * params.dt;

% Handle acquired resistance
if tumor_cell.exposure_time > params.texp
    tumor_cell.death_threshold = tumor_cell.death_threshold + params.Deltadeath * params.dt;
end

% Update drug exposure time
if tumor_cell.drug > params.dexp
    tumor_cell.exposure_time = tumor_cell.exposure_time + params.dt;
end
end
=======
function [tumor_cells, div_log] = updateTumorCells(tumor_cells, chem_field, params, div_log)

div_log_array = cell(1, 2* length(tumor_cells));  % 每个worker的局部日志

% processTumorCells.m
% Preallocate a cell array to store the updated tumor cells
updated_tumor_cells = cell(size(tumor_cells));

% Preallocate daughter_cells
daughter_cells = cell(length(tumor_cells), 2);

parfor i = 1:length(tumor_cells)
    div_log_temp = initDivLogEntry();

    temp_cell = tumor_cells{i};   % Create a temporary variable to store the current cell
    daughter_cells_temp = cell(1, 2);

    if isempty(temp_cell)
        updated_tumor_cells{i} = temp_cell;
        % Ensure consistent assignment for daughter_cells for this iteration
        daughter_cells_temp = {[], []};
    end

    % handle mutation
    if ~isempty(temp_cell)
        temp_cell = handleMutation(temp_cell,params)
    end

    if ~isempty(temp_cell)
        % Update cellular drug and oxygen uptake (pro2)
        temp_cell. oxygen = updateCellularOxygenUptake(temp_cell, chem_field, params)
        
        % Check oxygen levels and handle cell fate
        if temp_cell.oxygen <= params.oapop
            % Remove cell (apoptosis) (pro4a)
            temp_cell = [];
        end

        if ~isempty(temp_cell)
            % handle mutation
            temp_cell = handleMutation(temp_cell,params)

            if temp_cell.oxygen <= params.ohyp
                % Cell is hypoxic - no age increase
                temp_cell.type = 'hypoxic';
                temp_cell.oxygen_consumption = temp_cell.oxygen_consumption_nor / 2;
            else
                % Cell is normoxic - increase age (pro1a)
                temp_cell.type = 'normoxic';
                temp_cell.oxygen_consumption = temp_cell.oxygen_consumption_nor;
                temp_cell.age = temp_cell.age + params.dt;
            end

            % Check for cell maturity and proliferation
            if strcmp(temp_cell.type, 'normoxic') && ...
                    temp_cell.age >= temp_cell.maturation_time

                % Handle proliferation if space available
                if hasAvailableSpace(temp_cell, tumor_cells, params)
                    [daughter1, daughter2] = performProliferation(temp_cell, chem_field, params);
                    div_log_temp.pre_damage(end+1)     = temp_cell.damage;
                    div_log_temp.death_thresh(end+1)   = temp_cell.death_threshold;
                    div_log_temp.post_damage1(end+1)   = daughter1.damage;
                    div_log_temp.post_damage2(end+1)   = daughter2.damage;
                    temp_cell = [];   % Delete parent cell after proliferation
                    daughter_cells_temp = {daughter1, daughter2};
                end
            end

            % Update damage and resistance (pro5, pro6)
            if ~isempty(temp_cell)
                temp_cell = updateDamageAndResistance(temp_cell, chem_field, params);
                % Check cell survival
                if temp_cell.damage > temp_cell.death_threshold
                    temp_cell = [];
                end
            end
        end
        updated_tumor_cells{i} = temp_cell;      % Assign the updated cell to the preallocated cell array
    end
    % Assign the entire row for iteration i with a consistent indexing pattern
    daughter_cells(i, :) = daughter_cells_temp;
    div_log_array{i} = div_log_temp;
end

div_log.pre_damage   = [];
div_log.death_thresh = [];
div_log.post_damage1 = [];
div_log.post_damage2 = [];

for i = 1:length(div_log_array)
    if isempty(div_log_array{i})
        continue;
    end
    div_log.pre_damage   = [div_log.pre_damage,   div_log_array{i}.pre_damage];
    div_log.death_thresh = [div_log.death_thresh, div_log_array{i}.death_thresh];
    div_log.post_damage1 = [div_log.post_damage1, div_log_array{i}.post_damage1];
    div_log.post_damage2 = [div_log.post_damage2, div_log_array{i}.post_damage2];
end

% Combine the results from all workers
for j = 1:length(tumor_cells)
    tumor_cells{j} = updated_tumor_cells{j};
end

% Add new daughter cells to tumor_cells after the parfor loop
for k = 1:size(daughter_cells,1)
    for kk = 1:2
        if ~isempty(daughter_cells{k, kk})
            empty_tumor = cellfun('isempty', tumor_cells);
            idx_tumor = find(empty_tumor, 1);
            if ~isempty(idx_tumor)
                tumor_cells{idx_tumor} = daughter_cells{k, kk};
            else
                tumor_cells{end + 1} = daughter_cells{k, kk};
            end
        end
    end
end
end

function log = initDivLogEntry()
log.pre_damage = [];
log.death_thresh = [];
log.post_damage1 = [];
log.post_damage2 = [];
end

% handleMutation.m
function cell = handleMutation(cell, params)
if rand() <= params.mutation_rate * params.dt
    cell = applyMutation(cell, params);
else
    cell.oxygen_consumption_nor = cell.oxygen_consumption_nor;
    cell.death_threshold = cell.death_threshold;
    cell.proliferation_rate = cell.proliferation_rate;
end
end

% applyLinearMutation.m
function cell = applyMutation(cell, params)

oxygen_factor = 0.7 + 1 * rand();
threshold_factor = 0.7 + 1 * rand();
proliferation_factor = 0.7 + 1 * rand();

% Apply mutations with constraints
cell.oxygen_consumption_nor = constrainValue(cell.oxygen_consumption_nor * oxygen_factor, ...
    0.5 * params.rhoo, 4 * params.rhoo);
cell.death_threshold = constrainValue(cell.death_threshold * threshold_factor, ...
    0.5 * params.base_death_threshold, 4 * params.base_death_threshold);
cell.proliferation_rate = constrainValue(cell.proliferation_rate * proliferation_factor, ...
    0.5 * log(2) / params.age, 4 * log(2) / params.age);
end

function value = constrainValue(value, min_val, max_val)
value = min(max(value, min_val), max_val);
end

function available = hasAvailableSpace(cell, tumor_cells, params)
% Check neighboring positions for available space
% Count the number of neighbouring tumor cells
valid_tumor = ~cellfun('isempty', tumor_cells);
% valid_vessel = ~cellfun('isempty', vessel_agents);
pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
pos_tumor = cat(1, pos_tumor{:});
nbhd_number = sum(vecnorm(cell.position - pos_tumor, 2, 2) <= params.Rc);

if nbhd_number < params.max_proliferation_density + 1 % the calculation of nbhd number involving the current tumor cell
    available = true;
else
    available = false;
end
end

function tumor_cell = updateDamageAndResistance(tumor_cell, chem_field, params)
% Update drug-induced damage
position = tumor_cell.position;
drug_level = interp2(params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2, ...
    (params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2)', ...
    chem_field.drug, ...
    position(1), position(2), 'spline');
damage_change = drug_level;
tumor_cell.damage = tumor_cell.damage + (damage_change - params.pr * tumor_cell.damage) * params.dt;
tumor_cell.drug = tumor_cell.drug + damage_change * params.dt;

% Handle acquired resistance
if tumor_cell.exposure_time > params.texp
    tumor_cell.death_threshold = tumor_cell.death_threshold + params.Deltadeath * params.dt;
end

% Update drug exposure time
if tumor_cell.drug > params.dexp
    tumor_cell.exposure_time = tumor_cell.exposure_time + params.dt;
end
end
>>>>>>> d5569c1 (Initial commit)
