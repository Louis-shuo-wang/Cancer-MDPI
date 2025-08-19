function [tumor_cells] = updateTumorCells(tumor_cells, chem_field, params)

% processTumorCells.m
% Preallocate a cell array to store the updated tumor cells
updated_tumor_cells = cell(size(tumor_cells));

% Preallocate daughter_cells
daughter_cells = cell(length(tumor_cells), 2);

parfor i = 1:length(tumor_cells)

    temp_cell = tumor_cells{i};   % Create a temporary variable to store the current cell
    daughter_cells_temp = cell(1, 2);

    if isempty(temp_cell)
        updated_tumor_cells{i} = temp_cell;
        % Ensure consistent assignment for daughter_cells for this iteration
        daughter_cells_temp = {[], []};
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
end
