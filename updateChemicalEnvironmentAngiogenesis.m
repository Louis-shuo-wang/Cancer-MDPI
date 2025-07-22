<<<<<<< HEAD
function chem_field = updateChemicalEnvironmentAngiogenesis(tumor_cells, vessel_agents, chem_field, params)
% Update TAF concentration (equation 2)
% chem_field.endo = ADI(chem_field, params);

chem_field.TAF = updateTAF(chem_field.TAF, tumor_cells, vessel_agents, params);
end

function TAF = updateTAF(TAF, tumor_cells, vessel_agents, params)

grid_size = params.grid_size;

valid_tumor = ~cellfun('isempty', tumor_cells);
valid_vessel = ~cellfun('isempty', vessel_agents);

production = zeros(size(TAF));
uptake = zeros(size(TAF));

% Process tumor_cells (hypoxic production)
if any(valid_tumor)
    % Get positions and oxygen values from tumor cells
    pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
    pos_tumor = cat(1, pos_tumor{:});
    oxygen_tumor = cellfun(@(tc) tc.oxygen, tumor_cells(valid_tumor), 'UniformOutput', true);

    % Select hypoxic cells (oxygen <= ohyp)
    idx_hypoxic = oxygen_tumor <= params.ohyp;
    if any(idx_hypoxic)
        pos_hypoxic = pos_tumor(idx_hypoxic, :);
        pos_hypoxic = fliplr(pos_hypoxic) .* grid_size;
        hypoxic_pos_idx = floor(pos_hypoxic) + [1, 1];
        production_change = accumarray(hypoxic_pos_idx, 1, grid_size);
        production = production + params.eta * production_change;
    end
end

% Process vessel_agents (uptake)
% The position of an vessel agent is the center of a grid, then we assume
% this vessel agent will influence drug and oxygen at four corners of the
% grid: upper_left, upper_right, lower_left, lower_right
% So we respectively take flooring and ceiling of the first and second
% coordinates to get the four corner coordinates.
if any(valid_vessel)
    pos_vessel = cellfun(@(vc) vc.position, vessel_agents(valid_vessel), 'UniformOutput', false);
    pos_vessel = cat(1, pos_vessel{:});
    pos_vessel = fliplr(pos_vessel) .* grid_size;
    vessel_pos_idx = floor(pos_vessel) + [1, 1];
    uptake_change = accumarray(vessel_pos_idx, 1, grid_size);
    uptake = uptake + params.lambda * uptake_change;
end

% Update TAF using the finite difference (diffusion) scheme
rhs_TAF = -params.xic * TAF + production - uptake .* TAF;
TAF = ADImethod(TAF, rhs_TAF, params.Dc, params);
=======
function chem_field = updateChemicalEnvironmentAngiogenesis(tumor_cells, vessel_agents, chem_field, params)
% Update TAF concentration (equation 2)
% chem_field.endo = ADI(chem_field, params);

chem_field.TAF = updateTAF(chem_field.TAF, tumor_cells, vessel_agents, params);
end

function TAF = updateTAF(TAF, tumor_cells, vessel_agents, params)

grid_size = params.grid_size;

valid_tumor = ~cellfun('isempty', tumor_cells);
valid_vessel = ~cellfun('isempty', vessel_agents);

production = zeros(size(TAF));
uptake = zeros(size(TAF));

% Process tumor_cells (hypoxic production)
if any(valid_tumor)
    % Get positions and oxygen values from tumor cells
    pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
    pos_tumor = cat(1, pos_tumor{:});
    oxygen_tumor = cellfun(@(tc) tc.oxygen, tumor_cells(valid_tumor), 'UniformOutput', true);

    % Select hypoxic cells (oxygen <= ohyp)
    idx_hypoxic = oxygen_tumor <= params.ohyp;
    if any(idx_hypoxic)
        pos_hypoxic = pos_tumor(idx_hypoxic, :);
        pos_hypoxic = fliplr(pos_hypoxic) .* grid_size;
        hypoxic_pos_idx = floor(pos_hypoxic) + [1, 1];
        production_change = accumarray(hypoxic_pos_idx, 1, grid_size);
        production = production + params.eta * production_change;
    end
end

% Process vessel_agents (uptake)
% The position of an vessel agent is the center of a grid, then we assume
% this vessel agent will influence drug and oxygen at four corners of the
% grid: upper_left, upper_right, lower_left, lower_right
% So we respectively take flooring and ceiling of the first and second
% coordinates to get the four corner coordinates.
if any(valid_vessel)
    pos_vessel = cellfun(@(vc) vc.position, vessel_agents(valid_vessel), 'UniformOutput', false);
    pos_vessel = cat(1, pos_vessel{:});
    pos_vessel = fliplr(pos_vessel) .* grid_size;
    vessel_pos_idx = floor(pos_vessel) + [1, 1];
    uptake_change = accumarray(vessel_pos_idx, 1, grid_size);
    uptake = uptake + params.lambda * uptake_change;
end

% Update TAF using the finite difference (diffusion) scheme
rhs_TAF = -params.xic * TAF + production - uptake .* TAF;
TAF = ADImethod(TAF, rhs_TAF, params.Dc, params);
>>>>>>> d5569c1 (Initial commit)
end