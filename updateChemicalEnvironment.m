<<<<<<< HEAD
% updateChemicalEnvironment.m
function chem_field = updateChemicalEnvironment(tumor_cells, vessel_agents, chem_field, params)

% Update drug concentration (equation 4)
chem_field.drug = updateDrug(chem_field.drug, tumor_cells, vessel_agents, params);

% Update oxygen concentration (equation 5)
chem_field.oxygen = updateOxygen(chem_field.oxygen, tumor_cells, vessel_agents, params);
end

% updateDrug.m
function drug = updateDrug(drug, tumor_cells, vessel_agents, params)

grid_size = params.grid_size;

valid_tumor = ~cellfun('isempty', tumor_cells);
valid_vessel = ~cellfun('isempty', vessel_agents);

% Initialize source terms
uptake = zeros(size(drug));
supply = zeros(size(drug));

% Add drug uptake by tumor cells
if any(valid_tumor)
    pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
    pos_tumor = cat(1, pos_tumor{:});
    pos_tumor = fliplr(pos_tumor) .* grid_size;
    tumor_pos_idx = floor(pos_tumor) + [1, 1];
    uptake_change = accumarray(tumor_pos_idx, 1, grid_size);
    uptake = uptake + params.rhod * uptake_change;
end

% Add drug supply from vessels
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
    supply_change = accumarray(vessel_pos_idx, 1, grid_size);
    supply = supply + params.Sd * supply_change;
end

% Update using finite difference method
rhs_drug = - params.xid * drug - uptake .* drug + supply;
drug = ADImethod(drug, rhs_drug, params.Dd, params);
drug = max(0, drug);
end

% updateOxygen.m
function oxygen = updateOxygen(oxygen, tumor_cells, vessel_agents, params)

grid_size = params.grid_size;

valid_tumor = ~cellfun('isempty', tumor_cells);
valid_vessel = ~cellfun('isempty', vessel_agents);

% Initialize source terms
uptake = zeros(size(oxygen));
supply = zeros(size(oxygen));

% Add oxygen uptake by tumor cells
if any(valid_tumor)
    pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
    pos_tumor = cat(1, pos_tumor{:});
    pos_tumor = fliplr(pos_tumor) .* grid_size;
    tumor_pos_idx = floor(pos_tumor) + [1, 1];
    oxygen_cons = cellfun(@(tc) tc.oxygen_consumption, tumor_cells(valid_tumor), 'UniformOutput', true);
    uptake_change = accumarray(tumor_pos_idx, oxygen_cons, grid_size, [], [], true);
    uptake = uptake + uptake_change;
end

% Add oxygen supply from vessels
if any(valid_vessel)
    pos_vessel = cellfun(@(vc) vc.position, vessel_agents(valid_vessel), 'UniformOutput', false);
    pos_vessel = cat(1, pos_vessel{:});
    pos_vessel = fliplr(pos_vessel) .* grid_size;
    vessel_pos_idx = floor(pos_vessel) + [1, 1];
    supply_change = accumarray(vessel_pos_idx, 1, grid_size);
    supply = supply + supply_change;
end

% Update using finite difference method
rhs_oxygen = - params.xio * oxygen - uptake + params.So * (params.omax - oxygen) .* supply;
oxygen = ADImethod(oxygen, rhs_oxygen, params.Do, params);
oxygen = max(0, oxygen);
=======
% updateChemicalEnvironment.m
function chem_field = updateChemicalEnvironment(tumor_cells, vessel_agents, chem_field, params)

% Update drug concentration (equation 4)
chem_field.drug = updateDrug(chem_field.drug, tumor_cells, vessel_agents, params);

% Update oxygen concentration (equation 5)
chem_field.oxygen = updateOxygen(chem_field.oxygen, tumor_cells, vessel_agents, params);
end

% updateDrug.m
function drug = updateDrug(drug, tumor_cells, vessel_agents, params)

grid_size = params.grid_size;

valid_tumor = ~cellfun('isempty', tumor_cells);
valid_vessel = ~cellfun('isempty', vessel_agents);

% Initialize source terms
uptake = zeros(size(drug));
supply = zeros(size(drug));

% Add drug uptake by tumor cells
if any(valid_tumor)
    pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
    pos_tumor = cat(1, pos_tumor{:});
    pos_tumor = fliplr(pos_tumor) .* grid_size;
    tumor_pos_idx = floor(pos_tumor) + [1, 1];
    uptake_change = accumarray(tumor_pos_idx, 1, grid_size);
    uptake = uptake + params.rhod * uptake_change;
end

% Add drug supply from vessels
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
    supply_change = accumarray(vessel_pos_idx, 1, grid_size);
    supply = supply + params.Sd * supply_change;
end

% Update using finite difference method
rhs_drug = - params.xid * drug - uptake .* drug + supply;
drug = ADImethod(drug, rhs_drug, params.Dd, params);
drug = max(0, drug);
end

% updateOxygen.m
function oxygen = updateOxygen(oxygen, tumor_cells, vessel_agents, params)

grid_size = params.grid_size;

valid_tumor = ~cellfun('isempty', tumor_cells);
valid_vessel = ~cellfun('isempty', vessel_agents);

% Initialize source terms
uptake = zeros(size(oxygen));
supply = zeros(size(oxygen));

% Add oxygen uptake by tumor cells
if any(valid_tumor)
    pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
    pos_tumor = cat(1, pos_tumor{:});
    pos_tumor = fliplr(pos_tumor) .* grid_size;
    tumor_pos_idx = floor(pos_tumor) + [1, 1];
    oxygen_cons = cellfun(@(tc) tc.oxygen_consumption, tumor_cells(valid_tumor), 'UniformOutput', true);
    uptake_change = accumarray(tumor_pos_idx, oxygen_cons, grid_size, [], [], true);
    uptake = uptake + uptake_change;
end

% Add oxygen supply from vessels
if any(valid_vessel)
    pos_vessel = cellfun(@(vc) vc.position, vessel_agents(valid_vessel), 'UniformOutput', false);
    pos_vessel = cat(1, pos_vessel{:});
    pos_vessel = fliplr(pos_vessel) .* grid_size;
    vessel_pos_idx = floor(pos_vessel) + [1, 1];
    supply_change = accumarray(vessel_pos_idx, 1, grid_size);
    supply = supply + supply_change;
end

% Update using finite difference method
rhs_oxygen = - params.xio * oxygen - uptake + params.So * (params.omax - oxygen) .* supply;
oxygen = ADImethod(oxygen, rhs_oxygen, params.Do, params);
oxygen = max(0, oxygen);
>>>>>>> d5569c1 (Initial commit)
end