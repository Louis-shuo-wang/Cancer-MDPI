% performProliferation.m
function [daughter1, daughter2] = performProliferation(parent, chem_field, params)

% Find available space for second daughter
new_pos = parent.position + 0.1 * [cos(2 * pi * rand()), sin(2 * pi * rand())] .* params.dg;

% Create base daughter cell properties
daughter1 = struct();
daughter2 = struct();
daughters = {daughter1, daughter2};
daughters = handleMutation(daughters, parent, params);
daughter1 = daughters{1};
daughter2 = daughters{2};
daughter1.maturation_time = log(2) / daughter1.proliferation_rate;
daughter2.maturation_time = log(2) / daughter2.proliferation_rate;

% Position
daughter1.position = parent.position;
daughter2.position = new_pos;

% Cellular oxygen contents hypoxic cell & hypoxic cells consume less oxygen
daughter1.oxygen = updateCellularOxygenUptake(daughter1, chem_field, params);
if daughter1.oxygen <= params.oapop
    daughter1 = [];
elseif daughter1.oxygen <= params.ohyp
    daughter1.type = 'hypoxic';
    daughter1.oxygen_consumption = daughter1.oxygen_consumption_nor / 2;
else
    daughter1.type = 'normorxic';
    daughter1.oxygen_consumption = daughter1.oxygen_consumption_nor;
end

daughter2.oxygen = updateCellularOxygenUptake(daughter2, chem_field, params);
if daughter2.oxygen <= params.oapop
    daughter2 = [];
elseif daughter2.oxygen <= params.ohyp
    daughter2.type = 'hypoxic';
    daughter2.oxygen_consumption = daughter2.oxygen_consumption_nor/2;
else
    daughter2.type = 'normorxic';
    daughter2.oxygen_consumption = daughter2.oxygen_consumption_nor;
end

% Reset age and exposure time
daughter1.age = 0;
daughter2.age = 0;
daughter1.exposure_time = 0;
daughter2.exposure_time = 0;

% Split accumulated drug and damage
daughter1.drug = parent.drug / 2;
daughter2.drug = parent.drug / 2;
daughter1.damage = parent.damage / 2;
daughter2.damage = parent.damage / 2;

% Label
daughter1.label = [parent.label, 1];
daughter2.label = [parent.label, 2];
end

% handleMutation.m
function daughters = handleMutation(daughters, parent, params)
for j = 1:2
    if rand() <= params.mutation_rate * params.dt
        daughters{j} = applyMutation(daughters{j}, parent, params);
    else
        daughters{j}.oxygen_consumption_nor = parent.oxygen_consumption_nor;
        daughters{j}.death_threshold = parent.death_threshold;
        daughters{j}.proliferation_rate = parent.proliferation_rate;
    end
end
end

% applyLinearMutation.m
    function daughter_cell = applyMutation(daughter_cell, parent_cell, params)

        oxygen_factor = 0.7 + 1 * rand();
        threshold_factor = 0.7 + 1 * rand();
        proliferation_factor = 0.7 + 1 * rand();

        % Apply mutations with constraints
        daughter_cell.oxygen_consumption_nor = constrainValue(parent_cell.oxygen_consumption_nor * oxygen_factor, ...
            0.5 * params.rhoo, 4 * params.rhoo);
        daughter_cell.death_threshold = constrainValue(parent_cell.death_threshold * threshold_factor, ...
            0.5 * params.base_death_threshold, 4 * params.base_death_threshold);
        daughter_cell.proliferation_rate = constrainValue(parent_cell.proliferation_rate * proliferation_factor, ...
            0.5 * log(2) / params.age, 4 * log(2) / params.age);
    end

    function value = constrainValue(value, min_val, max_val)
        value = min(max(value, min_val), max_val);
    end
