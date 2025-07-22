<<<<<<< HEAD
% performProliferation.m
function [daughter1, daughter2] = performProliferation(parent, chem_field, params)

% Find available space for second daughter
new_pos = parent.position + 0.1 * [cos(2 * pi * rand()), sin(2 * pi * rand())] .* params.dg;

% Create base daughter cell properties
daughter1 = struct();
daughter2 = struct();
daughter1.oxygen_consumption_nor = parent.oxygen_consumption_nor;
daughter2.oxygen_consumption_nor = parent.oxygen_consumption_nor;
daughter1.death_threshold = parent.death_threshold;
daughter2.death_threshold = parent.death_threshold;
daughter1.proliferation_rate = parent.proliferation_rate;
daughter2.proliferation_rate = parent.proliferation_rate;
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
=======
% performProliferation.m
function [daughter1, daughter2] = performProliferation(parent, chem_field, params)

% Find available space for second daughter
new_pos = parent.position + 0.1 * [cos(2 * pi * rand()), sin(2 * pi * rand())] .* params.dg;

% Create base daughter cell properties
daughter1 = struct();
daughter2 = struct();
daughter1.oxygen_consumption_nor = parent.oxygen_consumption_nor;
daughter2.oxygen_consumption_nor = parent.oxygen_consumption_nor;
daughter1.death_threshold = parent.death_threshold;
daughter2.death_threshold = parent.death_threshold;
daughter1.proliferation_rate = parent.proliferation_rate;
daughter2.proliferation_rate = parent.proliferation_rate;
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
>>>>>>> d5569c1 (Initial commit)
end