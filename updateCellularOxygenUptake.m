<<<<<<< HEAD
% updateCellularUptake.m
function oxygen = updateCellularOxygenUptake(tumor_cell, chem_field, params)
% Add oxygen uptake by tumor cells
position = full(tumor_cell.position);
oxygen_level = interp2(params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2, ...
    (params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2)', ...
    chem_field.oxygen, position(1), position(2), 'spline');
oxygen = oxygen_level;
=======
% updateCellularUptake.m
function oxygen = updateCellularOxygenUptake(tumor_cell, chem_field, params)
% Add oxygen uptake by tumor cells
position = full(tumor_cell.position);
oxygen_level = interp2(params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2, ...
    (params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2)', ...
    chem_field.oxygen, position(1), position(2), 'spline');
oxygen = oxygen_level;
>>>>>>> d5569c1 (Initial commit)
end