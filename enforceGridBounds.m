<<<<<<< HEAD
function pos = enforceGridBounds(pos, params)
% Ensure nonzero position, so that the ceiling of position in Matlab is
% positive, then we can subtract 0.5 to translate the position to the
% grid center
    pos = max(pos, 0);
    pos = min(pos, params.grid_size .* params.dg - 10000 * [eps, eps]);
=======
function pos = enforceGridBounds(pos, params)
% Ensure nonzero position, so that the ceiling of position in Matlab is
% positive, then we can subtract 0.5 to translate the position to the
% grid center
    pos = max(pos, 0);
    pos = min(pos, params.grid_size .* params.dg - 10000 * [eps, eps]);
>>>>>>> d5569c1 (Initial commit)
end