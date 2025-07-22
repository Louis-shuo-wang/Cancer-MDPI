<<<<<<< HEAD
% visualizeSimulation.m
function plotTumorAndVessels(tumor_cells, vessel_agents, params)
hold on;

% Preallocate marker handles as empty arrays
h_vessel = [];
h_hypoxic = [];
h_normorxic = [];

% Treatment of vessel_agents cell may be wrong
valid_vessel = ~cellfun('isempty', vessel_agents);
valid_tumor = ~cellfun('isempty', tumor_cells);

% Plot vessels
if any(valid_vessel)
    pos_vessel = cellfun(@(vc) vc.position, vessel_agents(valid_vessel), 'UniformOutput', false);
    pos_vessel = cat(1, pos_vessel{:});

    % Plot each vessel as a filled rectangle
        % Create rectangle with vessel cell position as center
        % (x,y) is the bottom-left corner of the rectangle
        % Plot a dummy vessel (or save the first rectangle handle)
        h_vessel = scatter(pos_vessel(:, 1), pos_vessel(:, 2), 0.05, 'blue');
    % Create a dummy object (invisible plot) that will be used in the legend.
    % h_vessel = plot(nan, nan, 's', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
end

% Plot tumor cells

if any(valid_tumor)
    pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
    pos_tumor = cat(1, pos_tumor{:});
    oxygen_tumor = cellfun(@(tc) tc.oxygen, tumor_cells(valid_tumor), 'UniformOutput', true);
    % Select hypoxic cells (oxygen <= ohyp)
    idx_hypoxic = oxygen_tumor <= params.ohyp;
    if any(idx_hypoxic)
        pos_hypoxic = pos_tumor(idx_hypoxic, :);
        h_hypoxic = scatter(pos_hypoxic(:, 1), pos_hypoxic(:, 2), 4, 'magenta');
    end

    idx_normorxic = oxygen_tumor > params.ohyp;
    if any(idx_normorxic)
        pos_normorxic = pos_tumor(idx_normorxic, :);
        h_normorxic = scatter(pos_normorxic(:, 1), pos_normorxic(:, 2), 4, 'red');
    end
end


% Build legend handles conditionally
legend_handles = [];
legend_labels = {};

if ~isempty(h_vessel)
    legend_handles = [legend_handles, h_vessel];
    legend_labels = [legend_labels, {'Vessel'}];
end
if exist('h_hypoxic', 'var') && ~isempty(h_hypoxic) && isgraphics(h_hypoxic)
    legend_handles = [legend_handles, h_hypoxic];
    legend_labels = [legend_labels, {'Hypoxic'}];
end
if exist('h_normorxic', 'var') && ~isempty(h_normorxic) && isgraphics(h_normorxic)
    legend_handles = [legend_handles, h_normorxic];
    legend_labels = [legend_labels, {'Normorxic'}];
end

if ~isempty(legend_handles)
    legend(legend_handles, legend_labels, 'FontSize', 6, 'Location', 'best');
end

% title('Tumor Cells and Vessels');
hold off;

=======
% visualizeSimulation.m
function plotTumorAndVessels(tumor_cells, vessel_agents, params)
hold on;

% Preallocate marker handles as empty arrays
h_vessel = [];
h_hypoxic = [];
h_normorxic = [];

% Treatment of vessel_agents cell may be wrong
valid_vessel = ~cellfun('isempty', vessel_agents);
valid_tumor = ~cellfun('isempty', tumor_cells);

% Plot vessels
if any(valid_vessel)
    pos_vessel = cellfun(@(vc) vc.position, vessel_agents(valid_vessel), 'UniformOutput', false);
    pos_vessel = cat(1, pos_vessel{:});

    % Plot each vessel as a filled rectangle
        % Create rectangle with vessel cell position as center
        % (x,y) is the bottom-left corner of the rectangle
        % Plot a dummy vessel (or save the first rectangle handle)
        h_vessel = scatter(pos_vessel(:, 1), pos_vessel(:, 2), 0.05, 'blue');
    % Create a dummy object (invisible plot) that will be used in the legend.
    % h_vessel = plot(nan, nan, 's', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
end

% Plot tumor cells

if any(valid_tumor)
    pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
    pos_tumor = cat(1, pos_tumor{:});
    oxygen_tumor = cellfun(@(tc) tc.oxygen, tumor_cells(valid_tumor), 'UniformOutput', true);
    % Select hypoxic cells (oxygen <= ohyp)
    idx_hypoxic = oxygen_tumor <= params.ohyp;
    if any(idx_hypoxic)
        pos_hypoxic = pos_tumor(idx_hypoxic, :);
        h_hypoxic = scatter(pos_hypoxic(:, 1), pos_hypoxic(:, 2), 4, 'magenta');
    end

    idx_normorxic = oxygen_tumor > params.ohyp;
    if any(idx_normorxic)
        pos_normorxic = pos_tumor(idx_normorxic, :);
        h_normorxic = scatter(pos_normorxic(:, 1), pos_normorxic(:, 2), 4, 'red');
    end
end


% Build legend handles conditionally
legend_handles = [];
legend_labels = {};

if ~isempty(h_vessel)
    legend_handles = [legend_handles, h_vessel];
    legend_labels = [legend_labels, {'Vessel'}];
end
if exist('h_hypoxic', 'var') && ~isempty(h_hypoxic) && isgraphics(h_hypoxic)
    legend_handles = [legend_handles, h_hypoxic];
    legend_labels = [legend_labels, {'Hypoxic'}];
end
if exist('h_normorxic', 'var') && ~isempty(h_normorxic) && isgraphics(h_normorxic)
    legend_handles = [legend_handles, h_normorxic];
    legend_labels = [legend_labels, {'Normorxic'}];
end

if ~isempty(legend_handles)
    legend(legend_handles, legend_labels, 'FontSize', 6, 'Location', 'best');
end

% title('Tumor Cells and Vessels');
hold off;

>>>>>>> d5569c1 (Initial commit)
end