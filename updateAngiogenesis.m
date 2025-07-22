<<<<<<< HEAD
% processAngiogenesis.m
function [tip_cells, vessel_agents] = updateAngiogenesis(tip_cells, vessel_agents, chem_field, t, params)

% Step 1 (proe3): Update the angiogenic network based on current tip cell trajectories
% we have finished this step in updatePositions.m

% Step 2 (proe4): Handle anastomosis
% Check if any tip cell encounters an existing vessel. If so, randomly decide to remove the tip cell.
% we have finished this step in updatePositions.m

% Step 3 (proe5): Update tip cell ages
tip_cells = updateTipCellAges(tip_cells, params);

% Step 4 (dece1 & dece2, proe1a): Process branching
% For each tip cell, check if its age exceeds the minimal branching age (psi)
% and if local space is available (for example, no nearby vessels blocking the branch).
% If the branching conditions are satisfied, branch with probability proportional to TAF concentration.
[tip_cells, vessel_agents] = handleBranching(tip_cells, vessel_agents, chem_field, params);

% Step 5 (proe6): Process endothelial cell proliferation along the sprout to create new vessel agents.
vessel_agents = handleProliferation(tip_cells, vessel_agents, chem_field, t, params);
end

%% Helper Functions (remain in the same file or in their respective files)

% updateTipCellAges.m
function tip_cells = updateTipCellAges(tip_cells, params)

valid_tip = ~cellfun('isempty', tip_cells);
num_tip = sum(valid_tip);
ind_tip = find(valid_tip);

if any(valid_tip)
    age_tip = cellfun(@(tc) tc.age, tip_cells(valid_tip), 'UniformOutput', false);
    pro_time_tip = cellfun(@(tc) tc.proliferation_time, tip_cells(valid_tip), 'UniformOutput', false);
    age_tip = cat(1, age_tip{:}) + params.dt;
    pro_time_tip = cat(1, pro_time_tip{:}) + params.dt;
    for i = 1:num_tip
        tip_cells{ind_tip(i)}.age = age_tip(i);
        tip_cells{ind_tip(i)}.proliferation_time = pro_time_tip(i);
    end
end
end

% handleBranching.m
function [tip_cells, vessel_agents] = handleBranching(tip_cells, vessel_agents, chem_field, params)

for i = 1:length(tip_cells)
    if isempty(tip_cells{i})
        continue;
    end

    % Check branching: tip cell must be older than the minimum branching age (psi)
    if tip_cells{i}.age >= params.psi
        % Check local space condition
        % Determine which position to place new tip cells
        position = full(tip_cells{i}.position);
        direction = [-1, -1, -1, 1; ...
            1, -1, -1, 1; ...
            -1, 1, 1, 1;...
            -1, -1, 1, -1];

        % determine the tip cell branch direction along the gradient of TAF
        D = params.Dn;
        chi0 = params.chi0;
        alpha = params.alpha;
        dt = params.dt1;
        dx = params.dx;
        c = chem_field.TAF;
        pos = tip_cells{i}.position;
        pos_idx = floor(pos .* params.grid_size) + [1, 1];
        pos_ind = sub2ind(params.grid_size, pos_idx(2), pos_idx(1));

        c_pad_x = [c(:,1),c,c(:,end)];
        c_pad_y = [c(1,:);c;c(end,:)];
        chi = chi0./(1+alpha*c);
        chi_pad_x = [chi(:,1),chi,chi(:,end)];
        chi_pad_y = [chi(1,:);chi;chi(end,:)];

        % The second order derivative is updated by central difference operator,
        % which may be unstable in our implicit discretization
        P0 = 1-4*D*dt/dx^2 ...
            -(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1))...
            +(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2))...
            -(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:))...
            +(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));
        P0_neg_ind = P0 < 0;
        div = 2;

        while any(P0_neg_ind, 'all')
            dt = params.dt1 / div;
            P0 = 1-4*D*dt/dx^2 ...
                -(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1))...
                +(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2))...
                -(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:))...
                +(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));
            P0_neg_ind = P0 < 0;
            div = div + 1;
        end

        P1 = D*dt/dx^2-(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1));
        P2 = D*dt/dx^2+(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2));
        P3 = D*dt*dx^2-(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:));
        P4 = D*dt/dx^2+(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));

        P1(P1 < 0) = 0;
        P2(P2 < 0) = 0;
        P3(P3 < 0) = 0;
        P4(P4 < 0) = 0;

        P1 = P1(pos_ind);
        P2 = P2(pos_ind);
        P3 = P3(pos_ind);
        P4 = P4(pos_ind);
        P = [P1, P2, P3, P4];

        pos1_choose = [];
        pos2_choose = [];

        for j = 1:size(direction, 1)
            [~, ind] = max(P);
            P(ind) = [];

            valid_vessel = ~cellfun('isempty', vessel_agents);
            if ~(any(valid_vessel))
                vessel_network = [];
            else
                vessel_network = cellfun(@(vc) vc.position, vessel_agents(valid_vessel), 'UniformOutput', false);
                vessel_network = cat(1, vessel_network{:});
            end

            pos1 = position + direction(ind, 1:2) .* params.dg;
            pos2 = position + direction(ind, 3:4) .* params.dg;
            pos1_floor = floor(enforceGridBounds(pos1, params) .* params.grid_size) .* params.dg + [0.5, 0.5] .* params.dg;
            pos2_floor = floor(enforceGridBounds(pos2, params) .* params.grid_size) .* params.dg + [0.5, 0.5] .* params.dg;

            if (~(ismember(pos1_floor, vessel_network, 'rows'))) || (~(ismember(pos2_floor, vessel_network, 'rows')))
                if ~(ismember(pos1_floor, vessel_network, 'rows'))
                    pos1_choose = pos1_floor;
                end
                if ~(ismember(pos2_floor, vessel_network, 'rows'))
                    pos2_choose = pos2_floor;
                end
                break;
            end
            direction(ind, :) = [];
        end

        % Calculate local TAF concentration from the chemical field
        % TAF_conc = interp2(params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2, ...
        %     (params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2)', chem_field.TAF, ...
        %     position(1), position(2), 'makima');
        % Determine branching probability (scaled by TAF concentration)
        pos_vessel = floor(position .* params.grid_size) + [1, 1];
        TAF_con = chem_field.TAF(pos_vessel(2), pos_vessel(1));                % TAF concentration at current position
        TAF_norm = (max(chem_field.TAF, [], 'all') + 1e-6);                    % avoid division by 0
        con_normal = TAF_con / TAF_norm;                                       % normalized TAF concentration
        branch_prob = params.cbr * con_normal;

        % Attempt branching using a Poisson process approximation
        if rand() <= branch_prob * params.dt && (~isempty(pos1_choose) || ~isempty(pos2_choose))

            % Create a new tip cell (branch) at the same location with reset age
            new_tip.age = 0;
            new_tip.proliferation_time = tip_cells{i}.proliferation_time;
            % Add new tip cells to the existing tip cells array
            empty_ind = find(cellfun('isempty', tip_cells), 1);

            tip_cells{i}.age = 0;
            tip_cells{i}.position = pos1_choose;

            if ~isempty(pos1_choose)
                pos1_idx = floor(enforceGridBounds(pos1_choose, params) .* params.grid_size) + [1, 1];
                ind1 = sub2ind(params.grid_size, pos1_idx(2), pos1_idx(1));
                tip_cells{i}.position = pos1_choose;
                % tip_cells{i}.order = order_max + 1;
                vessel_agents{ind1} = struct('position', pos1_choose, ...
                    'order', tip_cells{i}.order);
                if ~isempty(pos2_choose)
                    pos2_idx = floor(enforceGridBounds(pos2_choose, params) .* params.grid_size) + [1, 1];
                    ind2 = sub2ind(params.grid_size, pos2_idx(2), pos2_idx(1));
                    new_tip.position = pos2_choose;
                    new_tip.order = tip_cells{i}.order;
                    if ~isempty(empty_ind)
                        tip_cells{empty_ind} = new_tip;
                    else
                        tip_cells{end+1} = new_tip;
                    end
                    vessel_agents{ind2} = struct('position', pos2_choose, ...
                        'order', tip_cells{i}.order);
                end
            else
                if ~isempty(pos2_choose)
                    pos2_idx = floor(enforceGridBounds(pos2_choose, params) .* params.grid_size) + [1, 1];
                    ind2 = sub2ind(params.grid_size, pos2_idx(2), pos2_idx(1));
                    tip_cells{i}.position = pos2_choose;
                    vessel_agents{ind2} = struct('position', pos2_choose, ...
                        'order', tip_cells{i}.order);
                end
            end
        end
    end
end
end

% handleProliferation.m
function vessel_agents = handleProliferation(tip_cells, vessel_agents, chem_field, t, params)

% At specified intervals, proliferate existing vessel agents along the sprout pathways
for i = 1:length(tip_cells)
    if isempty(tip_cells{i})
        continue;
    end
    % Determine proliferation direction based on current position, order
    % the four orthogonal directions by left, right, up, and down.
    direction = [[-1, 0]; ...
        [1, 0]; ...
        [0, -1]; ...
        [0, 1]] .* params.dg;
    pos_choose = determineValid(tip_cells{i}, vessel_agents, chem_field, direction, params);
    if (~isempty(pos_choose)) && (t >= params.proliferation_start) && (tip_cells{i}.proliferation_time >= params.doubling_time)
        tip_cells{i}.position = pos_choose;
        tip_cells{i}.proliferation_time = 0;
        pos_idx = floor(pos_choose .* params.grid_size + [1, 1]);
        ind = sub2ind(params.grid_size, pos_idx(2), pos_idx(1));
        vessel_agents{ind} = struct('position', pos_choose, ...
            'order', tip_cells{i}.order);
    end
end
end

% Helper functions
function pos_choose = determineValid(cell, vessel_agents, chem_field, direction, params)

D = params.Dn;
chi0 = params.chi0;
alpha = params.alpha;
% rho0 = params.rho0;
dt = params.dt1;
% dx = params.dx;
dx = params.dx;
c = chem_field.TAF;
pos = cell.position;
pos_idx = floor(pos .* params.grid_size) + [1, 1];
pos_ind = sub2ind(params.grid_size, pos_idx(2), pos_idx(1));

c_pad_x = [c(:, 2), c, c(:, end - 1)];
c_pad_y = [c(2, :); c; c(end - 1, :)];

chi = chi0 ./ (1 + alpha * c);
chi_pad_x = [chi(:,1),chi,chi(:,end)];
chi_pad_y = [chi(1,:);chi;chi(end,:)];

% The second order derivative is updated by central difference operator,
% which may be unstable in our implicit discretization
P0 = 1-4*D*dt/dx^2 ...
    -(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1))...
    +(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2))...
    -(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:))...
    +(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));
P0_neg_ind = P0 < 0;
div = 2;

while any(P0_neg_ind, 'all')
    dt = params.dt1/div;
    P0 = 1-4*D*dt/dx^2 ...
        -(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1))...
        +(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2))...
        -(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:))...
        +(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));
    P0_neg_ind = P0 < 0;
    div = div + 1;
end

P1 = D*dt/dx^2-(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1));
P2 = D*dt/dx^2+(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2));
P3 = D*dt*dx^2-(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:));
P4 = D*dt/dx^2+(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));

P1(P1 < 0) = 0;
P2(P2 < 0) = 0;
P3(P3 < 0) = 0;
P4(P4 < 0) = 0;

P1 = P1(pos_ind);
P2 = P2(pos_ind);
P3 = P3(pos_ind);
P4 = P4(pos_ind);
P = [P1, P2, P3, P4];
pos_choose = [];

for i = 1:4
    [~, ind] = max(P);
    P(ind) = [];

    vessel_network = cellfun(@(vc) vc.position, vessel_agents(~cellfun('isempty', vessel_agents)), 'UniformOutput', false);
    vessel_network = cat(1, vessel_network{:});
    pos_floor = floor(enforceGridBounds(pos + direction(ind, :), params) .* params.grid_size) .* params.dg + [0.5, 0.5] .* params.dg;

    if ~ismember(pos_floor, vessel_network, 'rows')
        pos_choose = pos_floor;
        break;
    end
end
=======
% processAngiogenesis.m
function [tip_cells, vessel_agents] = updateAngiogenesis(tip_cells, vessel_agents, chem_field, t, params)

% Step 1 (proe3): Update the angiogenic network based on current tip cell trajectories
% we have finished this step in updatePositions.m

% Step 2 (proe4): Handle anastomosis
% Check if any tip cell encounters an existing vessel. If so, randomly decide to remove the tip cell.
% we have finished this step in updatePositions.m

% Step 3 (proe5): Update tip cell ages
tip_cells = updateTipCellAges(tip_cells, params);

% Step 4 (dece1 & dece2, proe1a): Process branching
% For each tip cell, check if its age exceeds the minimal branching age (psi)
% and if local space is available (for example, no nearby vessels blocking the branch).
% If the branching conditions are satisfied, branch with probability proportional to TAF concentration.
[tip_cells, vessel_agents] = handleBranching(tip_cells, vessel_agents, chem_field, params);

% Step 5 (proe6): Process endothelial cell proliferation along the sprout to create new vessel agents.
vessel_agents = handleProliferation(tip_cells, vessel_agents, chem_field, t, params);
end

%% Helper Functions (remain in the same file or in their respective files)

% updateTipCellAges.m
function tip_cells = updateTipCellAges(tip_cells, params)

valid_tip = ~cellfun('isempty', tip_cells);
num_tip = sum(valid_tip);
ind_tip = find(valid_tip);

if any(valid_tip)
    age_tip = cellfun(@(tc) tc.age, tip_cells(valid_tip), 'UniformOutput', false);
    pro_time_tip = cellfun(@(tc) tc.proliferation_time, tip_cells(valid_tip), 'UniformOutput', false);
    age_tip = cat(1, age_tip{:}) + params.dt;
    pro_time_tip = cat(1, pro_time_tip{:}) + params.dt;
    for i = 1:num_tip
        tip_cells{ind_tip(i)}.age = age_tip(i);
        tip_cells{ind_tip(i)}.proliferation_time = pro_time_tip(i);
    end
end
end

% handleBranching.m
function [tip_cells, vessel_agents] = handleBranching(tip_cells, vessel_agents, chem_field, params)

for i = 1:length(tip_cells)
    if isempty(tip_cells{i})
        continue;
    end

    % Check branching: tip cell must be older than the minimum branching age (psi)
    if tip_cells{i}.age >= params.psi
        % Check local space condition
        % Determine which position to place new tip cells
        position = full(tip_cells{i}.position);
        direction = [-1, -1, -1, 1; ...
            1, -1, -1, 1; ...
            -1, 1, 1, 1;...
            -1, -1, 1, -1];

        % determine the tip cell branch direction along the gradient of TAF
        D = params.Dn;
        chi0 = params.chi0;
        alpha = params.alpha;
        dt = params.dt1;
        dx = params.dx;
        c = chem_field.TAF;
        pos = tip_cells{i}.position;
        pos_idx = floor(pos .* params.grid_size) + [1, 1];
        pos_ind = sub2ind(params.grid_size, pos_idx(2), pos_idx(1));

        c_pad_x = [c(:,1),c,c(:,end)];
        c_pad_y = [c(1,:);c;c(end,:)];
        chi = chi0./(1+alpha*c);
        chi_pad_x = [chi(:,1),chi,chi(:,end)];
        chi_pad_y = [chi(1,:);chi;chi(end,:)];

        % The second order derivative is updated by central difference operator,
        % which may be unstable in our implicit discretization
        P0 = 1-4*D*dt/dx^2 ...
            -(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1))...
            +(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2))...
            -(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:))...
            +(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));
        P0_neg_ind = P0 < 0;
        div = 2;

        while any(P0_neg_ind, 'all')
            dt = params.dt1 / div;
            P0 = 1-4*D*dt/dx^2 ...
                -(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1))...
                +(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2))...
                -(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:))...
                +(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));
            P0_neg_ind = P0 < 0;
            div = div + 1;
        end

        P1 = D*dt/dx^2-(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1));
        P2 = D*dt/dx^2+(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2));
        P3 = D*dt*dx^2-(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:));
        P4 = D*dt/dx^2+(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));

        P1(P1 < 0) = 0;
        P2(P2 < 0) = 0;
        P3(P3 < 0) = 0;
        P4(P4 < 0) = 0;

        P1 = P1(pos_ind);
        P2 = P2(pos_ind);
        P3 = P3(pos_ind);
        P4 = P4(pos_ind);
        P = [P1, P2, P3, P4];

        pos1_choose = [];
        pos2_choose = [];

        for j = 1:size(direction, 1)
            [~, ind] = max(P);
            P(ind) = [];

            valid_vessel = ~cellfun('isempty', vessel_agents);
            if ~(any(valid_vessel))
                vessel_network = [];
            else
                vessel_network = cellfun(@(vc) vc.position, vessel_agents(valid_vessel), 'UniformOutput', false);
                vessel_network = cat(1, vessel_network{:});
            end

            pos1 = position + direction(ind, 1:2) .* params.dg;
            pos2 = position + direction(ind, 3:4) .* params.dg;
            pos1_floor = floor(enforceGridBounds(pos1, params) .* params.grid_size) .* params.dg + [0.5, 0.5] .* params.dg;
            pos2_floor = floor(enforceGridBounds(pos2, params) .* params.grid_size) .* params.dg + [0.5, 0.5] .* params.dg;

            if (~(ismember(pos1_floor, vessel_network, 'rows'))) || (~(ismember(pos2_floor, vessel_network, 'rows')))
                if ~(ismember(pos1_floor, vessel_network, 'rows'))
                    pos1_choose = pos1_floor;
                end
                if ~(ismember(pos2_floor, vessel_network, 'rows'))
                    pos2_choose = pos2_floor;
                end
                break;
            end
            direction(ind, :) = [];
        end

        % Calculate local TAF concentration from the chemical field
        % TAF_conc = interp2(params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2, ...
        %     (params.dx / 2:params.dx:params.dx * params.grid_size(1) - params.dx / 2)', chem_field.TAF, ...
        %     position(1), position(2), 'makima');
        % Determine branching probability (scaled by TAF concentration)
        pos_vessel = floor(position .* params.grid_size) + [1, 1];
        TAF_con = chem_field.TAF(pos_vessel(2), pos_vessel(1));                % TAF concentration at current position
        TAF_norm = (max(chem_field.TAF, [], 'all') + 1e-6);                    % avoid division by 0
        con_normal = TAF_con / TAF_norm;                                       % normalized TAF concentration
        branch_prob = params.cbr * con_normal;

        % Attempt branching using a Poisson process approximation
        if rand() <= branch_prob * params.dt && (~isempty(pos1_choose) || ~isempty(pos2_choose))

            % Create a new tip cell (branch) at the same location with reset age
            new_tip.age = 0;
            new_tip.proliferation_time = tip_cells{i}.proliferation_time;
            % Add new tip cells to the existing tip cells array
            empty_ind = find(cellfun('isempty', tip_cells), 1);

            tip_cells{i}.age = 0;
            tip_cells{i}.position = pos1_choose;

            if ~isempty(pos1_choose)
                pos1_idx = floor(enforceGridBounds(pos1_choose, params) .* params.grid_size) + [1, 1];
                ind1 = sub2ind(params.grid_size, pos1_idx(2), pos1_idx(1));
                tip_cells{i}.position = pos1_choose;
                % tip_cells{i}.order = order_max + 1;
                vessel_agents{ind1} = struct('position', pos1_choose, ...
                    'order', tip_cells{i}.order);
                if ~isempty(pos2_choose)
                    pos2_idx = floor(enforceGridBounds(pos2_choose, params) .* params.grid_size) + [1, 1];
                    ind2 = sub2ind(params.grid_size, pos2_idx(2), pos2_idx(1));
                    new_tip.position = pos2_choose;
                    new_tip.order = tip_cells{i}.order;
                    if ~isempty(empty_ind)
                        tip_cells{empty_ind} = new_tip;
                    else
                        tip_cells{end+1} = new_tip;
                    end
                    vessel_agents{ind2} = struct('position', pos2_choose, ...
                        'order', tip_cells{i}.order);
                end
            else
                if ~isempty(pos2_choose)
                    pos2_idx = floor(enforceGridBounds(pos2_choose, params) .* params.grid_size) + [1, 1];
                    ind2 = sub2ind(params.grid_size, pos2_idx(2), pos2_idx(1));
                    tip_cells{i}.position = pos2_choose;
                    vessel_agents{ind2} = struct('position', pos2_choose, ...
                        'order', tip_cells{i}.order);
                end
            end
        end
    end
end
end

% handleProliferation.m
function vessel_agents = handleProliferation(tip_cells, vessel_agents, chem_field, t, params)

% At specified intervals, proliferate existing vessel agents along the sprout pathways
for i = 1:length(tip_cells)
    if isempty(tip_cells{i})
        continue;
    end
    % Determine proliferation direction based on current position, order
    % the four orthogonal directions by left, right, up, and down.
    direction = [[-1, 0]; ...
        [1, 0]; ...
        [0, -1]; ...
        [0, 1]] .* params.dg;
    pos_choose = determineValid(tip_cells{i}, vessel_agents, chem_field, direction, params);
    if (~isempty(pos_choose)) && (t >= params.proliferation_start) && (tip_cells{i}.proliferation_time >= params.doubling_time)
        tip_cells{i}.position = pos_choose;
        tip_cells{i}.proliferation_time = 0;
        pos_idx = floor(pos_choose .* params.grid_size + [1, 1]);
        ind = sub2ind(params.grid_size, pos_idx(2), pos_idx(1));
        vessel_agents{ind} = struct('position', pos_choose, ...
            'order', tip_cells{i}.order);
    end
end
end

% Helper functions
function pos_choose = determineValid(cell, vessel_agents, chem_field, direction, params)

D = params.Dn;
chi0 = params.chi0;
alpha = params.alpha;
% rho0 = params.rho0;
dt = params.dt1;
% dx = params.dx;
dx = params.dx;
c = chem_field.TAF;
pos = cell.position;
pos_idx = floor(pos .* params.grid_size) + [1, 1];
pos_ind = sub2ind(params.grid_size, pos_idx(2), pos_idx(1));

c_pad_x = [c(:, 2), c, c(:, end - 1)];
c_pad_y = [c(2, :); c; c(end - 1, :)];

chi = chi0 ./ (1 + alpha * c);
chi_pad_x = [chi(:,1),chi,chi(:,end)];
chi_pad_y = [chi(1,:);chi;chi(end,:)];

% The second order derivative is updated by central difference operator,
% which may be unstable in our implicit discretization
P0 = 1-4*D*dt/dx^2 ...
    -(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1))...
    +(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2))...
    -(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:))...
    +(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));
P0_neg_ind = P0 < 0;
div = 2;

while any(P0_neg_ind, 'all')
    dt = params.dt1/div;
    P0 = 1-4*D*dt/dx^2 ...
        -(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1))...
        +(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2))...
        -(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:))...
        +(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));
    P0_neg_ind = P0 < 0;
    div = div + 1;
end

P1 = D*dt/dx^2-(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1));
P2 = D*dt/dx^2+(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2));
P3 = D*dt*dx^2-(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:));
P4 = D*dt/dx^2+(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));

P1(P1 < 0) = 0;
P2(P2 < 0) = 0;
P3(P3 < 0) = 0;
P4(P4 < 0) = 0;

P1 = P1(pos_ind);
P2 = P2(pos_ind);
P3 = P3(pos_ind);
P4 = P4(pos_ind);
P = [P1, P2, P3, P4];
pos_choose = [];

for i = 1:4
    [~, ind] = max(P);
    P(ind) = [];

    vessel_network = cellfun(@(vc) vc.position, vessel_agents(~cellfun('isempty', vessel_agents)), 'UniformOutput', false);
    vessel_network = cat(1, vessel_network{:});
    pos_floor = floor(enforceGridBounds(pos + direction(ind, :), params) .* params.grid_size) .* params.dg + [0.5, 0.5] .* params.dg;

    if ~ismember(pos_floor, vessel_network, 'rows')
        pos_choose = pos_floor;
        break;
    end
end
>>>>>>> d5569c1 (Initial commit)
end