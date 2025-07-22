<<<<<<< HEAD
% updatePositions.m
% Update positions of all agents according to mechanical interactions
function [tumor_cells, tip_cells, vessel_agents, negative_count] = updatePositions(tumor_cells, tip_cells, vessel_agents, chem_field, negative_count, params)

valid_tumor = ~cellfun('isempty', tumor_cells);
valid_vessel = ~cellfun('isempty', vessel_agents);

% Update tumor cell positions
if any(valid_tumor)
    % Extract positions from tumor_cells into an N x 2 matrix.
    pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
    pos_tumor = cat(1, pos_tumor{:});
    num_tumor = size(pos_tumor, 1);
    idx_tumor = find(valid_tumor);

    % Add Brownian noise and update positions
    noise = sqrt(params.dt) * params.epsilon1 * randn(num_tumor, 2);
    pos_tumor = enforceGridBounds(pos_tumor + noise, params);

    % Write back updated positions to tumor_cells
    for i = 1:length(idx_tumor)
        tumor_cells{idx_tumor(i)}.position = pos_tumor(i, :);
    end
end

% Update tip cell positions
% update the five probabilities P_{i} for i in {0, 1, 2, 3, 4}
D = params.Dn;
chi0 = params.chi0;
alpha = params.alpha;
dt = params.dt1;
dx = params.dx;
c = chem_field.TAF;

c_pad_x = [c(:, 1), c, c(:, end)];
c_pad_y = [c(1, :); c; c(end, :)];
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
    dt = params.dt1/div;
    P0 = 1-4*D*dt/dx^2 ...
        -(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1))...
        +(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2))...
        -(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:))...
        +(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));
    P0_neg_ind = P0 < 0;
    div = div + 1;
end

run_time = round(params.dt1/dt);

P1 = D*dt/dx^2-(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1));
P2 = D*dt/dx^2+(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2));
P3 = D*dt*dx^2-(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:));
P4 = D*dt/dx^2+(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));

if any(P1 < 0, 'all') || any(P2 < 0, 'all') || any(P3 < 0, 'all') || any(P4 < 0, 'all')
    negative_count = negative_count+1;
end

P1(P1 < 0) = 0;
P2(P2 < 0) = 0;
P3(P3 < 0) = 0;
P4(P4 < 0) = 0;

% Normolize the probabilities
P = P0 + P1 + P2 + P3 + P4;
Q0 = P0 ./ P;
Q1 = (P0 + P1) ./ P;
Q2 = (P0 + P1 + P2) ./ P;
Q3 = (P0 + P1 + P2 + P3) ./ P;


for ii = 1:run_time
    for i = 1:length(tip_cells)
        if ~isempty(tip_cells{i})
            pos = tip_cells{i}.position;
            pos_idx = floor(pos .* params.grid_size) + [1, 1];
            pos_ind = sub2ind(params.grid_size, pos_idx(2), pos_idx(1));
            tip_Q0 = Q0(pos_ind);
            tip_Q1 = Q1(pos_ind);
            tip_Q2 = Q2(pos_ind);
            tip_Q3 = Q3(pos_ind);
            ran_move = rand();
            if ran_move <= tip_Q0
                continue;
            elseif ran_move <= tip_Q1
                tip_cells{i}.position = enforceGridBounds(pos + [-dx, 0], params);
            elseif ran_move <= tip_Q2
                tip_cells{i}.position = enforceGridBounds(pos + [dx, 0], params);
            elseif ran_move <= tip_Q3
                tip_cells{i}.position = enforceGridBounds(pos + [0, -dx], params);
            else
                tip_cells{i}.position = enforceGridBounds(pos + [0, dx], params);
            end
            pos_new = tip_cells{i}.position;
            if (~ismember(pos_new, pos, 'rows')) && (any(valid_vessel))
                pos_vessel = cellfun(@(vc) vc.position, vessel_agents(valid_vessel), 'UniformOutput', false);
                pos_vessel = cat(1, pos_vessel{:});
                if ismember(pos_new, pos_vessel, 'rows')
                    vessel_idx = find(valid_vessel);
                    idx = vessel_idx(ismember(pos_vessel, pos_new, 'rows'));
                    order = vessel_agents{idx}.order;
                    if tip_cells{i}.order ~= order
                        tip_cells{i} = [];
                    end
                else
                    pos_idx = floor(pos_new .* params.grid_size) + [1, 1];
                    pos_ind = sub2ind(params.grid_size, pos_idx(2), pos_idx(1));
                    vessel_agents{pos_ind} = struct('position', pos_new, ...
                        'order', tip_cells{i}.order);
                end
            end
        end
    end
end
=======
% updatePositions.m
% Update positions of all agents according to mechanical interactions
function [tumor_cells, tip_cells, vessel_agents, negative_count] = updatePositions(tumor_cells, tip_cells, vessel_agents, chem_field, negative_count, params)

valid_tumor = ~cellfun('isempty', tumor_cells);
valid_vessel = ~cellfun('isempty', vessel_agents);

% Update tumor cell positions
if any(valid_tumor)
    % Extract positions from tumor_cells into an N x 2 matrix.
    pos_tumor = cellfun(@(tc) tc.position, tumor_cells(valid_tumor), 'UniformOutput', false);
    pos_tumor = cat(1, pos_tumor{:});
    num_tumor = size(pos_tumor, 1);
    idx_tumor = find(valid_tumor);

    % Add Brownian noise and update positions
    noise = sqrt(params.dt) * params.epsilon1 * randn(num_tumor, 2);
    pos_tumor = enforceGridBounds(pos_tumor + noise, params);

    % Write back updated positions to tumor_cells
    for i = 1:length(idx_tumor)
        tumor_cells{idx_tumor(i)}.position = pos_tumor(i, :);
    end
end

% Update tip cell positions
% update the five probabilities P_{i} for i in {0, 1, 2, 3, 4}
D = params.Dn;
chi0 = params.chi0;
alpha = params.alpha;
dt = params.dt1;
dx = params.dx;
c = chem_field.TAF;

c_pad_x = [c(:, 1), c, c(:, end)];
c_pad_y = [c(1, :); c; c(end, :)];
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
    dt = params.dt1/div;
    P0 = 1-4*D*dt/dx^2 ...
        -(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1))...
        +(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2))...
        -(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:))...
        +(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));
    P0_neg_ind = P0 < 0;
    div = div + 1;
end

run_time = round(params.dt1/dt);

P1 = D*dt/dx^2-(dt/(4*dx^2))*(chi_pad_x(:,3:end)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,3:end)-c_pad_x(:,2:end-1));
P2 = D*dt/dx^2+(dt/(4*dx^2))*(chi_pad_x(:,1:end-2)+chi_pad_x(:,2:end-1)).*(c_pad_x(:,2:end-1)-c_pad_x(:,1:end-2));
P3 = D*dt*dx^2-(dt/(4*dx^2))*(chi_pad_y(3:end,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(3:end,:)-c_pad_y(2:end-1,:));
P4 = D*dt/dx^2+(dt/(4*dx^2))*(chi_pad_y(1:end-2,:)+chi_pad_y(2:end-1,:)).*(c_pad_y(2:end-1,:)-c_pad_y(1:end-2,:));

if any(P1 < 0, 'all') || any(P2 < 0, 'all') || any(P3 < 0, 'all') || any(P4 < 0, 'all')
    negative_count = negative_count+1;
end

P1(P1 < 0) = 0;
P2(P2 < 0) = 0;
P3(P3 < 0) = 0;
P4(P4 < 0) = 0;

% Normolize the probabilities
P = P0 + P1 + P2 + P3 + P4;
Q0 = P0 ./ P;
Q1 = (P0 + P1) ./ P;
Q2 = (P0 + P1 + P2) ./ P;
Q3 = (P0 + P1 + P2 + P3) ./ P;


for ii = 1:run_time
    for i = 1:length(tip_cells)
        if ~isempty(tip_cells{i})
            pos = tip_cells{i}.position;
            pos_idx = floor(pos .* params.grid_size) + [1, 1];
            pos_ind = sub2ind(params.grid_size, pos_idx(2), pos_idx(1));
            tip_Q0 = Q0(pos_ind);
            tip_Q1 = Q1(pos_ind);
            tip_Q2 = Q2(pos_ind);
            tip_Q3 = Q3(pos_ind);
            ran_move = rand();
            if ran_move <= tip_Q0
                continue;
            elseif ran_move <= tip_Q1
                tip_cells{i}.position = enforceGridBounds(pos + [-dx, 0], params);
            elseif ran_move <= tip_Q2
                tip_cells{i}.position = enforceGridBounds(pos + [dx, 0], params);
            elseif ran_move <= tip_Q3
                tip_cells{i}.position = enforceGridBounds(pos + [0, -dx], params);
            else
                tip_cells{i}.position = enforceGridBounds(pos + [0, dx], params);
            end
            pos_new = tip_cells{i}.position;
            if (~ismember(pos_new, pos, 'rows')) && (any(valid_vessel))
                pos_vessel = cellfun(@(vc) vc.position, vessel_agents(valid_vessel), 'UniformOutput', false);
                pos_vessel = cat(1, pos_vessel{:});
                if ismember(pos_new, pos_vessel, 'rows')
                    vessel_idx = find(valid_vessel);
                    idx = vessel_idx(ismember(pos_vessel, pos_new, 'rows'));
                    order = vessel_agents{idx}.order;
                    if tip_cells{i}.order ~= order
                        tip_cells{i} = [];
                    end
                else
                    pos_idx = floor(pos_new .* params.grid_size) + [1, 1];
                    pos_ind = sub2ind(params.grid_size, pos_idx(2), pos_idx(1));
                    vessel_agents{pos_ind} = struct('position', pos_new, ...
                        'order', tip_cells{i}.order);
                end
            end
        end
    end
end
>>>>>>> d5569c1 (Initial commit)
end