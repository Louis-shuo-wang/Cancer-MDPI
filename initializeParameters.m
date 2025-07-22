<<<<<<< HEAD
% initializeParameters.m
% Store all necessary parameter values in this file
function params = initializeParameters()

% Simulation grid parameters
params.grid_size = [100, 100];        % Size of simulation domain
params.dg = 1 ./ params.grid_size;
params.dy = params.dg(1);
params.dx = params.dg(2);
params.Rc = 0.5 * params.dx;         % Cell radius
params.nbins = 100;                   % number of bins to create histogram

% Time step (hours) not to small
% different time step for different agents
% not to big for molecules
params.dt = 0.1;                      % Time scale for tumor and tip agents
params.dt1 = 0.005;                   % Time scale for molecule diffusion
params.div_time = params.dt / params.dt1;

% TAF-related parametersinitial_oxygen_concentration
params.Dc = 0.12;                     % TAF diffusion coefficient    0.1
params.xic = 0.002;                   % TAF decay rate
params.eta = 500;                     % TAF production rate by hypoxic cells ---- η    5
params.lambda = 0.1;                  % TAF uptake rate by vessels -------------- λ

% Drug-related parameters
params.Dd = 0.5;                      % Drug diffusion coefficient
params.xid = 0.01;                    % Drug decay rate ---------------- ξd
params.rhod = 0.5;                    % Drug uptake rate ---------------- ρd
params.Sd = 2;                        % Drug supply rate from vessels  2
params.dexp = 0.01;                   % Drug exposure threshold
params.texp = 0.0026;                 % Time threshold for resistance (hours)
params.Deltadeath = 0;                % Death threshold increment rate ------------ Δdeath
params.pr = 0.2;                      % Drug-induced damage repair rate

% Oxygen-related parameters
params.Do = 0.64;                     % Oxygen diffusion coefficient
params.xio = 0.025;                   % Oxygen decay rate --------------- ξo
params.rhoo = 0.57;                   % Oxygen uptake rate -------------- ρo
params.So = 3.5;                      % Oxygen supply rate from vessels
params.omax = 1;                      % Maximum oxygen concentration
params.ohyp = 0.25;                   % Hypoxic threshold
params.oapop = 0.03;                  % Apoptosis threshold

% Tip cell parameters
params.Dn = 4.608e-4;                 % Tip cell random motion coefficient   0.15
params.chi0 = 0.38;                   % Base chemotaxis coefficient ------------ χ0
params.alpha = 0.6;                   % TAF sensitivity parameter
params.psi = 0.5;                     % Minimum age for branching (hours) ------ ψ
params.cbr = 1;                       % Branching rate coefficient
params.proliferation_start = 9/8;     % Time that proliferation can start 
params.doubling_time = 0.5;           % Endothelial cell doubling time

% Cell mechanical parameters
params.epsilon1 = 1e-2;               % Random motion coefficient ----------- ε1

% Tumor cell parameters
params.age = 0.5;                                                          % this equals to 8 hrs (cell cycle time)
params.resis_fraction = 0.01;                                              % fraction of pre_existing resistant cells
params.base_death_threshold = 0.5;                                         % Base death threshold
params.max_proliferation_density = 10;                                     % Maximum local cell density
params.resistance_mechanism = 'acquired';                                  % 'pre_existing' or 'acquired'
params.resistance_factor = 3;                                              % Ratio of death threshold of resistant over sensitive cells

% Mutation parameters
params.mutation_rate = 0;             % Probability of mutation during division
% params.mutation_type = 'random';      % 'linear' or 'random'               % we only consider random mutation algorithm

% Linear mutation factors (Table 1 from document)
params.mutation_oxygen_factors = [1.0, 4/3, 2.0, 4.0];
params.mutation_threshold_factors = [1.0, 2.0, 3.5, 5.0];
params.mutation_proliferation_factors = [1.0, 3/4, 3/5, 1/2];

% Random mutation constraints
params.mutation_min_factor = 0.7;
params.mutation_max_factor = 1.7;

% Vessel parameters
params.proliferation_interval = 18;   % Hours between vessel proliferation

% Initial condition parameters
params.tumor_center = [0.5, 0.5];
params.initial_tumor_number = 0;    % Initial tumor cell number   100
params.initial_tip_number = 0;        % Initial tip cell numbeer
params.initial_drug_concentration = 0; % Initial drug concentration
params.initial_oxygen_concentration = 0.5; % Initial oxygen concentration
params.initial_TAF_concentration = 1; % Initial TAF concentration
=======
% initializeParameters.m
% Store all necessary parameter values in this file
function params = initializeParameters()

% Simulation grid parameters
params.grid_size = [100, 100];        % Size of simulation domain
params.dg = 1 ./ params.grid_size;
params.dy = params.dg(1);
params.dx = params.dg(2);
params.Rc = 0.5 * params.dx;         % Cell radius
params.nbins = 100;                   % number of bins to create histogram

% Time step (hours) not to small
% different time step for different agents
% not to big for molecules
params.dt = 0.1;                      % Time scale for tumor and tip agents
params.dt1 = 0.005;                   % Time scale for molecule diffusion
params.div_time = params.dt / params.dt1;

% TAF-related parametersinitial_oxygen_concentration
params.Dc = 0.12;                     % TAF diffusion coefficient    0.1
params.xic = 0.002;                   % TAF decay rate
params.eta = 500;                     % TAF production rate by hypoxic cells ---- η    5
params.lambda = 0.1;                  % TAF uptake rate by vessels -------------- λ

% Drug-related parameters
params.Dd = 0.5;                      % Drug diffusion coefficient
params.xid = 0.01;                    % Drug decay rate ---------------- ξd
params.rhod = 0.5;                    % Drug uptake rate ---------------- ρd
params.Sd = 2;                        % Drug supply rate from vessels  2
params.dexp = 0.01;                   % Drug exposure threshold
params.texp = 0.0026;                 % Time threshold for resistance (hours)
params.Deltadeath = 0;                % Death threshold increment rate ------------ Δdeath
params.pr = 0.2;                      % Drug-induced damage repair rate

% Oxygen-related parameters
params.Do = 0.64;                     % Oxygen diffusion coefficient
params.xio = 0.025;                   % Oxygen decay rate --------------- ξo
params.rhoo = 0.57;                   % Oxygen uptake rate -------------- ρo
params.So = 3.5;                      % Oxygen supply rate from vessels
params.omax = 1;                      % Maximum oxygen concentration
params.ohyp = 0.25;                   % Hypoxic threshold
params.oapop = 0.03;                  % Apoptosis threshold

% Tip cell parameters
params.Dn = 4.608e-4;                 % Tip cell random motion coefficient   0.15
params.chi0 = 0.38;                   % Base chemotaxis coefficient ------------ χ0
params.alpha = 0.6;                   % TAF sensitivity parameter
params.psi = 0.5;                     % Minimum age for branching (hours) ------ ψ
params.cbr = 1;                       % Branching rate coefficient
params.proliferation_start = 9/8;     % Time that proliferation can start 
params.doubling_time = 0.5;           % Endothelial cell doubling time

% Cell mechanical parameters
params.epsilon1 = 1e-2;               % Random motion coefficient ----------- ε1

% Tumor cell parameters
params.age = 0.5;                                                          % this equals to 8 hrs (cell cycle time)
params.resis_fraction = 0.01;                                              % fraction of pre_existing resistant cells
params.base_death_threshold = 0.5;                                         % Base death threshold
params.max_proliferation_density = 10;                                     % Maximum local cell density
params.resistance_mechanism = 'acquired';                                  % 'pre_existing' or 'acquired'
params.resistance_factor = 3;                                              % Ratio of death threshold of resistant over sensitive cells

% Mutation parameters
params.mutation_rate = 0;             % Probability of mutation during division
% params.mutation_type = 'random';      % 'linear' or 'random'               % we only consider random mutation algorithm

% Linear mutation factors (Table 1 from document)
params.mutation_oxygen_factors = [1.0, 4/3, 2.0, 4.0];
params.mutation_threshold_factors = [1.0, 2.0, 3.5, 5.0];
params.mutation_proliferation_factors = [1.0, 3/4, 3/5, 1/2];

% Random mutation constraints
params.mutation_min_factor = 0.7;
params.mutation_max_factor = 1.7;

% Vessel parameters
params.proliferation_interval = 18;   % Hours between vessel proliferation

% Initial condition parameters
params.tumor_center = [0.5, 0.5];
params.initial_tumor_number = 0;    % Initial tumor cell number   100
params.initial_tip_number = 0;        % Initial tip cell numbeer
params.initial_drug_concentration = 0; % Initial drug concentration
params.initial_oxygen_concentration = 0.5; % Initial oxygen concentration
params.initial_TAF_concentration = 1; % Initial TAF concentration
>>>>>>> d5569c1 (Initial commit)
end