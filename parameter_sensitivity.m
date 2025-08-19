clear all; close all;

filename = 'Figures/vascular/(epsilon = 1e-2)/vascular0/info.mat';
info = load(filename, 'info').info;
time_init = info.time;
params = initializeParameters();

mutation_rates = [1e-4, 1e-3, 1e-2];
prs = [0.1, 0.2, 0.3];
Sds = [1, 2, 4];
Dds = [0.25, 0.5, 1];
Sos = [2.8, 3.5, 4.2];

initsetting = struct();
initsetting.info = info;

mainFolder = 'parameter_sensitivity';
if ~exist(mainFolder, 'dir')
    mkdir(mainFolder);
end

max_time = 400;

% ----- Run tests -----
% run_and_save(mainFolder, 'mutation_rate', mutation_rates, max_time, initsetting);
% run_and_save(mainFolder, 'pr', prs, max_time, initsetting);
% run_and_save(mainFolder, 'Sd', Sds, max_time, initsetting);
run_and_save(mainFolder, 'Dd', Dds, max_time, initsetting);
% run_and_save(mainFolder, 'So', Sos, max_time, initsetting);


% ----- Helper function -----
function run_and_save(mainFolder, paramName, values, max_time, initsetting)
    for i = 1:numel(values)
        params = initializeParameters();
        params.(paramName) = values(i);
        initsetting.params = params;

        folderName = sprintf('%s/%s/%s_%.1g', mainFolder, paramName, paramName, values(i));
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end

        [snapshots, numerics, info] = main(max_time, folderName, [], initsetting);
        
        save(fullfile(folderName, 'snapshots.mat'), 'snapshots');
        save(fullfile(folderName, 'info.mat'), 'info');
        save(fullfile(folderName, 'params.mat'), 'params');
        save(fullfile(folderName, 'numerics.mat'), 'numerics');
    end
end
