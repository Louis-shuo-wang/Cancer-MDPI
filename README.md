# Cancer-MDPI
This is the hybrid discrete-continuous (HDC) platform constructed for the project in:
https://doi.org/10.20944/preprints202507.2175.v1,
and
https://doi.org/10.20944/preprints202508.0408.v1.

(all math notations and formulas in this readme file is typed in latex format).

The main.m is the main function in our HDC platform:
[all_snapshots, numerics, info] = main(max_time, folderName, treatment, initsetting)

This function has three arguments:
max_time is the maximum time that the main function will run (ND time = max_time \times 0.1, where 0.1 is the time step);
folderName is the file path that store the simulation results;
treatment=[a,b], where treatment(1) = treatment on period (0~50, where 5 is one treatment cycle time), and treatment(2) = S_d (drug infusion rate, and we explored S_d=2,5,10 in our project for continuous treatment plans, and Sd=2, 10/4, 10/3, 10/2 for pulsed therapies).
initsetting is a struct consisting of "info" and "params". "info" is a struct storing all information of the last simulation results based on which you want to continue your current simulation, and if you start from t=0 (no simulation history, then set info=[]; "params" is struct storing all parameter values, and you can set params=[] if you have no idea, which then will be initialized and no error occurs). (Note that Each simulation will start from the initial time t=initial_time, and spans the time interval until t=initial_time+max_time).

The first output is "all_snapshots", which is the tumor and vessel agents, TAF field, oxygen field, TAF field information at four times: (initial_time + max_time/4):max_time/4:(initial_time + max_time);
numerics is a struct consisting of all quantative metrics considered in the project, or all time series data appearing in that project.
info is the information consisting of all information of discrete and continuous parts for HDC at the lat simulation time t=initial_time+max_time. 
