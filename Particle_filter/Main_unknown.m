%% Main script for estimateing the refractory period and conduction delay 
%% of the AV node based on the ECG. Used for paper ''.
% The following toolboxes are needed to run this code:
% distrib_computing_toolbox, signal_toolbox, statistics_toolbox, parallel Computing Toolbox

clc, clear all, close all
% Obs, Estimation of f-waves and RR interval series from ECG is needed before

%% load f-waves and RR

% load data
load('Data_f_waves_RR.mat')

%% calculates psi from f-waves

% parameters for estimating psi, information is found in the function
l_short = 0.5; 
l_long = 5; 
delta_w = 1.5;
Fs = 50;

% OBS, to run this you need to have the CVX-solver installed on matlab, 
% information about this are found at https://cvxr.com/cvx/
% You can instead run this demo without this by loading the file Data_quick.mat
% load('Data_quick.mat')
[SQI, mu_2, std_1] = f_waves_to_psi(f_waves, l_short, l_long, delta_w, Fs, RR/1000);

% clear f-waves, not needed anymore
clear f_waves l_short l_long delta_w Fs

%% calculates R_last based on the fastest RR-intervals per minute

% information is found in the function
[RR_min] = RR_to_RR_min(RR);

%% Estimates theta and phi using a particle filter

% creates a folder to save the results per heartbeat
folderName = 'Results'; 

% Hyperparameters for the particle filter
N = 400; % 40000 in paper
Thr_w = .25; % 500 in paper
HP_scatter = 0.3; % 0.3 in paper
N_AA = 25; % 25 in paper
Hyp_like = 30; % 30 in paper
Type = 2; % 2 in paper
N_jump = 10; % 40 in paper

% the covariance matrix based on data from paper ''
load('Covariance.mat')
Cov_orig_new = Cov_orig*0.75;
for i = 1:12
    Cov_orig_new(i,i) = Cov_orig(i,i);
end
Cov_orig = Cov_orig_new;


% run particle filter (it saves the results in a folder)
PF_Adaptive_Unknown(RR, RR_min, N, HP_scatter, Cov_orig, Thr_w, N_jump, N_AA, Hyp_like, SQI, mu_2, std_1, Type, folderName);


%% Estimates phi using a smoothing algorithm

% Hyper parameters
M = 1000; % 10000 in paper
sigma = Cov_orig;
T = length(RR);

% run smoothing algorithm
disp('Running the smoothing algorithm')
[phi, B] = PF_smoothing_segment(sigma, M, T, folderName);

%% Removes and saves

% removes all unnecessary data to save memory
rmdir(folderName, 's')

% saves only phi and the trejectory B
save('Smoothing_results', 'phi', 'B')

%% Plot the resutls

PF_plot_smooth(phi, T)

