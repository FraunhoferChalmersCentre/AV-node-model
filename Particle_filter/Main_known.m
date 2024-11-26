%% Main script for estimateing the refractory period and conduction delay 
%% of the AV node based on invasive measures. Used for paper ''.
% The following toolboxes are needed to run this code:
% distrib_computing_toolbox, signal_toolbox, statistics_toolbox, parallel Computing Toolbox

clc, clear all, close all

%% data

 load('iafdb_short_patient_1.mat')
 % Obs! You will need the cumulative RR and AA series in ms, as in impulsees_true and RR_all 

 RR = RR_all;
%% calculates R_last based on the fastest RR-intervals per minute

% information is found in the function
[RR_min] = RR_to_RR_min(RR);

%% Estimates theta and phi using a particle filter


% Hyperparameters for the particle filter
N = 100000; % The maximum number of particles (have > 100000 for good convergence)
Thr_w = 10; % The amout of added likelihood needed to start the next heartbeat (have > 500 for good convergence).
HP_scatter = 0.3; % The amout of scatter the parameters do between each heatbeat. Higher number means more scattering.
Hyp_like = 15; % The standard deviation of the weighting function. Should be large if you are uncertain about the RR interval series.
N_jump = 1000; % How many particles the algorithm runs before checking if Thr_w is meet. Good to have around 1/20 of N. 

% the covariance matrix based on data from my third paper.
load('Covariance.mat')
Cov_orig_new = Cov_orig*0.75;
for i = 1:12
    Cov_orig_new(i,i) = Cov_orig(i,i);
end
Cov_orig = Cov_orig_new;


[P_all, P_4_all, w_all, N_part_per_itt, t_sim, N_uniq, RR_est] = ...
    PF_Adaptive_known(RR, RR_min, impulses_true, N, HP_scatter, Cov_orig, Thr_w, N_jump, Hyp_like);

%% Estimates phi using a smoothing algorithm

% Hyper parameters
M = 100; % 20000 in paper
sigma = Cov_orig;
T = length(RR);

% run smoothing algorithm
disp('Running the smoothing algorithm')
[P_4_all_smooth, B_all, P_new_all] = PF_smoothing(P_all, P_4_all, w_all, sigma, M, T);


%% Plot the resutls

PF_plot_smooth(P_4_all_smooth, T)

