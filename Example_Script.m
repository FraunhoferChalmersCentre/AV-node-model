%% Example code for running the model
% Needed for the first time to create a version of the model that Matlab
% can use
mex AV_node_model.cpp -largeArrayDims

%%
clc, clear all, close all

% The minimum values, maximum prolongations, and time constants. 
% R refers to the refractory period, D to the conduction delay, and FP and
% SP to the fast and slow pathways. The values are just examples

R_FP = [300 400 250]; % Refractory period for FP
R_SP = [200 300 250]; % Refractory period for SP
D_FP = [5 7 250]; % Conduction delay for FP
D_SP = [15 7 250]; % Conduction delay for SP

Lambda = 8;

L_RR = 5000; % Length of the resulting RR interval series
R_last = 250; % Minimum RR interval from data

% 25 here makes sure that the AA series is long enough to result in a RR 
% interval series of length L_RR
AA = cumsum(-log(rand(25*L_RR, 1))/Lambda*1000+50); 

% out is the time stamps of ventricular activation from the model. 
% pathway_ind is a indicator of which pathway the impulse that activated
% the ventricle originated from
[out, diasInter, ~, pathway_ind] = AV_node_model(AA, R_FP, R_SP, R_last,...
    D_FP, D_SP, zeros(21,1), L_RR );

% Removes zeros from the created data
pathway_ind = pathway_ind(pathway_ind > 0);
diasInter = diasInter(diasInter~=0);
out = out(out>0);

% Calculating the RR series from the time stamps of ventricular activation
% in out
RR_model = diff(out);

% Plots the histogram of the simulated RR interval series
histogram(RR_model, 200:50:2000)
xlabel('RR (ms)')
