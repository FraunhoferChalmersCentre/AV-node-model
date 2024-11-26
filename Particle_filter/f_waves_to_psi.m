function [SQI, mu_2, std_1] = f_waves_to_psi(f_waves, l_short, l_long, delta_w, Fs, RR)

% Inputs:
% f_waves = the f-waves sampled at 1000 Hz
% RR = The RR interval series in seconds
% l_short = length of the local segment, in seconds
% l_long = length of the global segment, in seconds
% delta_w = Maximum range for the local segment, in Hz
% RR = The RR interval series in seconds

% Outputs:
% SQI = Signal quality index per heartbeat
% mu_2 = mean time between AA impulses per heartbeat
% std_1 = average standard deviation for AA impulses per heartbeat


% initilize outputs
SQI = zeros(1, length(RR));
mu_2 = zeros(1, length(RR));
std_1 = zeros(1, length(RR));

% downsample the f_waves to 50 Hz, according to [1]
s3 = decimate(f_waves, 1000/Fs);

% Mirror the end of the f_waves so it fits into 5 min segments.
% The mirrored data is not used in the particle filter, but it it used
% to be able to estimate the frequency for all segments.
% Example, if we have 59 s f_waves, we need 60 s for the last segment
% to be whole, otherwise we would only use 55 s. 
N_missing = Fs*l_long - mod(length(s3), Fs*l_long);
if N_missing == Fs*l_long
    s3_fixed = s3;
else
    s3_fixed = [s3, flip(s3(end-N_missing+1:end))];
end


% Runs the function
[Estimated_signal, ~, ~, ~, ~, Estimated_freq, ~, ~]=...
signalquality_updated_mikeal_fixed_MK(s3_fixed', l_short, l_long, delta_w, Fs);

% Unpackages the results
F_est = cell2mat(Estimated_freq);
sig_est = cell2mat(Estimated_signal');


% Loops through all heartbeats
for i_RR = 1:length(RR)
    
    if i_RR == 1
        Selected = 1:round((RR(i_RR)*Fs));
    else
        Selected = round((RR(i_RR-1)*Fs)):round((RR(i_RR)*Fs));
    end
    
    % Calculate the SQI, mu, std for each heartbeat
    SQI(i_RR) = 1 - ( std(s3_fixed(Selected) - sig_est(Selected)') / std(s3_fixed(Selected)) );
    mu_2(i_RR) = mean(1./F_est(Selected));
    std_1(i_RR) = std(1./F_est(Selected));
    
end

end

