function AA_in = PF_Generate_AA_from_f_waves(SQI, mu_2, std_1, AA_Start_time_curr)

% A function to generate an AA series based on the previous
% AA intervall that activated the previous heartbeat, and the freqency of
% the f_waves

% SQI - Signal quality index
% mu_2 - the mean value used to draw numbers from a normal distribution
% std_1 - the standard diviation used to draw numbers from another normal distribution
% AA_Start_time_curr - The known starting time of the generated AA series



% changes the SQI so that low values means low noise, and makes ever signal
% with a quality index over 0.3 have zero 'noise', aka we trust the signal
% if is is over 0.3.
SQI_changed = (1 - SQI) - 0.7;
SQI_changed = max(SQI_changed, 0);

% scale the reversed SQI to make it useable as a std for our normal distribution
SQI_changed = SQI_changed^2;


% draws a mu_1 based on mu_2
mu_1 = mu_2 + randn(1,1)*SQI_changed;

% mu can not be lower than 4 Hz or higher than 10 Hz.
mu_1 = min(mu_1, 0.25);
mu_1 = max(mu_1, 0.1);

std_1 = std_1*4.5;

r_norm = [];
N_ok = 0;
N_loop = 0;
while N_ok < 30
    % generates 100 impulses with the calculated normal distribution
    r_norm = [r_norm; randn(100, 1)*std_1 + mu_1];
    
    % removes all impulses under 50 ms
    r_norm(r_norm < 0.05) = [];

    N_ok = length(r_norm);

    N_loop = N_loop + 1;

    % if we cant generate an AA series, we create an empty vector
    if N_loop == 100
        break
    end
end

if N_loop == 100
    r_norm = [];
else
    r_norm = r_norm(1:30);
end

AA_in = [AA_Start_time_curr; AA_Start_time_curr+cumsum(r_norm)*1000];

end
