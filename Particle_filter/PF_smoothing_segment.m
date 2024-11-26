function [P_new_4_all, B_all] = PF_smoothing_segment(sigma, M, T, folderName)

% Implementation of the Forward Filtering Backward Sampling algorithm, as
% described in Algorithm 12.2 in https://link.springer.com/book/10.1007/978-3-030-47845-2

% sigma - The covariance matrix used in the particle filter
% M - Number of smoothing trajectories
% T - Number of heartbeats
% foldername - The name of the folder with the results from the particle filter


% For calculating the pdf
const = -12/2 * log(2 * pi) - 0.5 * log(det(sigma));
Sigma_inv = inv(sigma);

% Loads the T:th iteration of the particle filter
LoadName = [folderName, '\Section_', num2str(T)];
load(LoadName)

w_all = w_all; % For parafor

B_all = zeros(M, T);

% Starts the smoothing algorithm
for t = flip(1:T)

    P_4_all = P_4_all; % For parafor
    
    % samples M times from the weights in the T:th iteration of the particle filter 
    if t == T
        parfor m = 1:M
            B_all(m, T) = randsample(length(w_all), 1, true, w_all);
            P_new_4_all(t, m, :) = P_4_all(B_all(m, T), :);
        end
    else

        P_all_old = P_all;

        % Loads the data from the particle filter from the next heartbeat  
        LoadName = [folderName, '\Section_', num2str(t)];
        load(LoadName)
    
        w_all_log = log(w_all);
        x_mu = P_all';

        B_all_t = B_all(:, t+1); 
        parfor m = 1:M

            % Gets the particle from the previous iteration
            X_t1_B = P_all_old(B_all_t(m), :);
    
            % Calculates the difference between the new particles and the previous
            diff = X_t1_B' - x_mu;
            last = (Sigma_inv)*(diff);
            log_mvn = const - 0.5*dot(diff,last);
    
            % Updates the weights
            w_new = w_all_log + log_mvn';
            
            W_new = exp( w_new - max(w_new) );
            W_new = W_new/sum(W_new);
    
            % Draw new samples based on the updates weights
            B_all(m, t) = randsample(length(W_new), 1, true, W_new);

            P_new_4_all(t, m, :) = P_4_all(B_all(m, t), :);
        end

    end

end


end