function [P_new_4_all, B_all, P_new_all] = PF_smoothing(P_all, P_4_all, w_all, sigma, M, T)

% Implementation of the Forward Filtering Backward Sampling algorithm, as
% described in Algorithm 12.2 in https://link.springer.com/book/10.1007/978-3-030-47845-2

%P_new_all = zeros(M, T, 12);
w_all_T = w_all{T};

for i = 1:T
    w_all_log{i} = log(w_all{i});
    w_all{i} = [];
end
clear w_all

% For calculating the pdf
const = -12/2 * log(2 * pi) - 0.5 * log(det(sigma));
Sigma_inv = inv(sigma);



parfor m = 1:M

    B = zeros(1, T);
    B(T) = randsample(length(w_all_T), 1, true, w_all_T);

    for t = flip(1:T-1)

        X_t1_B = P_all{t+1}(B(t+1), :);

        x_mu = P_all{t}';

        diff = X_t1_B' - x_mu;
        last = (Sigma_inv)*(diff);
        log_mvn = const - 0.5*dot(diff,last);

        w_new = w_all_log{t} + log_mvn';
        
        W_new = exp( w_new - max(w_new) );
        W_new = W_new/sum(W_new);

        B(t) = randsample(length(W_new), 1, true, W_new);
    end

    B_all(m,:) = B;

    for i = 1:T
        P_new_all(m, i, :) = P_all{i}(B(i), :);
        P_new_4_all(i, m, :) = P_4_all{i}(B(i), :);
    end
    
end


end