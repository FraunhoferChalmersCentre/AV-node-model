function [P_all, P_4_all, w_all, N_part_per_itt, t_sim, N_uniq, RR_est] = ...
          PF_Adaptive_known(RR_all, RR_min_all, impulses_true, N, HP_scatter, Cov_orig, Thr_w, N_jump, Hyp_like)

% A particle filter using adaptive number of particles


% Calculates the number of iterations
T = length(RR_all);

% lower and upper bounds
lb = [100 0 25 100 0 25 2 0 25 2 0 25];
ub = [1000 1000 500 1000 1000 500 50 100 500 50 100 500];

% initials the particles
P = PF_init(N, ub, lb); % OBS
P_all{T} = [];
w_all{T} = [];
P_4_all{T} = [];

% initials vectors and matrices
RR_old = zeros(N, 1);
RR_est{T} = [];

S_big = zeros(21,N);
RR_big = zeros(N,1);
Q_big{N} = []; 
FP_D = zeros(N,1); FP_R = zeros(N,1); FP_Num = zeros(N,1);
SP_D = zeros(N,1); SP_R = zeros(N,1); SP_Num = zeros(N,1);
Q_last_AA = zeros(N,1);

Sum_weights = 0;

N_uniq = zeros(T, 1);
N_part_per_itt = zeros(T, 1);
t_sim = zeros(T, 1);

 % OBS
for i = 1:N
    Q_big{i} = zeros(5,0);
end

% starts the particle filter
for t = 1:T
    
    % Uses the simulated R_last
    R_l = RR_min_all(t);

    % While-loop to run the PF until x amount of particles have been used, or until
    % the sum of weights are over a threshold
    tic
    ip_start = 1;

    if t > 1
        P = [];
        S_big = [];
        Q_last_AA = [];

        Q_big = [];
        for i = 1:N_jump
            Q_big{i} = zeros(5,0);
        end
    end

    while (ip_start < N) && (Sum_weights < Thr_w)
        % Evaluate particle

       % Creates new particles
        if t > 1

            % Create new particles and update all other inputs corresponding to the newly drawn particles
            % The current segment
            Curr_seg = ip_start:ip_start+N_jump-1;
            
            % The new particles, drawn from the index of the last iteration and jittered
            P_addon = PF_jitter(P_all{t-1}(ind_for_all(Curr_seg),:), N_jump, ub, lb, HP_scatter, Cov_orig);
            P = [P; P_addon];
              
            % Updates the states and queues in the model
            S_big = [S_big, S_big_old(:, ind(Curr_seg))];
            Q_last_AA = [Q_last_AA; Q_last_AA_old(ind(Curr_seg))];
            
            for i = ip_start:(ip_start+N_jump-1)
                Q_big{i} = Q_big_old{ind(i)};
            end

        end

        parfor i = ip_start:ip_start+N_jump-1
            
            % extracts the current particles (x), queue (Q_in), and refractory status (S_in)
            x = P(i, :);
            Q_in = Q_big{i};
            S_in = S_big(:, i);
            
            % extracted all impulses that arrives after the queue
            AA_in = impulses_true(impulses_true >= Q_last_AA(i));
            if length(AA_in) > 50
                AA_in(51:end) = [];
            end
    
            % Runs the model
            [out, o1, o2, o3, o4, o5, o6, o7, o8, o9, o10, o11, RT_out, o13, o14, o15, Q_out, FP, SP] = Model_PF(AA_in, x(1:3), x(4:6), [R_l, 0, 0],...
                x(7:9), x(10:12), [60, 0, 0], S_in, Q_in );
    
            FP_D(i) = median(FP(1,FP(1,:)>0));
            FP_R(i) = median(FP(2,FP(1,:)>0));
            FP_Num(i) = length(FP(1,FP(1,:)>0));
    
            SP_D(i) = median(SP(1,SP(1,:)>0));
            SP_R(i) = median(SP(2,SP(1,:)>0));
            SP_Num(i) = length(SP(1,SP(1,:)>0));
    
            % Saves the queue
            % Removes all AA impulses not used from the queue
            Q_ind = find(Q_out(3,:), 1, 'last');
            if isempty(Q_ind)
                Q_ind = 0;
            end
    
            Q_big{i} = Q_out(:,1:Q_ind);
    
            % Saves the arrival time for the first non-used AA interval
            Q_last_AA(i) = Q_out(1, Q_ind+1);
    
            % Saves the Refractory times used last for that particle
            S_big(:, i) = RT_out;
    
            % Saves the time of RR output (actual time, not since last RR interval) 
            RR_big(i) = out(1);
    
        end
    
        % Calculates the new starts and ends for the adaptive particle filter
        ip_start = ip_start + N_jump;      
        ip_last = ip_start - 1;
    
        % Calculates the sum of all weights
        LPr_curr = -((RR_all(t) - RR_big(1:ip_last)).^2)/(Hyp_like^2);
        Sum_weights = sum(exp(LPr_curr));
                   
    end

    % resets the sum of all weights
    Sum_weights = 0;

    % Sets the RR interval for all non-calculated particles to inf
    RR_big(ip_start:end) = inf;

    % Calulates the weights and cretes a new index based on the weights
    LPr = (-(RR_all(t) - RR_big).^2/Hyp_like^2); 
    Pr = exp(LPr-logsum(LPr));
    ind = randsample(1:N, N, true, Pr); % OBS

    [~, ind_short, ind_for_all] = unique(ind);

    % Saves the weights and particles
    P_all{t} = P(ind(ind_short),:);
    w_all{t} = Pr(ind(ind_short));

    % Saves the RP and CD for the updated particles, viewed as samples from
    % the posterior of the RP and CD in the two pathways
    FP_R = FP_R(ind(ind_short)); FP_D = FP_D(ind(ind_short)); FP_Num = FP_Num(ind(ind_short));
    SP_R = SP_R(ind(ind_short)); SP_D = SP_D(ind(ind_short)); SP_Num = SP_Num(ind(ind_short));
    P_4_all_curr = [];
    P_4_all_curr(:,1) = FP_R;   P_4_all_curr(:,2) = SP_R;   P_4_all_curr(:,3) = FP_D;   P_4_all_curr(:,4) = SP_D;
    P_4_all{t} = P_4_all_curr;

    % saves the numer of particles, the time, and the number of unique particles for one iteration
    N_part_per_itt(t) = ip_last;
    N_uniq(t) = length(ind_short);

    % Saves copies of the last iteration 
    S_big_old = S_big;
    Q_big_old = Q_big;
    Q_last_AA_old = Q_last_AA;

    % Updates and saves all RR interval information
    RR_est{t} = RR_big(ind(ind_short)) - RR_old(ind(ind_short)); 
    RR_old = RR_big(ind);

    t_sim(t) = toc;

    ['HB: ', num2str(t), ' of ', num2str(T), ' Time: ', num2str(t_sim(t)), ', # part: ', num2str(N_part_per_itt(t)), ', # unique: ', num2str(N_uniq(t))]

end

end