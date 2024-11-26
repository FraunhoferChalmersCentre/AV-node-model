function [N_part_per_itt, t_sim, N_uniq, AA_out, count_bad] = ...
          PF_Adaptive_Unknown(RR_all, RR_min_all, N, HP_scatter, Cov_orig, Thr_w, N_jump, N_AA, Hyp_like, SQI, mu_2, std_1, Type, folderName)

% A particle filter using adaptive number of particles, based on paper ''

% Outputs:      
% For each heartbeat, the function saves 'P_all', 'w_all', 'P_4_all', and RR_est
% in the folder 'foldername'.

% P_all             The particles
% P_4_all           The estimated refractory period and conduction delay in the fast and slow pathway
% w_all             The weights for each particle
% RR_est            The estimated ventricular activation time

% Information about the run is saved in the following vectors:
% N_part_per_itt    The number of particles per heartbeat
% t_sim             The time of the timulation
% N_uniq            The number of unique particles per heartbeat
% AA_out            The used AA series
% count_bad         The number of times the model did not give the correct
%                    first ventricular activation time for the first heart beat

% Inputs:
% RR_all            The measured ventricular activation times
% RR_min_all        The refractory period in the coupling node
% N                 The number of particles
% HP_scatter        The hyper parameter determining how much the particles are scattered
% Cov_orig          The covarience matrix, determining the direction of the scatter
% Thr_w             The threshold for the summation of weights
% N_jump            The number of particles between each evaluation of Thr_w
% N_AA              The number of different atrial activations each particle does
% Hyp_like          The hyperparameter for the likelihood when calculating the weights
% SQI               The signal quality index (for generating AA series)
% mu_2              The mean value used to draw numbers from a normal distribution (for generating AA series)
% std_1             The standard diviation used to draw numbers from another normal distribution (for generating AA series)
% Type              Indicates if we should marginilize over P (Type == 1), or jitter for all new AA variants (Type == 2)
% folderName        The name of the folder for saving the data


if ~((Type == 1) || (Type == 2))
    error('Type needs to either be == 1 or 2. If Type == 1, we do not jitter P for each new AA variant. If Type == 2, we do.')
end


mkdir(folderName);


% Calculates the number of iterations
T = length(RR_all);

% lower and upper bounds
lb = [100 0 25 100 0 25 2 0 25 2 0 25];
ub = [1000 1000 500 1000 1000 500 50 100 500 50 100 500];
N_RR = 1;


% initials the particles
P = PF_init(N*N_AA, ub, lb);
P_old = P;
P_all = [];
w_all = [];
P_4_all = [];

% initials vectors and matrices
RR_est = [];
AA_start_list = zeros(N, 1);
ind_for_all = zeros(N,1);

S_big = zeros(21, N);
RR_big = zeros(N*N_AA, 1); % % OBS
Q_IN{N*N_AA} = [];
S_IN = zeros(21, N*N_AA);
AA_Start_time = zeros(N*N_AA, 1);
% 
% lambda_last_act = zeros(N,1);
% lambda_last_act_old = zeros(N,1);
% AA_start_time = zeros(N,1);
Q_big{N} = [];
FP_D = zeros(N*N_AA,1); FP_R = zeros(N*N_AA,1); FP_Num = zeros(N*N_AA,1);
SP_D = zeros(N*N_AA,1); SP_R = zeros(N*N_AA,1); SP_Num = zeros(N*N_AA,1);

Sum_weights = 0;
% AA_impulses_curr{N,1} = [];

N_uniq = zeros(T, 1);
N_part_per_itt = zeros(T, 1);
t_sim = zeros(T, 1);

for i = 1:N
    Q_big{i} = zeros(5,0);
end

count_bad = 0;
count_bad2 = 0;

% starts the particle filter
for t = 1:T
    % Gets the current R in the last node and lambda
    R_l = RR_min_all(t);

    % While-loop to run the PF until x amount of particles have been used, or until
    % the sum of weights are over a threshold
    tic
    ip_start = 1;

    if t > 1
        P = [];
        P_old = [];
        N_RR = 2;
    end

    while (ip_start < N) && (Sum_weights < Thr_w)
        % Evaluate particle

        if t > 1

            % Create new particles and update all other inputs corresponding to the newly drawn particles
            % The current segment
            Curr_seg = ip_start:ip_start+N_jump-1;
            

            % The new particles, drawn from the index of the last iteration and jittered
            if Type == 1
                P_addon_curr = PF_jitter(P_all(ind_for_all(Curr_seg),:), N_jump, ub, lb, HP_scatter, Cov_orig);
            end

            P_addon = zeros(N_jump*N_AA, 12);
            for i_AA = 1:N_AA
                if Type == 2
                    P_addon_curr = PF_jitter(P_all(ind_for_all(Curr_seg),:), N_jump, ub, lb, HP_scatter, Cov_orig);
                end

                fill_list = i_AA:(N_AA):(N_jump*N_AA);
                P_addon(fill_list, :) = P_addon_curr;
            end
            P = [P; P_addon];
            P_old = [P_old; P_all(ind_for_all(Curr_seg),:)];
              
        end
   
        % Gets the start and end index for the large current group
        i_weight_start = (ip_start-1)*N_AA+1;
        i_weight_end = i_weight_start+N_jump*N_AA-1;

        AA_save{i_weight_end} = [];

        parfor i_weight = i_weight_start:i_weight_end 
            
            % extracts the current queue (Q_in), and refractory status
            % (S_in) based on i, which is the index for the small group            
            i = ceil((i_weight)/N_AA);
            Q_in = Q_big{i};
            S_in = S_big(:, i);
            AA_Start_time_curr = AA_start_list(i);

            x = P(i_weight, :);
            x_old = P_old(i,:);

            % AA_in = PF_Generate_AA_from_f_waves(impulses_est, AA_Start_time_curr);
            AA_in = PF_Generate_AA_from_f_waves(SQI(t), mu_2(t), std_1(t), AA_Start_time_curr);

            % Runs the model
            [out, FP, SP, o19, o20, o21] = Model_PF_Unknown_short(AA_in, x, x_old, [R_l, 0, 0], [60, 0, 0], S_in, Q_in, N_RR);

            if t>1
                RR_prev = RR_est(ind_for_all(i));
                if out(1) ~= RR_prev
                    count_bad2 = count_bad2 + 1;
                end
            end

            FP_D(i_weight) = median(FP(1,FP(1,:)>0));
            FP_R(i_weight) = median(FP(2,FP(1,:)>0));
            FP_Num(i_weight) = length(FP(1,FP(1,:)>0));
    
            SP_D(i_weight) = median(SP(1,SP(1,:)>0));
            SP_R(i_weight) = median(SP(2,SP(1,:)>0));
            SP_Num(i_weight) = length(SP(1,SP(1,:)>0));

            % Checks if the AA impulse which activated the ventricles comes
            % from the queue, and that we got an output from the model
            if (o20(N_RR)>0) && (~isempty(o19)) && ~((length(o19) == 1) && (N_RR == 2))

                % Gets the queue (o19) at the time when the second activation of the ventricles have occured (o20(N_RR)) 
                Q_of_activ_AA = o19{o20(N_RR)};
                Q_of_activ_AA = Q_of_activ_AA( Q_of_activ_AA(:,end)==0, :);
                S_IN(:, i_weight) = o21(:, o20(N_RR));
                
                % Changes the queue to match the format for an imputed queue in the model
                Q_IN{i_weight} = Q_of_activ_AA(:, 1:5)';

                % Gets the time for the AA impulse which activated the ventricles
                AA_Start_time(i_weight) = o19{o20(N_RR)}(1);

                % Saves the AA impulses used in the model
                ind_AA = find(AA_in == AA_Start_time(i_weight));
                AA_save_curr = AA_in(1:ind_AA);
                AA_save_curr(1) = [];

                AA_save{i_weight} = AA_save_curr;
        
                % Saves the time of RR output (actual time, not since last RR interval) 
                RR_big(i_weight) = out(N_RR);

            else
                RR_big(i_weight) = inf;
                count_bad = count_bad + 1;
            end
        end
    
        % Calculates the new starts and ends for the adaptive particle filter
        ip_start = ip_start + N_jump;      
        ip_last = ip_start - 1;
    
        % Calculates the sum of all weights
        LPr_curr = -((RR_all(t) - RR_big(1:ip_last,:)).^2)/(Hyp_like^2);
        Sum_weights = sum(sum(exp(LPr_curr)));
                   
    end

    % resets the sum of all weights
    Sum_weights = 0;

    % Sets the RR interval for all non-calculated particles to inf
    RR_big(i_weight_end+1:end, :) = inf;

    % Calulates the weights and cretes a new index based on the weights
    LPr = (-(RR_all(t) - RR_big).^2/Hyp_like^2);
    Pr = exp(LPr-logsum(LPr));
    ind = randsample(1:N*N_AA, N, true, Pr);
    
    AA_start_list = AA_Start_time(ind);
    S_big = S_IN(:, ind);
    for i = 1:N
        Q_big{i} = Q_IN{ind(i)};
    end

    [~, ind_short, ind_for_all] = unique(ind);

    % Saves the weights, particles, and the AA impulses
    P_all = P(ind(ind_short),:);
    w_all = Pr(ind(ind_short));

    for i = 1:length(ind_short)
        AA_out{t}{i} = AA_save{ind(ind_short(i))};
    end

    % Saves the RP and CD for the updated particles, viewed as samples from
    % the posterior of the RP and CD in the two pathways
    FP_R = FP_R(ind(ind_short)); FP_D = FP_D(ind(ind_short)); FP_Num = FP_Num(ind(ind_short));
    SP_R = SP_R(ind(ind_short)); SP_D = SP_D(ind(ind_short)); SP_Num = SP_Num(ind(ind_short));
    P_4_all = [];
    P_4_all(:,1) = FP_R;   P_4_all(:,2) = SP_R;   P_4_all(:,3) = FP_D;   P_4_all(:,4) = SP_D;

    % saves the numer of particles, the time, and the number of unique particles for one iteration
    N_part_per_itt(t) = ip_last;
    N_uniq(t) = length(ind_short);

    % Updates and saves all RR interval information
    RR_est = RR_big(ind(ind_short)); 


    % RR_big_OLD = RR_big(ind(ind_short));
    t_sim(t) = toc;
        
    ['HB: ', num2str(t), ' of ', num2str(T), ' Time: ', num2str(t_sim(t)), ', # part: ', num2str(N_part_per_itt(t)), ', # unique: ', num2str(N_uniq(t))]

SaveName = [folderName, '\Section_', num2str(t)];
save(SaveName, 'RR_est', 'P_all', 'w_all', 'P_4_all')

end

count_bad(2) = count_bad2;

end