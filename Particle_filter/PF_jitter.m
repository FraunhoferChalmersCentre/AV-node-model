function P_new = PF_jitter(P, N, ub, lb, HP_scatter, Cov)

% Jitters the particles of the particle filter

% P - The particles
% N - Number of particles
% ub - upper bounds
% lb - lower bounds
% HP_scatter - The hyper parameter determining how much the particles are scattered
% Cov - The covarience matrix, determining the direction of the scatter

sigma = Cov*HP_scatter;

P_new = P + randn(N,12)*chol(sigma);

x_plot = 0:50:2000;

% Make sure that the fast pathway is faster than the slow, and has a higher
% refractory period, as well as that the new particles are inside the bounds
List_ok = 1:N;
while ~isempty(List_ok) 
    
    Rf = P_new(List_ok,1) + P_new(List_ok,2).*(1-exp(-x_plot./P_new(List_ok,3)));
    Rs = P_new(List_ok,4) + P_new(List_ok,5).*(1-exp(-x_plot./P_new(List_ok,6)));
    Df = P_new(List_ok,7) + P_new(List_ok,8).*(exp(-x_plot./P_new(List_ok,9)));
    Ds = P_new(List_ok,10) + P_new(List_ok,11).*(exp(-x_plot./P_new(List_ok,12)));
    
    List_ok2 = (sum(Rs>Rf,2)>8) | (sum(Ds<Df,2)>8) | (sum(P_new(List_ok,:)<ub, 2) ~= 12) | (  sum(P_new(List_ok,:)>lb, 2) ~= 12  );
    List_ok = List_ok(List_ok2);

    P_new(List_ok,:) = P(List_ok, :) + randn(length(List_ok),12)*chol(sigma);

end

end