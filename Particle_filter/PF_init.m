function P = PF_init(N, ub, lb)

% Initilizes the particles in the particle filter

% N - Number of particles
% ub - upper bounds
% lb - lower bounds

% Generates random particles
P = rand(N, 12).*(ub-lb)+lb;
x_plot = 0:50:2000;

% Make sure that the fast pathway is faster than the slow, and has a higher
% refractory period.
List_ok = 1:N;
while ~isempty(List_ok) 

    Rf = P(List_ok,1) + P(List_ok,2).*(1-exp(-x_plot./P(List_ok,3)));
    Rs = P(List_ok,4) + P(List_ok,5).*(1-exp(-x_plot./P(List_ok,6)));

    eq1 = (sum(Rs>Rf,2)>8);
    clear Rf Rs

    Df = P(List_ok,7) + P(List_ok,8).*(exp(-x_plot./P(List_ok,9)));
    Ds = P(List_ok,10) + P(List_ok,11).*(exp(-x_plot./P(List_ok,12)));

    eq2 = (sum(Ds<Df,2)>8);
    clear Df Ds

    List_ok2 = (eq1) | (eq2) | (sum(P(List_ok,:)<ub, 2) ~= 12) | (  sum(P(List_ok,:)>lb, 2) ~= 12  );
    List_ok = List_ok(List_ok2);

    P(List_ok,:) = rand(length(List_ok), 12).*(ub-lb)+lb;
end

end