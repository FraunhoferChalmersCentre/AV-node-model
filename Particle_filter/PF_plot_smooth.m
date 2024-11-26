function PF_plot_smooth(P_4_all_smooth, T)

% Function for plotting the results of the particle filter and the
% smoothing

lb = [100 0 25 100 0 25 2 0 25 2 0 25];
ub = [1000 1000 500 1000 1000 500 50 100 500 50 100 500];

lb4 = [lb(1)+lb(2), lb(4)+lb(5), lb(7)+lb(8), lb(10)+lb(11)];
ub4 = [ub(1)+ub(2), ub(4)+ub(5), ub(7)+ub(8), ub(10)+ub(11)];

lb4(3:4) = lb4(3:4)*10;
ub4(3:4) = ub4(3:4)*10;

% To calculate the mode in 4 dim
N_bins = 200;
H_PF = [];

H_smooth = zeros(N_bins, T, 4);

for ik = 1:T

    for i = 1:4
        if i > 2
            [a_smooth, ~] = hist(P_4_all_smooth(ik,:,i)*10, linspace(lb4(i), ub4(i), N_bins) );
        else
            [a_smooth, ~] = hist(P_4_all_smooth(ik,:,i), linspace(lb4(i), ub4(i), N_bins) );
        end

        H_smooth(:, ik, i) = a_smooth;
    end

end

% To plot smooth results in 4 dim

Yax = {{'$R^{FP}$ (ms)'},{'$R^{SP}$ (ms)'},{'$D^{FP}$ (ms)'},{'$D^{SP}$ (ms)'}};
fig = figure;

tile = tiledlayout(2,2,'TileSpacing','Compact');

for i = 1:4

    % subplot(2,2,i)
    nexttile
    imagesc(1:T, linspace(lb4(i), ub4(i), N_bins), H_smooth(:,1:T,i)./max(H_smooth(:,1:T,i)));
    set(gca, 'YDir', 'normal');  % Flip the y-axis direction

    % xlabel('Heartbeat','Interpreter','latex')
    ylabel(Yax{i},'Interpreter','latex')

end

xlabel(tile,'Heartbeat', 'Interpreter','latex')

tile.TileSpacing = 'compact';
tile.Padding = 'compact';
set(gca,'color','none')

end
