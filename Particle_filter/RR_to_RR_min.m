function [RR_min] = RR_to_RR_min(RR)

% Description:
% Calculated the refractory period used for the coupling node in the AV
% node model by finding the minimum RR interval for each minute.


% calculates how many 1-minute segments there are
N_one_min_seg = ceil(RR(end)/1000/60);

if N_one_min_seg == 1
    RR_min = ones(length(RR),1)*min(diff(RR));
else

% finds the minimum RR interval per minute
RR_min = zeros(1, N_one_min_seg);
for i = 1:N_one_min_seg

    % selecting the RR intervals in current minute
    ind_curr = (RR/1000 > (i-1)*60) & (RR/1000 < i*60);

    % find first RR interval
    first = find(ind_curr, 1);

    if first > 1
        ind_curr(first-1) = true;
    end

    % finding the minimum
    RR_min(i) = min(diff(RR(ind_curr)));
end

% interpolates to get a RR_min for each heartbeat
time_points = linspace(RR(1)/1000, RR(end)/1000, N_one_min_seg);
RR_min = interp1(time_points, RR_min, RR/1000, 'cubic');

end

end

