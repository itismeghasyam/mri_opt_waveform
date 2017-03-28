function [ delta_t, g_next, psblty ] = traversal(k_curr, k_next, g_curr, slew_rates, g_max, gamma_c);

if k_curr == k_next
    disp('The current and next points in k-space are the same');
    delta_t = 0;
    g_next = [0 0];
    return
end


delta_kx = k_next(1)-k_curr(1);
delta_ky = k_next(2)-k_curr(2);

%% For k-x

% Min time with slew rate limitation
if delta_kx >0
    delta_gx = slew_rates(1);
else
    delta_gx = -slew_rates(1); % assuming the change is within the slew rates;
end

if delta_kx > 0 && g_curr(1) == g_max(1) % moving at max speed
    delta_gx = 0;
end
if delta_kx < 0 && g_curr(1) == -g_max(1) % moving at max speed in opp dir
    delta_gx = 0;
end

if delta_gx == 0
    delta_tx_s = max(0,delta_kx/(gamma_c*g_curr(1)));
else
    delta_tx_s = (-g_curr(1) + sqrt((g_curr(1))^2 + 2*gamma_c*delta_kx*delta_gx))/delta_gx;
    delta_tx_s = min(delta_tx_s, (-g_curr(1) - sqrt( (g_curr(1))^2 + 2*gamma_c*delta_kx*delta_gx))/delta_gx);
    
    if delta_tx_s <0
        delta_tx_s = max(0,(-g_curr(1) + sqrt( (g_curr(1))^2 + 2*gamma_c*delta_kx*delta_gx))/delta_gx);
    end
end


% Min time with max gradient limitation
delta_tx_g = 2*delta_kx/(gamma_c*(g_curr(1) + sign(delta_kx)*g_max(1)));
delta_tx_g = min( delta_tx_g ,2*delta_kx/(gamma_c*(g_curr(1) - sign(delta_kx)*g_max(1))));

 if delta_tx_g <0
        delta_tx_g = max(0,2*delta_kx/(gamma_c*(g_curr(1) + sign(delta_kx)*g_max(1))));
 end

% The time that satisfies both the max gradient and max slew rate
% limitation
delta_tx = max(delta_tx_g, delta_tx_s);

if abs(delta_tx) == Inf
    delta_tx = min(delta_tx_g, delta_tx_s);
end


if round(abs(delta_kx)) == 0
    delta_tx = 0;
end


%% For k-y

% Min time with slew rate limitation
if delta_ky >0
    delta_gy = slew_rates(2);
else
    delta_gy = -slew_rates(2); % assuming the change is within the slew rates;
end

if delta_ky > 0 && g_curr(2) == g_max(2) % moving at max speed
    delta_gy = 0;
end
if delta_ky < 0 && g_curr(2) == -g_max(2) % moving at max speed in opp dir
    delta_gy = 0;
end

if delta_gy == 0
    delta_ty_s = max(0,delta_ky/(gamma_c*g_curr(2)));
else
    delta_ty_s = (-g_curr(2) + sqrt( (g_curr(2))^2 + 2*gamma_c*delta_ky*delta_gy))/delta_gy;
    delta_ty_s = min(delta_ty_s, (-g_curr(2) - sqrt( (g_curr(2))^2 + 2*gamma_c*delta_ky*delta_gy))/delta_gy);
    
    if delta_ty_s <0
        delta_ty_s = max(0,(-g_curr(2) + sqrt( (g_curr(2))^2 + 2*gamma_c*delta_ky*delta_gy))/delta_gy);
    end
    
end


% Min time with max gradient limitation
delta_ty_g = 2*delta_ky/(gamma_c*(g_curr(2) + sign(delta_ky)*g_max(2)));
delta_ty_g = min( delta_ty_g ,2*delta_ky/(gamma_c*(g_curr(2) - sign(delta_ky)*g_max(2))));

 if delta_ty_g <0
        delta_ty_g = max(0,2*delta_ky/(gamma_c*(g_curr(2) + sign(delta_ky)*g_max(2))));
 end

% The time that satisfies both the max gradient and max slew rate
% limitation
delta_ty = max(delta_ty_g, delta_ty_s);

if abs(delta_ty) == Inf
    delta_ty = min(delta_ty_g, delta_ty_s);
end

if round(abs(delta_ky)) == 0
    delta_ty = 0;
end



%% Finding common travel time

% Now delta_ty and delta_tx give the minimum times it takes for the
% traversal in k-space, if it were to be different points, but both the
% times are for the same point in k-space, so we look at the maximum of
% both the times and use that to recalculate the delta_gx and delta_gy

if delta_ty >= delta_tx 
    delta_t = delta_ty;
    delta_gx = (2/delta_t^2)*(delta_kx/gamma_c - g_curr(1)*delta_t);
    delta_gy = (2/delta_t^2)*(delta_ky/gamma_c - g_curr(2)*delta_t);
else
    delta_t = delta_tx;
    delta_gx = (2/delta_t^2)*(delta_kx/gamma_c - g_curr(1)*delta_t);
    delta_gy = (2/delta_t^2)*(delta_ky/gamma_c - g_curr(2)*delta_t);
end


%% Possibility check

% Check if delta_gx and delta_gy confirm with norms of slew rate

if ((abs(delta_gx) > slew_rates(1)) || (abs(delta_gy) > slew_rates(2)))
    psblty = 1;
else
    psblty = 0;
end

% Other conditions can be added which check the current position with
% respect to the farthest point in the grid and then loop back to give
% gradient waveforms and optimal time wrt the current position in the grid,
% not just the two points.



%% Outputs
% delta_t = max(delta_tx, delta_ty);
delta_g = [delta_gx, delta_gy];
g_next = g_curr(:) + delta_t*delta_g(:);

end