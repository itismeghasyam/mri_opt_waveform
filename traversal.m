function [ delta_t, gx_wave, gy_wave] = traversal(k_curr, k_next, g_curr, slew_rates, g_max, gamma_c, t_step)

delta_kx = k_next(1) - k_curr(1);
delta_ky = k_next(2) - k_curr(2);


% To minimize total time = N*t_step, i.e to minimize N

%% For x

thresh = 10^-4;
area_diff = delta_kx;
net_area = 0;
gx_wave = [g_curr(1)];
N = 0;
while abs(area_diff) > thresh
       N = N+1;
    g_next = g_curr(1) + sign(delta_kx)*slew_rates(1)*t_step;
    
    if abs(g_next) >g_max(1)
        g_next = sign(g_next)*g_max(1);
    end
    
    curr_area = 0.5*(g_next+g_curr(1))*t_step;
    net_area = curr_area + net_area;
    
    
    if (net_area)*delta_kx > delta_kx^2
        net_area = net_area - curr_area;
        curr_area = delta_kx - net_area;
        g_next = 2*curr_area/t_step - g_curr(1);
    end
    
    gx_wave = [gx_wave; g_next];
    area_diff = abs(net_area-delta_kx);
    

end
    

%% For y
area_diff = delta_ky;
net_area = 0;
gy_wave = [g_curr(2)];
M = 0;
while abs(area_diff) > thresh
        
    M = M+1;
    g_next = g_curr(2) + sign(delta_ky)*slew_rates(2)*t_step;
    
    if abs(g_next) >g_max(2)
        g_next = sign(g_next)*g_max(2);
    end
    
    curr_area = 0.5*(g_next+g_curr(2))*t_step;
    net_area = curr_area + net_area;
    
    
    if (net_area)*delta_ky > delta_ky^2
        net_area = net_area - curr_area;
        curr_area = delta_ky - net_area;
        g_next = 2*curr_area/t_step - g_curr(2);
    end
    
    gy_wave = [gy_wave; g_next];
    area_diff = abs(net_area-delta_ky);

end

% Find the largest of M and N, to satisty the slew rates and limits of kx
% and ky

if M>N
    gx_wave = zeros(size(gy_wave));
        
    k = floor(abs(delta_kx/(t_step*g_max(1))));
    k_l = M-k;
    g_vec = linspace(g_curr(1),  sign(delta_kx)*g_max(1), k_l);
    
    gx_wave = [g_vec(:); sign(delta_kx)*g_max(1)*ones(k,1)];
    
    delta_t = M*t_step;
end
    
if N>=M
    gy_wave = zeros(size(gx_wave));
        
    k = floor(abs(delta_ky/(t_step*g_max(2))));
    k_l = N-k;
    g_vec = linspace(g_curr(2),  sign(delta_ky)*g_max(2), k_l);
    
    gy_wave = [g_vec(:); sign(delta_ky)*g_max(2)*ones(k,1)];
    
    delta_t = N*t_step;
end
    
    


end