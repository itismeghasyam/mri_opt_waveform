%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
This code was possible because of the optimization routine written by Miki
Lustig (Stanford MRI) to find optimal gradient waveforms based on the slew rate and max
gradient constraints
%}

% Meghasyam Tummalacherla

% Parameter limits

gx_max = 15*10^-3; % 15 mT/m
gy_max = 15*10^-3; % 15 mT/m
slew_x_max = 10; % 100 T/m/s
slew_y_max = 10; % 100 T/m/s


g_max = [gx_max, gy_max];
slew_rates = [slew_x_max, slew_y_max];

gamma_c = 42.57*10^6; % 42.57*10^6 /s/T

% Full grid pattern
grid_size = 4;
full_grid = ones(grid_size, grid_size);
figure(1);
subplot(1,2,1);
imagesc(full_grid);


% The postions of the points on the grid in k-space
fovx = 10*10^-2; % 10cm
fovy = 10*10^-2; % 10cm
delta_kx = 1/fovx;
delta_ky = 1/fovy;

% Only the relative positions of the points in the k-space matters to find
% the optimal path as we can center the k-space later
points_pos = cell(1, grid_size^2);
for i =1:grid_size
    for j=1:grid_size
        temp_pos = [delta_kx*i, delta_ky*j];
        points_pos{(i-1)*grid_size + j} = temp_pos;
    end
end


% Undersampling grid
% N = round(grid_size^2 /2); % %50 percent undersampling
% U = zeros(1, grid_size^2);
% U(randi(grid_size^2,[1 N])) = 1;
% U = reshape(U, [grid_size, grid_size]);
% U = full_grid;
% N = sum(sum(U)); % Number of points on the grid

U = [1 0 0 1; 1 1 1 0 ;0 1 1 1;1 0 1 0];
N = sum(sum(U));
figure(1);
subplot(1,2,2);
imagesc(U);
colormap('gray');

points_pos_U = cell(1, N);
count = 1;
for i =1:grid_size
    for j=1:grid_size
        if U(i,j) > 0
            temp_pos = [delta_kx*i, delta_ky*j].*U(i,j);
            points_pos_U{count} = temp_pos;
            count = count+1;
        end
    end
end

% P_mat = perms([1:N]); % The matrix that contains all possible paths amongst the N points in the undersampled grid

% time limits
fs_adc = 1*10; % The Machine can move at a max of 10hz
t_sampling = 1/fs_adc;
% t_max = traversal(full_grid);
t_step = t_sampling; % Usually the sampling time is much much faster than the time it takes the gradients to move

% We are not using these values for time limits as our assumption of the
% algorithm is such that the ADC is much faster than the gradient waveform
% change

% assuming the gradient changes from one sample in k-space to the other in
% linear steps only. I.e the Gx curve is piece wise linear and so is the Gy
% curve.



t_opt = Inf; % Initialization
% t_path_all = zeros(size(P_mat,1),1);


% key is to minimize delta_t i.e the travel time it takes to move from the current
% point in k-space to the next. delta_g is also limited by the fact that
% gx_next cannot exceed the bounds on the maximum gradient magnitude set
% machine, and also on the max slew rates.

size_path = N;

paths = [9 8 10 7 6 3 4 5 2 1;9 10 8 7 6 3 4 5 2 1; 9 6 7 10 8 2 1 3 4 5]; 
g_path = cell(1,3);
t_path = zeros(1,3);

for index = 1:3
    index
    path_curr = paths(index,:);
      
    C = zeros(size_path, 3);
    
    for i=1:size_path
        C(i,1:2) = cell2mat(points_pos_U(path_curr(i)));
    end
      
    g_curr = [0, 0];
    
    %% Using Lustig's Optimization code
    [C,t,g,s,k,phi,sta,stb] = minTimeGradient(C,0,g_curr,[],g_max(1),slew_rates(1),4*10^-3,[],[]);
    %%
    
    g_path{index} = g;
    t_path(index) = t;
    
    if t <= t_opt
        t_opt = t;
        index_opt = index;
    end
end

% 
% gx_wave = zeros(1, N);
% gy_wave = zeros(1, N);
% t_path = zeros(1, N);
% 
% path_curr = P_mat(index_opt, :);
% g_curr = [0 0];
% k_curr = cell2mat(points_pos_U(path_curr(1)));
% for point_i = 2:N
%     k_next = cell2mat(points_pos_U(path_curr(point_i)));
%    [C,t,g,s,k,phi,sta,stb] = minTimeGradient(C,0,g_curr,[],g_max(1),slew_rates(1),4*10^-3,[],[]);
% end
% 

for index=1:3
    figure(index+1);
    subplot(2,1,1);
    plot(g_path{index}(:,1));
    title(['gx, time = ' num2str(t_path(index))]);
    subplot(2,1,2);
    plot(g_path{index}(:,2));
    title(['gy, path = ' num2str(paths(index,:))]);
end    
