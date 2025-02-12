% ib2D.m
% This script is the main program.
% Original Code by Charlie Peskin:
% https://www.math.nyu.edu/faculty/peskin/ib_lecture_notes/index.html
% Vectorized and commented by Tristan Goodwill,2019.4

% This new version of the code is written by Yifei Zhu for the NYU ib course.
% add visualization options and video generation

%% Video setup
currentTime = datetime('now');
postfix = sprintf('_%02d%02d-%02d%02d%02d', currentTime.Month, currentTime.Day, currentTime.Hour, currentTime.Minute, floor(currentTime.Second));
videoName = ['./video/ib_final', postfix, '.avi'];
vidObj = VideoWriter(videoName, 'Motion JPEG AVI');
vidObj.FrameRate = 40;
vidObj.Quality = 95;
open(vidObj);

%% Initialize simulation
global dt Nb N h rho mu ip im a;
global kp km dtheta K;
global f0;
global viz_options particles;
global theta omega I0 L0 Lp m0 s ds Px Py pivot;
global analysis_data viz_gap record_interval vorticity_snapshots is_upper is_lower X_grid Y_grid;
initialize;
rod=s;
init_a;

frame = getframe(gcf); 
writeVideo(vidObj, frame);
constantFlow = false;  % 开启常数流场测试模式
record = true;

num_records = floor(clockmax/record_interval);
t_rec      = zeros(num_records, 1);   % 记录时间
theta_rec  = zeros(num_records, 1);   % 记录角度 theta (rad)
omega_rec  = zeros(num_records, 1);   % 记录角速度 omega (rad/s)
tau_rec    = zeros(num_records, 1);   % 记录净扭矩 tau (N.m)
error_rec  = zeros(num_records, 1);   % 记录刚性误差（例如：目标与实际位置的最大差值）
record_index = 0;  % 数据记录的索引

% Steady-state detection parameters
tol_omega = 1e-4;         % 角速度阈值 (rad/s)
tol_tau   = 1e-4;         % 净转矩阈值 (N.m)
steady_count = 0;         % 连续满足稳态条件的步数计数器
required_steady_steps = 100;  % 连续满足要求的步数（可根据问题尺度调整）

%% Run simulation
for clock=1:clockmax
% prelimnary substep from n to n+ 1/2
    XX = X + (dt/2)*vec_interp(u,X);

    theta_mid = theta + 0.5*dt*omega;

    global Z_mid;
    Z_mid = zeros(Nb,2);
    Z_mid(:,1) = Px + (rod - Lp)*cos(theta_mid);
    Z_mid(:,2) = Py + (rod - Lp)*sin(theta_mid);

    F_mid = Force(XX);              % = K(Z_mid - XX)
    ff = vec_spread(F_mid,XX);
    ff(:,:,1) = ff(:,:,1) + f0;

    tau_mid = - sum( (rod - Lp)' .* ( -sin(theta_mid) * F_mid(:,1) + cos(theta_mid) * F_mid(:,2) ) * ds );

    omega_mid = omega + (dt/(2*I0))*tau_mid;

    % solve for fluid field
    [u,uu] = fluid(u,ff);

    % full timestep update
    X = X + dt*vec_interp(uu,XX);
    theta = theta + dt*omega_mid;
    omega = omega + dt/I0*tau_mid;

    % 数据记录（每隔 record_interval 步记录一次）
    if record && mod(clock, record_interval) == 0
        record_index = record_index + 1;
        t_rec(record_index)     = clock * dt;
        theta_rec(record_index) = theta;
        omega_rec(record_index) = omega;
        tau_rec(record_index)   = tau_mid;  % 此处记录半步计算的净转矩
        % 记录刚性误差：目标位置 Z_mid 与中间位置 XX 的最大差值
        error_rec(record_index) = norm(Z_mid - XX, inf);

        %% debug
        analysis_data.time(record_index)   = clock * dt;
        analysis_data.angle(record_index)  = theta;
        analysis_data.torque(record_index) = tau_mid;
        % 当前流场速度幅值
        flow_speed = sqrt(u(:,:,1).^2 + u(:,:,2).^2);
        
        % 上下半区速度差异
        analysis_data.upper_speed(record_index) = mean(flow_speed(is_upper), 'all');
        analysis_data.lower_speed(record_index) = mean(flow_speed(is_lower), 'all');
        
        % 压力场差异
        if exist('p','var') % 若压力场已计算
            analysis_data.pressure_diff(record_index) = mean(p(is_upper), 'all') - mean(p(is_lower), 'all');
        end
        
        % 涡量场特征提取
        [du_dy, du_dx] = gradient(u(:,:,1), h, h);
        [dv_dy, dv_dx] = gradient(u(:,:,2), h, h);
        vorticity = dv_dx - du_dy;
       
        % 提取棒后方涡量极值
        vort_region = vorticity(N/2:N, 1:N/2);
        analysis_data.vorticity_peak(record_index) = max(abs(vort_region(:)));

        %-----------------------------%
        % 保存涡量场快照（用于后续 GIF 动画）
        %-----------------------------%
        % 这里采用与可视化中类似的计算方法：
        vort_field = (u(ip,:,2) - u(im,:,2) - u(:,ip,1) + u(:,im,1)) / (2*h);
        vorticity_snapshots{record_index} = vort_field;
    end

    % Steady-state check: 当角速度和净转矩均低于阈值时，增加计数，否则重置
    if abs(omega) < tol_omega && abs(tau_mid) < tol_tau
        steady_count = steady_count + 1;
    else
        steady_count = 0;
    end

    % 如果连续满足稳态条件的步数达到要求，则退出仿真
    if steady_count >= required_steady_steps
        fprintf('Steady state detected at time %.4f s (after %d steps)\n', clock*dt, clock);
        break;
    end

    % visualize only at fixed interval of time, not every timestep
    if mod(clock*dt/viz_gap, 1) == 0
        switch viz_option
            case 'vorticity'
                vorticity=(u(ip,:,2)-u(im,:,2)-u(:,ip,1)+u(:,im,1))/(2*h);
                contour(xgrid,ygrid,vorticity,values)
                hold on
                plot(X(:,1),X(:,2),'ko')
                axis([0,L,0,L])
                caxis(valminmax)
                axis equal
                axis manual
                drawnow
                hold off
                frame = getframe(gcf); 
                writeVideo(vidObj, frame);
            case "lagrangian particles"
                particle_vel1 = particle_interp(u, particles);
                particles_temp = particles + 0.5*dt*particle_vel1;
                particle_vel2 = particle_interp(uu, particles_temp);
                particles = particles + dt*particle_vel2;
                particles = mod(particles, L);
                speed = sqrt(particle_vel2(:,1).^2 + particle_vel2(:,2).^2);
                speed_min = 0; speed_max = 2;
                avg_speed = mean(speed);
                grid_speed = sqrt(u(:, :, 1).^2 + u(:, :, 2).^2);
                max_grid_speed = max(grid_speed(:));
                CFL_val = (max_grid_speed * dt) / h;
    
                cla
                scatter(particles(:,1), particles(:,2), 10, speed, 'filled', 'MarkerEdgeColor', 'none') 
                hold on
                plot(X(:,1),X(:,2),'ko')

                % Add force visualization
                %quiver(X(:,1), X(:,2), F_mid(:,1), F_mid(:,2), 4.5, 'g-', 'LineWidth', 1)
                % 根据杆上参数 s 和给定的 pivot 参数位置 Lp 分割杆的两部分
                front_indices = find(rod < Lp);   % pivot 后侧：s 值小于 Lp 的点
                back_indices = find(rod >= Lp); % pivot 前侧：s 值大于等于 Lp 的点
                % 计算前后两段的合力
                front_force = sum(F_mid(front_indices,:), 1);
                back_force = sum(F_mid(back_indices,:), 1);
                % 计算前后两段的力学中心
                front_center = mean(X(front_indices,:), 1);
                back_center = mean(X(back_indices,:), 1);
                % 在力学中心绘制合力箭头
                scale_factor = 0.002*L;
                quiver(front_center(1), front_center(2), front_force(1), front_force(2), scale_factor, 'b-', 'LineWidth', 2)
                quiver(back_center(1), back_center(2), back_force(1), back_force(2), scale_factor, 'g-', 'LineWidth', 2)

                plot(pivot(1), pivot(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2)
                axis([0 L 0 L])
                axis square;
                axis manual
                %colormap(jet)
                colormap(parula)
                colorbar
                clim([speed_min speed_max])
                
                title(sprintf('Time = %.2f, Avg Speed = %.2f, CFL = %.1f', clock * dt, avg_speed, CFL_val))
                drawnow
                hold off
                frame = getframe(gcf); 
                writeVideo(vidObj, frame);
            otherwise
                error('Invalid visualization mode.');
        end
    end
end


close(vidObj);
disp(['Video saved as: ' videoName]);



% Data analysis: Plot recorded data after simulation
figure;
subplot(2,2,1);
plot(t_rec, theta_rec, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\theta (rad)');
title('Rod Angle vs. Time');
grid on;

subplot(2,2,2);
plot(t_rec, omega_rec, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\omega (rad/s)');
title('Angular Velocity vs. Time');
grid on;

subplot(2,2,3);
plot(t_rec, tau_rec, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\tau (N.m)');
title('Net Torque vs. Time');
grid on;

subplot(2,2,4);
plot(t_rec, error_rec, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Max Rigidity Error (m)');
title('Max Rigidity Error vs. Time');
grid on;

fprintf('Final stable rod angle: %.4f degrees\n', rad2deg(theta_rec(end)));


%% 新增物理量分析
figure;
% 速度不对称性演化
subplot(2,2,1);
plot(analysis_data.time, analysis_data.upper_speed, 'r-', ...
     analysis_data.time, analysis_data.lower_speed, 'b-');
xlabel('Time (s)'); ylabel('Avg Speed (m/s)');
legend('Upper Half', 'Lower Half');
title('Flow Speed Asymmetry');

% 压力差与角度关系
subplot(2,2,2);
scatter(rad2deg(analysis_data.angle), analysis_data.pressure_diff, [], analysis_data.time, 'filled');
xlabel('Angle (deg)'); ylabel('Upper-Lower Pressure Diff (Pa)');
colorbar; title('Pressure Difference vs Angle');

% 涡量峰值与扭矩关系
subplot(2,2,3);
scatter(analysis_data.torque, analysis_data.vorticity_peak, [], analysis_data.time, 'filled');
xlabel('Torque (N·m)'); ylabel('Vorticity Peak (1/s)');
title('Vorticity vs Torque Correlation');

% 角度-扭矩相位图
subplot(2,2,4);
plot(rad2deg(analysis_data.angle), analysis_data.torque, 'k-o');
xlabel('Angle (deg)'); ylabel('Torque (N·m)');
hold on; 
yline(0, '--r'); 
plot(rad2deg(theta_rec(end)), 0, 'ro', 'MarkerSize', 10);
title('Torque-Angle Phase Diagram');

%% 涡量场动态可视化（生成GIF）
% vort_min = min(vorticity(:));
% vort_max = max(vorticity(:));
% %%
% figure;
% for n = 1:10:length(analysis_data.time)
%     pcolor(X_grid, Y_grid, vorticity);
%     shading interp; colormap(jet); caxis([vort_min vort_max]);
%     hold on;
%     plot(X(:,1), X(:,2), 'w-', 'LineWidth', 2);
%     hold off;
%     title(sprintf('Vorticity at t=%.2fs', analysis_data.time(n)));
%     drawnow;
% 
%     % 保存为GIF
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [A, map] = rgb2ind(im, 256);
%     if n == 1
%         imwrite(A, map, 'vorticity_evolution.gif', 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
%     else
%         imwrite(A, map, 'vorticity_evolution.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
%     end
% end

% 计算所有快照中的涡量全局最小值和最大值（用于固定色标范围）
vort_min = inf;
vort_max = -inf;
for k = 1:record_index
    current_vort = vorticity_snapshots{k};
    vort_min = min(vort_min, min(current_vort(:)));
    vort_max = max(vort_max, max(current_vort(:)));
end

% 生成 GIF 动画
figure;
gif_filename = 'vorticity_evolution.gif';
for n = 1:10:record_index  % 例如每隔10个记录点生成一帧，可根据需要调整
    pcolor(X_grid, Y_grid, vorticity_snapshots{n});
    shading interp;
    colormap(jet);
    caxis([vort_min vort_max]);
    hold on;
    % 如果希望在 GIF 中显示刚棒位置，也可以绘制 X 或其它相关数据
    plot(X(:,1), X(:,2), 'w-', 'LineWidth', 2);
    hold off;
    title(sprintf('Vorticity at t = %.2fs', analysis_data.time(n)));
    drawnow;

    % 捕获当前帧并写入 GIF 文件
    frame = getframe(gcf);
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256);
    if n == 1
        imwrite(A, map, gif_filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
    else
        imwrite(A, map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end