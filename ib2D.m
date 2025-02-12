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
global viz_gap;
initialize;
rod=s;
init_a;

frame = getframe(gcf); 
writeVideo(vidObj, frame);
constantFlow = false;  % 开启常数流场测试模式

%% Run simulation
for clock=1:clockmax
    if ~constantFlow
        % 这里是原来的刚棒和外力更新代码
    % 1. 计算theta_mid, 使用当前theta和omega半步前进
    theta_mid = theta + 0.5*dt*omega;
    global Z_mid;
    Z_mid = zeros(Nb,2);
    Z_mid(:,1) = Px + (rod - Lp)*cos(theta_mid);
    Z_mid(:,2) = Py + (rod - Lp)*sin(theta_mid);

    XX = X + (dt/2)*vec_interp(u,X); % 用中点速度更新X的中点位置
    F_mid = Force(XX);              % = K(Z_mid - XX)
    ff = vec_spread(F_mid,XX);
    ff(:,:,1) = ff(:,:,1) + f0;
    [u,uu] = fluid(u,ff);
    X_new = X + dt*vec_interp(uu,XX); % 全步X更新（使用中点u）

    % 用X_new和F_mid计算中点扭矩tau_mid
    Rx = XX(:,1) - Px;
    Ry = XX(:,2) - Py;
    tau_mid = sum((Rx.*F_mid(:,2) - Ry.*F_mid(:,1)) * ds);

    % 全步更新omega和theta
    omega_new = omega + (dt/I0)*tau_mid;
    theta_new = theta + dt*omega_new;

    % 更新X,theta,omega用于下个时间步
    X = X_new;
    theta = theta_new;
    omega = omega_new;
    else
        % constantFlow 模式：保持流场不变
        uu = u;  % 直接使用常数流场
    end

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