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
initialize;
init_a;

frame = getframe(gcf); 
writeVideo(vidObj, frame);
%% Run simulation
for clock=1:clockmax
    XX=X+(dt/2)*vec_interp(u,X); % Euler step to midpoint
    ff=vec_spread(Force(XX),XX); % Force at midpoint
    ff(:,:,1) = ff(:,:,1) + f0;
    [u,uu]=fluid(u,ff); % Step Fluid Velocity
    X=X+dt*vec_interp(uu,XX); % full step using midpoint velocity

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
            speed_min = 0; speed_max = 1;
            avg_speed = mean(speed);

            cla
            scatter(particles(:,1), particles(:,2), 10, speed, 'filled', 'MarkerEdgeColor', 'none') 
            hold on
            plot(X(:,1),X(:,2),'ko')
            axis([0 L 0 L])
            axis square;
            axis manual
            colormap(jet)
            colorbar
            clim([speed_min speed_max])
            
            title(sprintf('Time = %.2f, Avg Speed = %.4f', clock * dt, avg_speed))
            drawnow
            hold off
            frame = getframe(gcf); 
            writeVideo(vidObj, frame);
        otherwise
            error('Invalid visualization mode.');
    end
end


close(vidObj);
disp(['Video saved as: ' videoName]);