%initialize.m
%% Initialize Parameters and special indices

L=1.0 % Box size
N=64 % Number of grid cells
h=L/N % Grid spacing
ip=[(2:N),1] % Grid index shifted left
im=[N,(1:(N-1))] % Grid index shifted right
Nb=ceil(pi*(L/2)/(h/2)) % Number of IB points
dtheta=2*pi/Nb % IB point spacing
kp=[(2:Nb),1] % IB index shifted left
km=[Nb,(1:(Nb-1))] % IB index shifted right
K=1 % Elastic stiffness
rho=1 % Fluid density
mu=0.01 % viscosity
tmax=10 % Run until time
dt=0.01 % Time step
clockmax=ceil(tmax/dt)

%% Initialize the visualization options
viz_option = "lagrangian particles"; % all options: "vorticity", "lagrangian particles"

%% Initialize boundary and velocity
f0=1e-3;
k=0:(Nb-1);
theta = k'*dtheta;
X = (L/2) + (L/4)*[cos(theta), sin(theta)];

u=zeros(N,N,2);
j1=0:(N-1); % Initialize fluid velocity as (0,sin(2*pi*x/L))
x=j1'*h;
u(j1+1,:,2)=sin(2*pi*x/L)*ones(1,N);

%% Initialize animation
switch viz_option
    case 'vorticity'
        vorticity=(u(ip,:,2)-u(im,:,2)-u(:,ip,1)+u(:,im,1))/(2*h);
        dvorticity=(max(max(vorticity))-min(min(vorticity)))/5;
        values= (-10*dvorticity):dvorticity:(10*dvorticity); % Get vorticity contours
        valminmax=[min(values),max(values)];
        xgrid=zeros(N,N);
        ygrid=zeros(N,N);
        for j=0:(N-1)
          xgrid(j+1,:)=j*h;
          ygrid(:,j+1)=j*h;
        end
        contour(xgrid,ygrid,vorticity,values)
        caxis(valminmax)

        set(gcf,'double','on')
        contour(xgrid,ygrid,vorticity,values)
        hold on
        plot(X(:,1),X(:,2),'ko')
        axis([0,L,0,L])
        caxis(valminmax)
        axis equal
        axis manual
        drawnow
        hold off
    case "lagrangian particles"
        % visualize by Lagrangian Particle Tracking
        np=40;
        [particles_x, particles_y] = meshgrid(linspace(0,L,np),linspace(0,L,np));
        particles = [particles_x(:), particles_y(:)];
        [X_grid, Y_grid] = meshgrid(linspace(0, L, N), linspace(0, L, N));
        particles_vel_x = interp2(X_grid, Y_grid, u(:,:,1), particles(:,1), particles(:,2), 'linear', 0);
        particles_vel_y = interp2(X_grid, Y_grid, u(:,:,2), particles(:,1), particles(:,2), 'linear', 0);
        particles_vel2 = [particles_vel_x, particles_vel_y];
        speed = sqrt(particles_vel2(:,1).^2 + particles_vel2(:,2).^2);
        max(speed)
        speed_min = 0; speed_max = 1;
        avg_speed = mean(speed);
        
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
    otherwise
        error('Invalid visualization mode.');
end