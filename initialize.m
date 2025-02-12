%initialize.m
%% Initialize Parameters and special indices

L=1.0 % Box size                                        % m
N=64 % Number of grid cells
h=L/N % Grid spacing                                    % m
ip=[(2:N),1] % Grid index shifted left
im=[N,(1:(N-1))] % Grid index shifted right
K = 1e6;                                                % N/m²
rho=1; % Fluid density                                  % kg/m³
mu=0.01; % viscosity                                    % kg/(m·s)
tmax=10; % Run until time                               % s
dt=2e-4; % Time step                                    % s
clockmax=ceil(tmax/dt);

%% 数据记录参数
viz_gap = 0.0025; % 可视化时间间隔 (s)
record_interval = round(viz_gap / dt); % 数据记录步长间隔

%% Initialize the visualization options
viz_option = "lagrangian particles"; % all options: "vorticity", "lagrangian particles"

%% Initialize boundary and velocity
%f0=1e-3;
%k=0:(Nb-1);
%theta = k'*dtheta;
%X = (L/2) + (L/4)*[cos(theta), sin(theta)];

% Assume the 假设线段长度 L0，pivot点位于 Lp
%f0=1e2;
f0=2e0;
L0 = L/2;
Px=L/2;
Py=L/2;
Nb = ceil(L0/(h/2));
ds = L0/(Nb-1);
s = (0:(Nb-1))*ds;
Lp = (Nb-1)/2*ds; %L0/2;
m0 = 5e0;            % 线密度(可根据需要调整)
theta0 = pi/3;
X = zeros(Nb,2);
X(:,1) = Px + (s - Lp)*cos(theta0);
X(:,2) = Py + (s - Lp)*sin(theta0);
pivot = [Px;Py];
% I0 = m0 * integral_0^L0 (s - Lp)^2 ds = m0 * [ ( (L0 - Lp)^3 + (Lp)^3 ) / 3 ]
I0 = m0*((L0 - Lp)^3 + Lp^3)/3; 
theta = theta0;   % 初始角度
omega = 0;   % 初始角速度


u=zeros(N,N,2);
U0 = 0;  % 设定水平向右的常数速度，比如1 m/s，可根据需要调整
u(:,:,1) = U0;  % x方向速度全部设为 U0
u(:,:,2) = 0;   % y方向速度为 0
%j1=0:(N-1); % Initialize fluid velocity as (0,sin(2*pi*x/L))
%x=j1'*h;
%u(j1+1,:,2)=sin(2*pi*x/L)*ones(1,N);

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
        plot(pivot(1), pivot(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2)
        axis([0 L 0 L])
        axis square;
        axis manual
        
        %colormap(jet)
        colormap(parula)
        colorbar
        clim([speed_min speed_max])
        title(sprintf('Time = %.2f, Avg Speed = %.4f', clock * dt, avg_speed))
        drawnow
        hold off
    otherwise
        error('Invalid visualization mode.');
end