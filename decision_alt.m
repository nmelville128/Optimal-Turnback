function [dec_alt, g_m, bank_0, bank_m, touchdown] = decision_alt(v_stall, v_ldmax, v_climb, thrust, ld_max, weight_max, landing_distance, takeoff_distance, wind, weight, runway, c_grad, t_reaction, bank_max, aoa_max, alt_buffer, plots)

%% Optimal Turnback

%% Variables (you need to change these)
% Environment variables
g = 32.174;
S = 1; %This cancels out, so the value does not matter
rho = 1; %This cancels out, so the value does not matter


%% Aerodynamic calculations
[cl0,cl1,cd0,cd1,cd2] = AeroForcesSetup(v_stall, v_ldmax, ld_max, weight_max, S, rho);

climb_angle = atan(c_grad);
L_climb = (weight - thrust*sin(climb_angle))/(cos(climb_angle));
cl_climb = 2*L_climb/(rho*S*v_climb^2);
a_climb = (cl_climb-cl0)/cl1;

%% Problem setup
% Number of states, control inputs, and collocation points
Nx = 6; % x, y, z, vx, vy, vz. X is departure direction, Y is right, Z is down.
Nu = 3; % Bank, aoa, and thrust
Nv = Nx + Nu;
No = 2; % extra states
N = 51; % How many collocation points
Nr = 6; % Number of collocation points dedicated to reaction time (aoa and bank are constant)


% Time
%tf = 20;
s_max = 120;
dt = 1/(N-Nr-1);

% Initial guess for Z0
% x0 = [0; 0;-300; 140; -20; -25];
% u0 = [1;0.5;0];
% xf = [0; 0; 0; -90; 0; 0];
% Z0 = [x0 + (xf-x0).*linspace(0,1,N);ones(Nu,N).*u0];
% Z0 = [s_max/2;0;reshape(Z0,[N*Nv,1])];
%Z0 = Zf;
Z0 = load('decision_alt_init.mat','Zf');
Z0 = Z0.Zf;

% Inequality constraints
Ni = N + 3;
Aineq = zeros(Ni,N*Nv+No);
bineq = zeros(Ni,1);
for i = 1:N
   % Aineq(i,2)=-1;
    Aineq(i,Nv*(i-1)+5)=1;
    bineq(i,1)=-alt_buffer;
end
Aineq(Ni-2,end-5)=1;
bineq(Ni-2)=1;
Aineq(Ni-1,end-8)=1;
bineq(Ni-1)=runway;
Aineq(Ni,end-8)=-1;
bineq(Ni)=-landing_distance;

% Equality Constraints
Ne = 9; %19
%Ne = 9; %Number of equality constraints. 6 for a TCCC. More for a turnback.
Aeq = zeros(Ne + 2*Nr,N*Nv+No);

Aeq(1:Nx,1+No:No+Nx) = zeros(Nx);
Aeq(1,2) = 1/c_grad;
Aeq(1,No+1) = 1;
Aeq(2,No+2) = 1;
Aeq(3,No+3) = 1;
Aeq(3,2) = -1;
Aeq(4,No+4) = 1;
Aeq(5,No+5) = 1;
Aeq(6,No+6) = 1;
for i = 1:Nr
Aeq(Ne + 2*i - 1,No+i*Nv-2) = 1;
Aeq(Ne + 2*i,No+i*Nv-1) = 1;
end

Aeq(7,end-7) = 1;  %y
Aeq(8,end-4) = 1; %vy
Aeq(9,end-3) = 1; %vz
%Aeq(4+Nx,end-8) = 1; 
%Aeq(5+Nx,end-5) = 1;
% Aeq(8,end - 5) = 1;
% Aeq(9,end - 3) = 1;

beq = [takeoff_distance;0;0;v_climb*cos(atan(c_grad));-1;-v_climb*sin(atan(c_grad));0;0;0;repmat([0,a_climb]',[Nr,1])];
%beq = x0;

% Upper and Lower Bounds
LB = [20;-5000;repmat(-[100000,100000,100000,1000,1000,1000,bank_max,0.0,0]',[N,1])];
UB = [s_max;-0;repmat([100000,100000,100000,1000,1000,1000,bank_max,aoa_max,0]',[N,1])];
%UB(end-8)=10000;

% Pack up variables for easier readability
P.Nx = Nx;
P.Nu = Nu;
P.Nv = Nv;
P.N = N;
%P.tf = tf;
P.dt = dt;
P.wind = wind;
P.rho = rho;
P.g = g;
P.cl0 = cl0;
P.cl1 = cl1;
P.cd0 = cd0;
P.cd1 = cd1;
P.cd2 = cd2;
P.S = S;
P.m = weight;
P.vf = v_ldmax;
P.ld_max = ld_max;
P.s_max = s_max;
P.v_stall = v_stall;
P.No = No;
P.Nr = Nr;
P.Tr = t_reaction;

% Objective Function
obj = @(Z)( myObjective(Z,P) );

% Nonlinear constraint
nonlcon = @(Z) ( myConstraint(Z,P) );

% Solver options
options = optimoptions('fmincon','MaxFunctionEvaluations',N*0.25e+05,'MaxIterations',N*1.000000e+04);

%%
% The key step
[Zf,fval,exitflag] = fmincon(obj,Z0,Aineq,bineq,Aeq,beq,LB,UB,nonlcon,options); 
%Zf
beep
toc
exitflag

fval
if exitflag <= 0
    dec_alt = nan;

    g_m = nan;

    bank_0 = nan;

    bank_m = nan;

    touchdown = nan;


else
    save('decision_alt_init.mat','Zf');

    %% Post processing
    %c = -myObjective(Zf,P)
    s=Zf(1);
    %alt_min = Zf(2)
    Z = reshape(Zf(No+1:end),[Nv,N]);
    ts = linspace(0,1,N-Nr+1);
    ts = ts*s;
    ts = [linspace(0,t_reaction,Nr),ts(2:end)+t_reaction];

    x = Z(1,:);
    y = Z(2,:);
    z = Z(3,:);
    vx = Z(4,:);
    vy = Z(5,:);
    vz = Z(6,:);
    bank = Z(7,:);
    aoa = Z(8,:);
    thrust = Z(9,:);
    v = sqrt(vx.^2 + vy.^2 + vz.^2);
    h = sqrt(x.^2 + y.^2);
    cl = cl0 + cl1.*aoa;
    cd = cd0 + cd1.*aoa + cd2.*aoa.^2;
    L = 0.5*rho*v.^2.*cl*S;
    %D = 0.5*rho*v.^2.*cd*S;
    gs = L/weight;

    dec_alt = fval;

    g_m = max(gs);

    bank_0 = bank(Nr+1);

    bank_m = max([max(bank),-min(bank)]);

    touchdown = x(end);

    if plots
        % Plotting
        x = Z(1,:);
        y = Z(2,:);
        z = Z(3,:);
        vx = Z(4,:);
        vy = Z(5,:);
        vz = Z(6,:);
        bank = Z(7,:);
        aoa = Z(8,:);
        thrust = Z(9,:);
        v = sqrt(vx.^2 + vy.^2 + vz.^2);
        h = sqrt(x.^2 + y.^2);

        tiledlayout(2,3)
        nexttile
        figure(1)
        scatter(ts,-z);
        xlabel('time (s)')
        ylabel('altitude (ft)')

        nexttile
        scatter(x,y);
        xlabel('x (ft)')
        ylabel('y (ft)')

        nexttile
        plot(ts,v)
        xlabel('time (s)')
        ylabel('velocity (KIAS)')

        nexttile
        plot(ts(2:end-1),aoa(2:end-1))
        xlabel('time (s)')
        ylabel('aoa')

        %g-forces
        nexttile
        cl = cl0 + cl1.*aoa;
        cd = cd0 + cd1.*aoa + cd2.*aoa.^2;
        L = 0.5*rho*v.^2.*cl*S;
        %D = 0.5*rho*v.^2.*cd*S;
        gs = L/weight;
        plot(ts(2:end-1),gs(2:end-1));
        xlabel('time (s)')
        ylabel('g-force')

        nexttile
        plot(ts(2:end-1),bank(2:end-1))
        xlabel('time (s)')
        ylabel('bank (s)')
    end
end

end
% % Send to Simulink for visualization
% %control = [bank;aoa;thrust];
% x0 = [0;0;0;100;0;0];
% control = repmat([0,0.5,0],N,1);
% control = timeseries(control,ts);

%% Functions
% Objective Function
function cost = myObjective(Z,P)
ld_max = P.ld_max;
%cost = Z(end-8) + ld_max*Z(end-6); %Use this for TCCC
cost = -Z(2); %Use this for turnback

end

% Dynamics
function [cineq,ceq] = myConstraint(Z, P)
% Unpack
N = P.N;
Nx = P.Nx;
Nu = P.Nu;
Nv = P.Nv;
%tf = P.tf;
dt = P.dt;
wind = P.wind;
rho = P.rho;
g = P.g;
cl0 = P.cl0;
cl1 = P.cl1;
cd0 = P.cd0;
cd1 = P.cd1;
cd2 = P.cd2;
S = P.S;
m = P.m;
vf = P.vf;
s_max = P.s_max;
v_stall = P.v_stall;
No = P.No;
Nr = P.Nr;
Tr = P.Tr;

% Constants
kn_fps = 1.6878;
fps_kn = 1/kn_fps;

% Intermediate
s = Z(1);
Z = Z(No+1:end);

bank = Z(Nx + 1:Nv:end,1);
aoa = Z(Nx + 2:Nv:end,1);
thrust = Z(Nx + 3:Nv:end,1);

x = Z(1:Nv:end,1);
y = Z(2:Nv:end,1);
z = Z(3:Nv:end,1);
vx = Z(4:Nv:end,1);
vy = Z(5:Nv:end,1);
vz = Z(6:Nv:end,1);

v = sqrt(vx.^2 + vy.^2 + vz.^2);
cl = cl0 + cl1*aoa;
cd = cd0 + cd1*aoa + cd2*aoa.^2;
L = 0.5*rho*v.^2.*cl*S;
D = 0.5*rho*v.^2.*cd*S;

% State derivatives
dx = kn_fps*(vx + wind(1));
dy = kn_fps*(vy + wind(2));
dz = kn_fps*(vz + wind(3));

% body frame
ax_thrust_b = (1/m).*thrust * g;

% wind frame
ax_w = cosd(20 * aoa - 5 ).*ax_thrust_b - (1/m).*D * g;
az_w = - sind(20 * aoa - 5 ).*ax_thrust_b - (1/m).*L * g;

% inertial frame
heading = (180/pi)*atan2(vy,vx);
climb_angle = (180/pi)*atan2(vz,sqrt(vx.^2 + vy.^2));


%az = cosd(bank).*cosd(climb_angle).*ax_w -sind(climb_angle).*az_w;
%ay = (sind(bank).*sind(climb_angle).*cosd(heading) - cosd(bank).*sind(heading)).*ax_w + sind(bank).*cosd(climb_angle).*az_w;
%ax = (cosd(bank).*sind(climb_angle).*cosd(heading) + sind(bank).*sind(heading)).*ax_w + cosd(bank).*cosd(climb_angle)*az_w + g;

az = fps_kn*(sind(climb_angle).*ax_w + cosd(climb_angle).*cosd(bank).*az_w + g);

ay = fps_kn*(sind(heading).*cosd(climb_angle).*ax_w + (-sind(heading).*sind(climb_angle).*cosd(bank) + cosd(heading).*sind(bank)).*az_w);
%ay = 0.*ax_w;

ax = fps_kn*(cosd(heading).*cosd(climb_angle).*ax_w + (-cosd(heading).*sind(climb_angle).*cosd(bank) - sind(heading).*sind(bank)).*az_w);


% reshaping
Xdot = [dx, dy, dz, ax, ay, az];
X = [x, y, z, vx, vy, vz];

Xdot = reshape(Xdot',[N*Nx,1]);
X = reshape(X',[N*Nx,1]);

ceq = X(Nx + 1:end) - X(1:end - Nx) - 0.5*s*dt*Xdot(Nx + 1:end) - 0.5*s*dt*Xdot(1:end - Nx);
ceq(1:Nx*Nr-Nx) = X(Nx + 1:Nx*Nr) - X(1:Nx*Nr - Nx) - 0.5*(Tr/Nr)*Xdot(Nx + 1:Nx*Nr) - 0.5*(Tr/Nr)*Xdot(1:Nx*Nr - Nx);

%ceq = X(Nx + 1:end) - X(1:end - Nx) - 0.5*dt*Xdot(Nx + 1:end) - 0.5*dt*Xdot(1:end - Nx);
%cineq = zeros(N,1);
cineq = [];
%cineq=  - (X(end-2))^2 - (X(end-1))^2 + (X(end))^2 + v_stall^2;

end