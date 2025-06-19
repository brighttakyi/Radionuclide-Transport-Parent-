%Solution of 1-D Transport Equation by Bright Takyi (10878402)
% Transport Equation; s_t = -v*s_x + D*s_xx - lambda
% Initial and Boundary conditions;
% s(x,0)=s0, s(0,t)=s1, ds/dx=0
%***********************************************************************************
% Initializing Parameters
 clear
 close all
 clc
L = 10000; % lenght of river [m]
T = 2*3600; % simulation time [sec]
%T = 172800; % simulation time (20 days)[s]
%oneday = 24*3600; % converting one-day into seconds
dt = 1;
dx = 20;
n = T/dt;
m = L/dx;
%m = 500; % number of subintervals of length of river
%n = 43200; % number of subintervals of simulation time
%n = 360; %1728; % number of subintervals of simulation time
%dx = L/m; % space step-size
%dt = %T/n; % time step-size
%oneday = 24*3600;
x = linspace(0,L,m+1); % space discretization
t = linspace(0,T,n+1); % time discretization
rx = dt/(2.*dx);
rxx = dt/(dx^2);
%halflife = 3600;
halflife = 3.8*24*3600; % halflife of Ra-222 in seconds
decay = log(2)/halflife; % decay rate of Ra-222
decayp = 0.0158;
save('decayp','decayp')
%decay = 0;
e = 1-decay*dt;
v = 3.*ones(m+1,1);% wave celerity [m/s]
%v1 = xlsread('C:\Users\G\Desktop\Velocity.xls','A4:B23');
%v=repmat(v1',m+1,1);
%v = 40.*ones(m+1,1); % wave celerity [m/s]                                        
%D1 = xlsread('C:\Users\G\Desktop\Diffusivity.xls','A4:B23');
%D=repmat(D1',m+1,1);% diffusion coefficient [m^3/s]
D = 0.003.*ones(m+1,1); % diffusion coefficient [m^3/s]
%D = 5.*ones(m+1,1);
% Setting up Q matrix, initial and boundary conditions
s = zeros(m+1,n+1);  
sa = 0; % initial concentration of radon-220 [mol/m^3]
sb = 100;% boundary condition at X = 0
%s(:,1) = 100;
s(1,:) = sb;
%sa = s(1);
s(2:end,1) = sa;
%legend('Analytical at distance = 1000 m')
%Q(:,1) = 850; % initial condition
%Q(1,:) = 850; % boundary condition
%Qa = Q(1);
 
% For loop 
for k = 2:n+1   % temporal loop
    a = zeros(m+1,1);
    b = zeros(m+1,1);
    d = zeros(m+1,1);
    h = zeros(m+1,1);
    A = zeros(m,m);
    %F = zeros(m+1,1);
    for i = 2:m+1   % spatial loop
          h(i) = e*s(i,k-1);
          a(i) = -0.5*(v(i-1)*rx+v(i)*rx)+0.5*(D(i-1)*rx+D(i)*rx);
          if i==m+1
              b(i) = 1+0.5*(v(i-1)*rx+v(i)*rx)+0.5*(D(i-1)+D(i)+D(i-1)+D(i))*rxx;
              d(i) = -0.5*(D(i-1)+D(i))*rxx;
          else   % for i < m+1
              b(i) = 1+0.5*(v(i+1)*rx+v(i)*rx)+0.5*(D(i+1)+D(i)+D(i-1)+D(i))*rxx;
 
              d(i) = -0.5*(D(i+1)+D(i))*rxx;
 
          end
          
        if i==2
            h(i) = h(i)-a(i)*s(1,k);
            %F(i) = s(i,k-1)-a(i)*sa;
            A(i-1,i-1) = b(i);
            A(i-1,i) = d(i);
        elseif i==m+1
            A(i-1,i-1) = b(i);
            A(i-1,i-2) = a(i)+d(i);
           % h(i) = s(i,k-1);
        else
            A(i-1,i-2) = a(i);
            A(i-1,i-1) = b(i);
            A(i-1,i) = d(i);
            %h(i) = s(i,k-1);
        
        end
    end
    
    
    
    s(2:end,k) = (A\h(2:end,1));
 
end
s1 = s; % concentrations to transfer to daughter
s = s./sb;
spN = s(51,:);
load spA
%sp2 = s(51,1:1001);
sR = spA(1:7201) - spN(1:7201);
save ('s1','s1')
%sR = spA - spN;
figure
plot(t(1:7201),spA(1:7201),'r','LineWidth',3)
hold on
plot(t(1:7201),spN(1:7201),'y','LineWidth',3)
xlabel('Time [s]')
ylabel('Relative Concentration [S/SO]')
%title('Graph of Radionuclide Concentration s/so [Bq]')
legend('Analytical','Numerical','Location','SouthEast')

%figure
%plot(t(1:7201),sR(1:7201),'k','LineWidth',3)
%xlabel('Time [s]')
%ylabel('Residual')
%title('Graph of Radionuclide Concentration Residual')
figure
%plot(t(1:7201),s(1,1:7201),'LineWidth',3);
plot(t(1:7201),s(51,1:7201),'LineWidth',3);
hold on
plot(t(1:7201),s(151,1:7201),'LineWidth',3);
plot(t(1:7201),s(301,1:7201),'LineWidth',3);
plot(t(1:7201),s(501,1:7201),'LineWidth',3);
xlabel('Time [s]')
ylabel('Relative Concentration [s/so]')
legend('x = 1000 m','x = 3000 m', 'x = 6000 m','x = 10000 m','location','SouthEast')
%legend('x = 0 m','x = 500 m', 'x = 1000 m', 'x = 3000 m','location','NorthEast')
 
 
figure
%plot(x,s(:,16),'LineWidth',3);
plot(x,s(:,601),'LineWidth',3);
hold on 
plot(x,s(:,1801),'LineWidth',3);
plot(x,s(:,3001),'LineWidth',3);
plot(x,s(:,7201),'LineWidth',3);
xlabel('Distance from initial section [m]')
ylabel('Relative Concentration [s/so]')
legend('After 600 s','After 1800 s', 'After 3000 s', 'After 7200 s','location','SouthWest')
 
 
%figure
%mesh(x,t,s','LineWidth',3)
%xlabel('Distance [m]')
%ylabel('Time [s]')
%zlabel('Radionuclide Concentration (s)[Bq]')
%title('A 3D Graph of Ra-222 concentration')

