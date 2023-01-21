% Matlab code to simulate some pouplar friction models
% Part of this code was adapted from the work of Auralius Manurung. 
% Check: LuGre friction model in MATLAB (https://github.com/auralius/LuGre), GitHub.
clear all;

% Initialization of parameters
M=1; % Mass in kg
K=20; % N/m
tspan=[0 5];  % This is the time range of the simulation and the solver pick its own sapmling rate;
x0=[0 0]; % initial condition
miu_c=0.1; % coulomb friction
miu_v=0.5; % coulomb frictio
v_b=0.1 % the velocity of the belt in m/s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coulomb friction model
% Numerical integration with ode45 and ode23s
[t1,x1]=ode23s(@(t,x) odefunc_coulomb(t,x,M,miu_c,K,v_b),tspan,x0);
F_f1 = nan(size(t1));
for i=1:size(t1)
F_n=M*9.81; % normal force in newton
F_f1(i)=miu_c*F_n* sign(x1(i,2));
end
% Plot the displacement and relative velocity of the body
figure(1)
subplot(2,1,1)
plot(t1,v_b*t1,'r','LineWidth',2)
hold on
plot(t1, x1(:,1),'b','LineWidth',2)
ylabel('Position (m)')
legend('reference displacement - x_{ref}', 'displacement- x', 'Location','best','interpreter','latex')
xlabel('Time (s)')
title('Simulation of stick-slip motion based on Coulomb friction model') 
grid on

subplot(2,1,2)
yyaxis left
plot(t1,F_f1,'LineWidth',2)
ylabel('Friction Force (N)')
ylim([0 1.2]);
hold on
yyaxis right
plot(t1, x1(:,2),'LineWidth',2);
ylabel('$\frac{dx}{dt}$ (m/s)', 'interpreter','latex')
legend('Friction force - F_f', 'Sliding velocity - x_{dot}', 'Location','best','interpreter','latex')
xlabel('Time (s)')
ylim([0 0.4]);
title('Friction force and sliding velocity')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coulomb + visocus friction model
% Numerical integration with ode45
[t2,x2]=ode23s(@(t,x) odefunc_coulomb_viscous(t,x,M,miu_c,miu_v,K,v_b),tspan,x0);

F_f2 = nan(size(t2));
for i=1:size(t2)
F_n=M*9.81; % normal force in newton
F_f2(i)=miu_c*F_n* sign(x2(i,2)) +(miu_v*x2(i,2));
end

% Plot the displacement and relative velocity of the body
figure(2)
subplot(2,1,1) 
plot(t2,v_b*t2,'r','LineWidth',2)
hold on
plot(t2, x2(:,1),'b','LineWidth',2)
ylabel('Position (m)')
legend('reference displacement - x_{ref}', 'displacement- x', 'Location','best','interpreter','latex')
xlabel('Time (s)')
title('Simulation of stick-slip motion based on Coulomb+viscous friction model') 
grid on

subplot(2,1,2)
yyaxis left
plot(t2,F_f2,'LineWidth',2)
ylabel('Friction Force (N)')
ylim([0 1.5]);
hold on
yyaxis right
plot(t2, x2(:,2),'LineWidth',2);
ylabel('$\frac{dx}{dt}$ (m/s)', 'interpreter','latex')
legend('Friction force - F_f', 'Sliding velocity - x_{dot}', 'Location','best','interpreter','latex')
xlabel('Time (s)')
ylim([0 0.4]);
title('Friction force and sliding velocity')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LuGre friction model

sigma_0 = 1e5;
sigma_1  = sqrt(1e5);
sigma_2  = 0.4;
Fc = miu_c*9.81;
Fs = Fc+0.5;
vs = 0.001;
%vs = 0.001;

time_span = [0 5];
q_initial = [0 0 0];
    
% Use ode23s or ode45
[t3, q] = ode23s(@odefunc_LuGre, tspan, q_initial, [], ...
                        M,K,v_b, Fs, Fc, sigma_0, sigma_1, sigma_2, vs);   

% The problem is, there is no clean way to pass out other results with
% the built-in solver. We have to recompute the friction force.
[~,zdot,F_f3] = odefunc_LuGre(t3, q', M,K,v_b, Fs, Fc, sigma_0, sigma_1, sigma_2, vs);

figure(4) 

subplot(2,1,1)

plot(t3,v_b*t3,'r','LineWidth',2)
hold on
plot(t3, q(:,1),'b','LineWidth',2)
ylabel('Position (m)')
legend('reference displacement - x_{ref}', 'displacement- x', 'Location','best','interpreter','latex')
xlabel('Time (s)')
title('Simulation of stick-slip motion based on Coulomb + viscous friction model')
grid on

subplot(2,1,2)
hold on
yyaxis left
plot(t3, F_f3, 'LineWidth',2);
ylabel('Friction Force (N)')
ylim([0 1.6]);
hold on
yyaxis right
plot(t3, q(:,2),'LineWidth',2);
ylabel('$\frac{dx}{dt}$ (m/s)', 'interpreter','latex')
legend('Friction force - F_f', 'Sliding velocity - x_{dot}', 'Location','best','interpreter','latex')
xlabel('Time (s)')
ylim([0 0.5]);
title('Friction force and sliding velocity')
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure(5)
figure('DefaultAxesFontSize',20)
set(gcf,'color','w');
set(gca,'FontWeight','bold');
extraInputs = {'interpreter','latex','fontsize',25,'FontWeight', 'bold'};
set(gcf,'color','w');
%subplot(2,1,1)
plot(t2,v_b*t2,'m','LineWidth',3)
hold on
plot(t1, x1(:,1),'r','LineWidth',3)
plot(t2, x2(:,1),'k','LineWidth',3)
plot(t3, q(:,1),'b','LineWidth',3)
ylabel('Displacement (m)',extraInputs{:})
legend('Reference','Coulomb', 'Coulomb + viscous','LuGre')
xlabel('Time (s)',extraInputs{:})
%title('Simulation of stick-slip motion based on different friction model',extraInputs{:})
ylim([0 0.65])
grid on

%subplot(2,1,2)
%figure(7)
figure('DefaultAxesFontSize',20)
set(gcf,'color','w');
set(gca,'FontWeight','bold');
extraInputs = {'interpreter','latex','fontsize',25,'FontWeight', 'bold'};

plot(t1, x1(:,2),'r','LineWidth',3)
hold on
plot(t2, x2(:,2),'k','LineWidth',3)
plot(t3, q(:,2),'b','LineWidth',3)
ylabel('Velocity (m/s)',extraInputs{:})
legend('Coulomb', 'Coulomb + viscous','LuGre')
xlabel('Time (s)',extraInputs{:})
%title('Simulation of stick-slip motion based on different friction model',extraInputs{:})
ylim([0 0.28])
grid on

%figure(8)
figure('DefaultAxesFontSize',20)
set(gcf,'color','w');
set(gca,'FontWeight','bold');
extraInputs = {'interpreter','latex','fontsize',25,'FontWeight', 'bold'};
set(gcf,'color','w');
%subplot(2,1,1)
plot(x1(:,2),F_f1,'r','LineWidth',3)
hold on
plot(x2(:,2),F_f2,'k','LineWidth',3)
plot(q(:,2),F_f3,'b','LineWidth',3)
ylabel('Force (N)',extraInputs{:})
legend('Coulomb', 'Coulomb + viscous','LuGre')
xlabel('Velocity (m/s)',extraInputs{:})
%title('Simulation of stick-slip motion based on different friction model',extraInputs{:})
%ylim([0 0.28])
grid on

%subplot(2,1,2)
%figure(10)
figure('DefaultAxesFontSize',20)
set(gcf,'color','w');
set(gca,'FontWeight','bold');
extraInputs = {'interpreter','latex','fontsize',25,'FontWeight', 'bold'};

plot(t1,F_f1,'r','LineWidth',3)
hold on
plot(t2,F_f2,'k','LineWidth',3)
plot(t3,F_f3,'b','LineWidth',3)
ylabel('Friction Force (N)',extraInputs{:})
legend('Coulomb', 'Coulomb + viscous','LuGre')
xlabel('Time (s)',extraInputs{:})
%title('Frictional force with different models',extraInputs{:})
ylim([0 1.8])
grid on