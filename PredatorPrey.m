%% Predator-Prey System

%system of first order differentials, x(t) denotes the population
%of the prey species, y(t) denotes the population of predator species

%a,c,alpha,gamma,r -> all positive constants

% presentation ideas:
% populations vs time, predator pop. vs. prey pop. quiver and trajectory
% use same initial conditions for every plot
% 1. explain the equations
% 2. show plots of the solution without capacity factor
% 3. show plots of the solution with capacity factor
% 4. talk about the special case where y == x (c * (r + alpha) == a * gamma)

clear;close all;
%solving numerically with ode45
tstart=0;tfinal=30;
global a; a=1.1; %growth rate of prey population
global alpha; alpha=0.4; %predation rate on prey population
global c; c=0.4; %death rate of predator population
global gamma; gamma=0.1; %growth rate of predator population
global r; r=0; %modifies the prey equation -> logistic growth expression
u0=zeros(2,1); %create u0 columns
u0(1)=1; u0(2)=1; %set up initial conditions for x(0) and y(0);



options=odeset('RelTol',1e-8);

[t,u]=ode45(@rhsPPS,[tstart,tfinal],u0,options); %solve for values
x1=u(:,1);y1=u(:,2); %get values for x and y

%now plot x and y as functions of t
plot(t,x1)
xlabel('time')
ylabel('population')
title('Predator-Prey Populations Over Time')
hold on
plot(t,y1)
legend('prey population','predator population')

%% QUIVER PLOTS
%% special equilibrium a=12,alpha=1,c=6,gamma=2
figure
x=0:.6:10;y=x;
[X,Y]=meshgrid(x,y);
dxdt=(a*X-r*X.^2)-alpha*X.*Y; %prey equation
dydt=(-c*Y)+gamma*X.*Y; %predator equation
quiver(X,Y,dxdt,dydt)
xx=0.1:.05:0.2; yy=9*ones(size(xx));
streamline(X,Y,dxdt,dydt,xx,yy) %do several initial conditions
hold on
yy1=1:1.5:5; xx1=9.*ones(size(yy1));
streamline(X,Y,dxdt,dydt,xx1,yy1)
streamline(X,Y,dxdt,dydt,0.8,9)
hold off
axis equal
xlabel('prey population')
ylabel('predator population')
title('Phase-Space for special equilibrium')
%% simple baboon vs. cheetah: a=1.1,alpha=0.4,c=0.4,gamma=0.1,r=0
figure
x=0:.6:20;y=0:0.3:10;
[X,Y]=meshgrid(x,y);
dxdt=(a*X-r*X.^2)-alpha*X.*Y; %prey equation
dydt=(-c*Y)+gamma*X.*Y; %predator equation
quiver(X,Y,dxdt,dydt)
yy=5:0.5:7; xx=9*ones(size(yy));
streamline(X,Y,dxdt,dydt,xx,yy,[1e-3,1e5]) %do several initial conditions
axis equal
xlabel('prey population')
ylabel('predator population')
title('Baboon population vs. Cheetah population phase space')
%% phase plots
figure
dxdt=(a.*x1-r.*(x1.^2))-alpha.*x1.*y1;
plot(x1,dxdt)
title('Predator population growth vs. predator population')
xlabel('Predator population')
ylabel('Predator population growth')

figure
dydt=-c.*y1+gamma.*x1.*y1;
plot(y1,dydt)
title('Prey population growth vs. prey population')
xlabel('Prey population')
ylabel('Prey population growth')

figure
plot(x1,y1)
title('Predator population vs. prey population')
xlabel('Prey population')
ylabel('Predator population')

figure
plot(dxdt,dydt)
title('Predator population growth vs. prey population growth')
xlabel('Prey population growth')
ylabel('Predator population growth')

%% equilibrium conditions
u0(1)=a/r; u0(2)=0; %set up initial conditions for x(0) and y(0);
[t,u]=ode45(@rhsPPS,[tstart,tfinal],u0,options); %solve for values
x1=u(:,1);y1=u(:,2); %get values for x and y

%now plot x and y as functions of t
plot(t,x1)
xlabel('time')
ylabel('population')
title('Predator-Prey Populations Over Time')
hold on
plot(t,y1)
legend('prey population','predator population')


figure
u0(1)=c/gamma; u0(2)=(a*gamma-c*r)/(alpha*gamma);
[t,u]=ode45(@rhsPPS,[tstart,tfinal],u0,options); %solve for values
x1=u(:,1);y1=u(:,2); %get values for x and y

%now plot x and y as functions of t
plot(t,x1)
xlabel('time')
ylabel('population')
ylim([0 5])
title('Predator-Prey Populations Over Time')
hold on
plot(t,y1)
legend('prey population','predator population')







