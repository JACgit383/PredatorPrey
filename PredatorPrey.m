%% Predator-Prey System
%by Jacob Christensen and Lucas Barnes

%system of first order differentials, x denotes the population
%of the prey species, y denotes the population of predator species

%a,c,alpha,gamma,r -> all positive constants

% presentation ideas:
% populations vs time, predator pop. vs. prey pop. quiver and trajectory
% use same initial conditions for every plot

%% Motivation
% talk about why we think this system is interesting
% wanted to learn more about the relationship between predator and prey through
% this model

%% Description of the Physical System

% The system describes the ebb and flow of the predator/prey population
% based on the rates at which each species grows, dies, and the rate at
% which the predators prey on the prey. If there are a lot of prey to eat,
% the predators will flourish and the prey will die off, whereas if there
% are not a lot of prey, the predators will not have enough to eat and they
% will die off and the prey will flourish.

%% Differential equations that govern the system

%Lotka-Volterra equations
%x'=a*x-alpha*x*y;
%y'=-c*y+gamma*x*y;

%x represents the prey population
%y represents the predator population

%a = growth rate of prey population
%alpha = predation rate on prey population
%c = death rate of predator population
%gamma = growth rate of predator population



%% 2. show plots of the solution without capacity factor
% talk about how we solved the equations
%a=1.1,alpha=0.4,c=0.4,gamma=0.1,r=0
clear;close all;
myFigureDefaultsTBN()
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
xlabel('Time')
ylabel('Population')
title('Predator-Prey Populations Over Time')
hold on
plot(t,y1)
legend('prey population','predator population')

figure
x=0:1.5:20;y=0:0.75:10;
[X,Y]=meshgrid(x,y);
dxdt=(a*X-r*X.^2)-alpha*X.*Y; %prey equation
dydt=(-c*Y)+gamma*X.*Y; %predator equation
quiver(X,Y,dxdt,dydt)
xlim([0,20])
ylim([0,10])
yy=5:0.5:7; xx=9*ones(size(yy));
streamline(X,Y,dxdt,dydt,xx,yy,[1e-3,1e5]) %do several initial conditions
xlim([0,20])
ylim([0,10])
xlabel('Prey population')
ylabel('Predator population')
title('Baboon population vs. Cheetah population phase space')


%% 3. show plots of the solution with capacity factor
% mention that this system represents an ideal situation in which the
% environment can sustain unlimited populations of prey or predators

%x'=a*x-r*(x^2)-alpha*x*y;
%r accounts for env. factors
tstart=0;tfinal=60;
r=0.05; %modifies the prey equation -> logistic growth expression
[t,u]=ode45(@rhsPPS,[tstart,tfinal],u0,options); %solve for values
x1=u(:,1);y1=u(:,2); %get values for x and y

%now plot x and y as functions of t
plot(t,x1)
xlabel('Time')
ylabel('Population')
title('Predator-Prey Populations Over Time')
hold on
plot(t,y1)
legend('prey population','predator population')
figure
x=0:.6:12;y=0:0.3:10;
[X,Y]=meshgrid(x,y);
dxdt=(a*X-r*X.^2)-alpha*X.*Y; %prey equation
dydt=(-c*Y)+gamma*X.*Y; %predator equation
quiver(X,Y,dxdt,dydt)
yy=5:0.5:7; xx=9*ones(size(yy));
streamline(X,Y,dxdt,dydt,xx,yy,[1e-3,1e5]) %do several initial conditions
axis equal
xlabel('Prey population')
ylabel('Predator population')
title('Baboon population vs. Cheetah population phase space')
ylim([0,8])
xlim([0,12])

%% 4. talk about the special case where y == x == c * (r + alpha) == a * gamma
%a=12,alpha=1,c=6,gamma=2
a=12; %growth rate of prey population
alpha=1; %predation rate on prey population
c=6; %death rate of predator population
gamma=2; %growth rate of predator population
r=3; %modifies the prey equation -> logistic growth expression

[t,u]=ode45(@rhsPPS,[tstart,tfinal],u0,options); %solve for values
x1=u(:,1);y1=u(:,2); %get values for x and y

%now plot x and y as functions of t
plot(t,x1)
xlabel('Time')
ylabel('Population')
title('Predator-Prey Populations Over Time')
hold on
plot(t,y1)
xlim([0,10])
legend('prey population','predator population')


figure
x=0:0.8:10;y=x;
[X,Y]=meshgrid(x,y);
dxdt=(a*X-r*X.^2)-alpha*X.*Y; %prey equation
dydt=(-c*Y)+gamma*X.*Y; %predator equation
quiver(X,Y,dxdt,dydt)
ylim([0,10])
xx=0.1:.05:0.2; yy=9*ones(size(xx));
streamline(X,Y,dxdt,dydt,xx,yy) %do several initial conditions
hold on
yy1=1:1.5:5; xx1=9.*ones(size(yy1));
streamline(X,Y,dxdt,dydt,xx1,yy1)
streamline(X,Y,dxdt,dydt,0.8,9)
hold off
xlabel('Prey population')
ylabel('Predator population')
title('Phase-Space for Equal Population Equilibrium')

%% references

% "The Predator Prey Equations," PSU, http://www.personal.psu.edu/sxt104/class/Math251/Notes-Predator-Prey.pdf
% "Lotka Volterra Equations", Wikipedia, https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations

%%
clear;close all;
%solving numerically with ode45
tstart=0;tfinal=30;
global a; a=12; %growth rate of prey population
global alpha; alpha=1; %predation rate on prey population
global c; c=6; %death rate of predator population
global gamma; gamma=2; %growth rate of predator population
global r; r=3; %modifies the prey equation -> logistic growth expression
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
x=0:0.9:10;y=x;
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
