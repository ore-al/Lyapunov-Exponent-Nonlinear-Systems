clear
clc
close all

%% Obtaining last vector in time series to use as initial condition for code to obtain Lyapunov exponents
% Obtaining system parameters
w=700;
[z1,z2,kr1,kr2,z3,krnl1,krnl2,ap,kr3,ap2,z4,omega]=nondim_p(w);
kr2=0.07;

%%  Defining forcing values to be used for Lyapunov Exponent chart in parameter space F-Omega
b=fliplr(0:10:500); % 
%%
s=length(b);
xx=zeros(s,7);

% Defining frequency vector 
om=linspace(0.8,2,60);
f=length(om);
% Defining parameters for numerical solver
xo= [0 0 0 0 0 0 0]';
opts = odeset('RelTol',1e-7,'AbsTol',1e-7); 

%% Finding and storing steady-state displacement for different variations of nonlinear system 
for mm=1:f
for j=1:s
    
F=b(j);
omega=om(mm);
tspan= [0 8*pi*500/omega];
for kk=1:3
pp=[z1 z2 kr1 kr2 z3 krnl1 krnl2 F ap kr3 ap2 z4 omega].'; 
[t,x]=ode45(@(t,x) Nondimensional_ODE(t,x,pp),tspan,xo,opts);
xo=x(end,:)'; 
end 

x1=x(end,1);
x2=x(end,2);
x3=x(end,3);
x4=x(end,4);
x5=x(end,5);
x6=x(end,6);
x7=x(end,7);

xx(j,:)=[x1 x2 x3 x4 x5 x6 x7];
end
filename=sprintf('IC_0p07kr2_%d',mm);            
save(filename,'xx','omega','F');
% keyboard
end

