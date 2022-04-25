% These script calculates Lyapunov exponents in a nonlinear system using
% the method proposed by Wolf et al. 
% Wolf, Alan, et al. "Determining Lyapunov exponents from a time series." Physica D: nonlinear phenomena 16.3 (1985): 285-317.

%% Obtaining Lyapunov exponents for Nonlinear "HIM-HAS-ABS" system using ODE45
clear
clc
close all
%% Defining nonlinear system parameters and parameters for numerical solver
w=700;
[z1,z2,kr1,kr2,z3,krnl1,krnl2,ap,kr3,ap2,z4,omega]=nondim_p(w);
% Defining parameters for numerical solver
% keyboard
xo= [0 0 0 0 0 0 0]';
t_step=0.01;                  
opts = odeset('RelTol',1e-7,'AbsTol',1e-7); 


%% Obtaining Lyapunov exponents for multiple values of system parameter "Omega" and "F"
% Defining system parameters and parameters for numerical solver
%%
b=fliplr(0:10:500);         % This parameter should be the same as used in "Initial_conditions_for_LEcode.m"
kr2=0.07;

%%
om=linspace(0.8,2,60);      % This parameter should be the same as used in "Initial_conditions_for_LEcode.m"
f=length(om);
s=length(b);

tic
for mm=1:f
filename=sprintf('IC_0p07kr2_%d',mm);  
load(filename)
parfor j=1:s 
% Defining system parameters corresponding to initial conditions
F=b(j);
omega=om(mm);

% Defining starting point for numerical simulation based on initial conditions obtained from "Initial_conditions_for_LEcode.m"
% "xx" is obtained when 'filename' is loaded 
x1=xx(j,1);x2=xx(j,2);x3=xx(j,3);x4=xx(j,4);x5=xx(j,5);x6=xx(j,6);x7=xx(j,7);
t_final= 8*pi*500/omega;
pp=[z1 z2 kr1 kr2 z3 krnl1 krnl2 F ap kr3 ap2 z4 omega].'; 

% Defining initial re-orthonormalization time 
t1=0;         
% Calling function to obtain Lyapunov Exponents's. Method to compute
% Lyapunov equations is stored in the function file LE_mex
% "LE_mex" contains same content as "LE.m"
[la1,la2,la3,la4,la5,la6,la7]=LE_mex(t1,t_step,t_final,pp,x1,x2,x3,x4,x5,x6,x7);
LE(j,:)=[la1,la2,la3,la4,la5,la6,la7,F];

end    

% Saving Lyapunov Exponent values calculated for one set of system parameters  
filename=sprintf('LEs_0p07kr2_%d',mm);  
save(filename,'LE');
mm
end
toc

