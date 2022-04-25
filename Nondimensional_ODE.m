function dxdt = Nondimensional_ODE(t,x,p)  %#codegen
% This function file contains the nonlinear differential equations for the
% system of equations evaluated
z1=p(1);
z2=p(2);
kr1=p(3);
kr2=p(4);
z3=p(5);
krnl1=p(6);
krnl2=p(7);
F=p(8);
ap=p(9);
kr3=p(10);
ap2=p(11);
z4=p(12);
omega=p(13);


x1=x(1);x2=x(2);x3=x(3);x4=x(4);x5=x(5);x6=x(6);x7=x(7);
dxdt=zeros(7,1);
dxdt(1)= x2;
dxdt(2)= -2*z1*x2 - x1 - 2*z2*(x2-x6) - kr1*(x1-x5) - kr2*(x1-x3) - 2*z3*(x2-x4) - krnl1*x1^3 -krnl2*(x1-x3)^3 + F*(omega^2)*sin(omega*x(7));
dxdt(3)= x4;
dxdt(4)= -2*z3*ap*(x4-x2) - kr2*ap*(x3-x1) - krnl2*ap*(x3-x1)^3 ; 
dxdt(5)= x6;
dxdt(6)= -kr3*ap2*x5 - 2*z4*ap2*x6 - kr1*ap2*(x5-x1) - 2*z2*ap2*(x6-x2);
dxdt(7)= 1;
dxdt = dxdt(:);
end