function dxdt = HHIM_all(t,x,p) % This function contains nonlinear differential equations used in the function file to calculate Lyapunov Exponents
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
d1x1=x(8);d1x2=x(9);d1x3=x(10);d1x4=x(11);d1x5=x(12);d1x6=x(13);d1x7=x(14);
d2x1=x(15);d2x2=x(16);d2x3=x(17);d2x4=x(18);d2x5=x(19);d2x6=x(20);d2x7=x(21);
d3x1=x(22);d3x2=x(23);d3x3=x(24);d3x4=x(25);d3x5=x(26);d3x6=x(27);d3x7=x(28);
d4x1=x(29);d4x2=x(30);d4x3=x(31);d4x4=x(32);d4x5=x(33);d4x6=x(34);d4x7=x(35);
d5x1=x(36);d5x2=x(37);d5x3=x(38);d5x4=x(39);d5x5=x(40);d5x6=x(41);d5x7=x(42);
d6x1=x(43);d6x2=x(44);d6x3=x(45);d6x4=x(46);d6x5=x(47);d6x6=x(48);d6x7=x(49);
d7x1=x(50);d7x2=x(51);d7x3=x(52);d7x4=x(53);d7x5=x(54);d7x6=x(55);d7x7=x(56);


dxdt=zeros(56,1);
dxdt(1)= x2;
dxdt(2)= -2*z1*x2 - x1 - 2*z2*(x2-x6) - kr1*(x1-x5) - kr2*(x1-x3) - 2*z3*(x2-x4) - krnl1*x1^3 -krnl2*(x1-x3)^3 + F*(omega^2)*sin(omega*x(7));
dxdt(3)= x4;
dxdt(4)= -2*z3*ap*(x4-x2) - kr2*ap*(x3-x1) - krnl2*ap*(x3-x1)^3 ; 
dxdt(5)= x6;
dxdt(6)= -kr3*ap2*x5 - 2*z4*ap2*x6 - kr1*ap2*(x5-x1) - 2*z2*ap2*(x6-x2);
dxdt(7)=1;

dxdt(8)=d1x2;
dxdt(9)=d1x5*kr1 + 2*d1x4*z3 + 2*d1x6*z2 + d1x3*(kr2 + 3*krnl2*(x1 - x3)^2) - d1x1*(kr1 + kr2 + 3*krnl1*x1^2 + 3*krnl2*(x1 - x3)^2 + 1) - d1x2*(2*z1 + 2*z2 + 2*z3) + F*d1x7*omega^3*cos(omega*x7);
dxdt(10)=d1x4;
dxdt(11)=d1x1*(ap*kr2 + 3*ap*krnl2*(x1 - x3)^2) - d1x3*(ap*kr2 + 3*ap*krnl2*(x1 - x3)^2) + 2*ap*d1x2*z3 - 2*ap*d1x4*z3;
dxdt(12)=d1x6;
dxdt(13)=ap2*d1x1*kr1 - d1x6*(2*ap2*z2 + 2*ap2*z4) - d1x5*(ap2*kr1 + ap2*kr3) + 2*ap2*d1x2*z2;
dxdt(14)=0;


dxdt(15)=d2x2;
dxdt(16)=d2x5*kr1 + 2*d2x4*z3 + 2*d2x6*z2 + d2x3*(kr2 + 3*krnl2*(x1 - x3)^2) - d2x1*(kr1 + kr2 + 3*krnl1*x1^2 + 3*krnl2*(x1 - x3)^2 + 1) - d2x2*(2*z1 + 2*z2 + 2*z3) + F*d2x7*omega^3*cos(omega*x7);
dxdt(17)=d2x4;
dxdt(18)=d2x1*(ap*kr2 + 3*ap*krnl2*(x1 - x3)^2) - d2x3*(ap*kr2 + 3*ap*krnl2*(x1 - x3)^2) + 2*ap*d2x2*z3 - 2*ap*d2x4*z3;
dxdt(19)=d2x6;
dxdt(20)=ap2*d2x1*kr1 - d2x6*(2*ap2*z2 + 2*ap2*z4) - d2x5*(ap2*kr1 + ap2*kr3) + 2*ap2*d2x2*z2;
dxdt(21)=0;

dxdt(22)=d3x2;
dxdt(23)=d3x5*kr1 + 2*d3x4*z3 + 2*d3x6*z2 + d3x3*(kr2 + 3*krnl2*(x1 - x3)^2) - d3x1*(kr1 + kr2 + 3*krnl1*x1^2 + 3*krnl2*(x1 - x3)^2 + 1) - d3x2*(2*z1 + 2*z2 + 2*z3) + F*d3x7*omega^3*cos(omega*x7);
dxdt(24)=d3x4;
dxdt(25)=d3x1*(ap*kr2 + 3*ap*krnl2*(x1 - x3)^2) - d3x3*(ap*kr2 + 3*ap*krnl2*(x1 - x3)^2) + 2*ap*d3x2*z3 - 2*ap*d3x4*z3;
dxdt(26)=d3x6;
dxdt(27)=ap2*d3x1*kr1 - d3x6*(2*ap2*z2 + 2*ap2*z4) - d3x5*(ap2*kr1 + ap2*kr3) + 2*ap2*d3x2*z2;
dxdt(28)=0;

dxdt(29)=d4x2;
dxdt(30)=d4x5*kr1 + 2*d4x4*z3 + 2*d4x6*z2 + d4x3*(kr2 + 3*krnl2*(x1 - x3)^2) - d4x1*(kr1 + kr2 + 3*krnl1*x1^2 + 3*krnl2*(x1 - x3)^2 + 1) - d4x2*(2*z1 + 2*z2 + 2*z3) + F*d4x7*omega^3*cos(omega*x7);
dxdt(31)=d4x4;
dxdt(32)=d4x1*(ap*kr2 + 3*ap*krnl2*(x1 - x3)^2) - d4x3*(ap*kr2 + 3*ap*krnl2*(x1 - x3)^2) + 2*ap*d4x2*z3 - 2*ap*d4x4*z3;
dxdt(33)=d4x6;
dxdt(34)=ap2*d4x1*kr1 - d4x6*(2*ap2*z2 + 2*ap2*z4) - d4x5*(ap2*kr1 + ap2*kr3) + 2*ap2*d4x2*z2;
dxdt(35)=0;

dxdt(36)=d5x2;
dxdt(37)=d5x5*kr1 + 2*d5x4*z3 + 2*d5x6*z2 + d5x3*(kr2 + 3*krnl2*(x1 - x3)^2) - d5x1*(kr1 + kr2 + 3*krnl1*x1^2 + 3*krnl2*(x1 - x3)^2 + 1) - d5x2*(2*z1 + 2*z2 + 2*z3) + F*d5x7*omega^3*cos(omega*x7);
dxdt(38)=d5x4;
dxdt(39)=d5x1*(ap*kr2 + 3*ap*krnl2*(x1 - x3)^2) - d5x3*(ap*kr2 + 3*ap*krnl2*(x1 - x3)^2) + 2*ap*d5x2*z3 - 2*ap*d5x4*z3;
dxdt(40)=d5x6;
dxdt(41)=ap2*d5x1*kr1 - d5x6*(2*ap2*z2 + 2*ap2*z4) - d5x5*(ap2*kr1 + ap2*kr3) + 2*ap2*d5x2*z2;
dxdt(42)=0;

dxdt(43)=d6x2;
dxdt(44)=d6x5*kr1 + 2*d6x4*z3 + 2*d6x6*z2 + d6x3*(kr2 + 3*krnl2*(x1 - x3)^2) - d6x1*(kr1 + kr2 + 3*krnl1*x1^2 + 3*krnl2*(x1 - x3)^2 + 1) - d6x2*(2*z1 + 2*z2 + 2*z3) + F*d6x7*omega^3*cos(omega*x7);
dxdt(45)=d6x4;
dxdt(46)=d6x1*(ap*kr2 + 3*ap*krnl2*(x1 - x3)^2) - d6x3*(ap*kr2 + 3*ap*krnl2*(x1 - x3)^2) + 2*ap*d6x2*z3 - 2*ap*d6x4*z3;
dxdt(47)=d6x6;
dxdt(48)=ap2*d6x1*kr1 - d6x6*(2*ap2*z2 + 2*ap2*z4) - d6x5*(ap2*kr1 + ap2*kr3) + 2*ap2*d6x2*z2;
dxdt(49)=0;

dxdt(50)=d7x2;
dxdt(51)=d7x5*kr1 + 2*d7x4*z3 + 2*d7x6*z2 + d7x3*(kr2 + 3*krnl2*(x1 - x3)^2) - d7x1*(kr1 + kr2 + 3*krnl1*x1^2 + 3*krnl2*(x1 - x3)^2 + 1) - d7x2*(2*z1 + 2*z2 + 2*z3) + F*d7x7*omega^3*cos(omega*x7);
dxdt(52)=d7x4;
dxdt(53)=d7x1*(ap*kr2 + 3*ap*krnl2*(x1 - x3)^2) - d7x3*(ap*kr2 + 3*ap*krnl2*(x1 - x3)^2) + 2*ap*d7x2*z3 - 2*ap*d7x4*z3;
dxdt(54)=d7x6;
dxdt(55)=ap2*d7x1*kr1 - d7x6*(2*ap2*z2 + 2*ap2*z4) - d7x5*(ap2*kr1 + ap2*kr3) + 2*ap2*d7x2*z2;
dxdt(56)=0;

dxdt = dxdt(:);