function [la1,la2,la3,la4,la5,la6,la7]=LE(t1,t_step,t_final,pp,x1,x2,x3,x4,x5,x6,x7) %#codegen
lam1=0;
lam2=0;
lam3=0;
lam4=0;
lam5=0;
lam6=0;
lam7=0;

a=length(t1:t_step:t_final);
time=zeros(a,1);
xx1=zeros(a,1);
xx2=zeros(a,1);
xx3=zeros(a,1);
xx4=zeros(a,1);
xx5=zeros(a,1);
xx6=zeros(a,1);
xx7=zeros(a,1);

d1x1=0;d1x2=0;d1x3=0;d1x4=0;d1x5=0;d1x6=0;d1x7=1;
d2x1=0;d2x2=0;d2x3=0;d2x4=0;d2x5=0;d2x6=1;d2x7=0;
d3x1=0;d3x2=0;d3x3=0;d3x4=0;d3x5=1;d3x6=0;d3x7=0;
d4x1=0;d4x2=0;d4x3=0;d4x4=1;d4x5=0;d4x6=0;d4x7=0;
d5x1=0;d5x2=0;d5x3=1;d5x4=0;d5x5=0;d5x6=0;d5x7=0;
d6x1=0;d6x2=1;d6x3=0;d6x4=0;d6x5=0;d6x6=0;d6x7=0;
d7x1=1;d7x2=0;d7x3=0;d7x4=0;d7x5=0;d7x6=0;d7x7=0;

lambda_1=zeros(1,a);lambda1=zeros(1,a);
lambda_2=zeros(1,a);lambda2=zeros(1,a);
lambda_3=zeros(1,a);lambda3=zeros(1,a);
lambda_4=zeros(1,a);lambda4=zeros(1,a);
lambda_5=zeros(1,a);lambda5=zeros(1,a);
lambda_6=zeros(1,a);lambda6=zeros(1,a);
lambda_7=zeros(1,a);lambda7=zeros(1,a);



%% Evaluating Jacobian for system (Euqations in function HHIM_i were obtained using this function i.e. i =1,2,...)
% [J]=Lyap_Jac(p);

tic
for i=1:a
    %% Storing time and position values from simulation
    time(i)=t1;
    xx1(i)=x1;
    xx2(i)=x2;
    xx3(i)=x3;
    xx4(i)=x4;
    xx5(i)=x5;
    xx6(i)=x6;
    xx7(i)=x7;
    
    %% defining initial conditions and running solver for re-orthonomormalization step
    tspan=t1+[0,t_step];           % increasing size of tspan for each re-orthonormalization step  
    ao=[d1x1,d1x2,d1x3,d1x4,d1x5,d1x6,d1x7];
    bo=[d2x1,d2x2,d2x3,d2x4,d2x5,d2x6,d2x7];
    co=[d3x1,d3x2,d3x3,d3x4,d3x5,d3x6,d3x7];
    do=[d4x1,d4x2,d4x3,d4x4,d4x5,d4x6,d4x7];
    eo=[d5x1,d5x2,d5x3,d5x4,d5x5,d5x6,d5x7];
    fo=[d6x1,d6x2,d6x3,d6x4,d6x5,d6x6,d6x7];
    go=[d7x1,d7x2,d7x3,d7x4,d7x5,d7x6,d7x7];
% end
    xo=[x1 x2 x3 x4 x5 x6 x7 ao bo co do eo fo go];
    [A]=ode5(@(t,x) HHIM_all(t,x,pp),tspan,xo);
%     keyboard
    % If you experience "too many output argument" error that probably has to
    % do with how you are calling function 
%     [t,A]=ode45(@(t,A) HHIM_all(t,A,pp),tspan,xo);
        
    %% setting initial conditions for time t+t_step based on end-value of 'A' from numerical analysis
    x1=A(size(A,1),1); 
    x2=A(size(A,1),2);
    x3=A(size(A,1),3); 
    x4=A(size(A,1),4);
    x5=A(size(A,1),5); 
    x6=A(size(A,1),6);
    x7=A(size(A,1),7);
    
    d1x1=A(size(A,1),8); 
    d1x2=A(size(A,1),9);
    d1x3=A(size(A,1),10); 
    d1x4=A(size(A,1),11);
    d1x5=A(size(A,1),12); 
    d1x6=A(size(A,1),13);
    d1x7=A(size(A,1),14);
    
    d2x1=A(size(A,1),15); 
    d2x2=A(size(A,1),16);
    d2x3=A(size(A,1),17); 
    d2x4=A(size(A,1),18);
    d2x5=A(size(A,1),19); 
    d2x6=A(size(A,1),20);
    d2x7=A(size(A,1),21);
    
    d3x1=A(size(A,1),22); 
    d3x2=A(size(A,1),23);
    d3x3=A(size(A,1),24); 
    d3x4=A(size(A,1),25);
    d3x5=A(size(A,1),26); 
    d3x6=A(size(A,1),27);
    d3x7=A(size(A,1),28);
    
    d4x1=A(size(A,1),29); 
    d4x2=A(size(A,1),30);
    d4x3=A(size(A,1),31); 
    d4x4=A(size(A,1),32);
    d4x5=A(size(A,1),33); 
    d4x6=A(size(A,1),34);
    d4x7=A(size(A,1),35);
    
    d5x1=A(size(A,1),36); 
    d5x2=A(size(A,1),37);
    d5x3=A(size(A,1),38); 
    d5x4=A(size(A,1),39);
    d5x5=A(size(A,1),40); 
    d5x6=A(size(A,1),41);
    d5x7=A(size(A,1),42);
       
    d6x1=A(size(A,1),43); 
    d6x2=A(size(A,1),44);
    d6x3=A(size(A,1),45); 
    d6x4=A(size(A,1),46);
    d6x5=A(size(A,1),47); 
    d6x6=A(size(A,1),48);
    d6x7=A(size(A,1),49);
    
    d7x1=A(size(A,1),50); 
    d7x2=A(size(A,1),51);
    d7x3=A(size(A,1),52); 
    d7x4=A(size(A,1),53);
    d7x5=A(size(A,1),54); 
    d7x6=A(size(A,1),55);
    d7x7=A(size(A,1),56);
     
       
    %% Calculating Magnitude of dx,dy,dz at t'', and finite and Instantaneous lambda1 for 1st-DOF
    m1_t2=sqrt((d1x1^2)+(d1x2^2)+(d1x3^2)+(d1x4^2)+(d1x5^2)+(d1x6^2)+(d1x7^2));
    lambda_1(i)=(1/t_step)*log(m1_t2);
    lam1=lam1+lambda_1(i);             %summing instantaneous lyapunov exponent for each iteration
    lambda1(i)=lam1/length(lambda_1);  %computing Finite time lyapunov exponent for each iteration
   
    %% Normalizing dx,dy and dz for 1st-DOF 
    V1=[d1x1 d1x2 d1x3 d1x4 d1x5 d1x6 d1x7];
    V1= V1/norm(V1);
    d1x1=V1(1); d1x2=V1(2); d1x3=V1(3);d1x4=V1(4); d1x5=V1(5); d1x6=V1(6); d1x7=V1(7);
    
    %% Calculating Magnitude of dx,dy,dz at t'', and finite and Instantaneous lambda1 for 2nd-DOF
    V2=[d2x1 d2x2 d2x3 d2x4 d2x5 d2x6 d2x7];
    V2= V2-dot(V2,V1)*V1;
    m2_t2=sqrt((V2(1)^2)+(V2(2)^2)+(V2(3)^2)+(V2(4)^2)+(V2(5)^2)+(V2(6)^2)+(V2(7)^2));
    lambda_2(i)=(1/t_step)*log(m2_t2);
    lam2=lam2+lambda_2(i);             %summing instantaneous lyapunov exponent for each iteration
    lambda2(i)=lam2/length(lambda_2);  %computing Finite time lyapunov exponent for each iteration
    
    %% Normalizing dx,dy and dz for 2nd-DOF 
    V2= V2/norm(V2);
    d2x1=V2(1); d2x2=V2(2); d2x3=V2(3);d2x4=V2(4); d2x5=V2(5); d2x6=V2(6);d2x7=V2(7);
    %% Calculating Magnitude of dx,dy,dz at t'', and finite and Instantaneous lambda1 for 3rd-DOF
    V3=[d3x1,d3x2,d3x3,d3x4,d3x5,d3x6,d3x7];
    V3=V3-dot(V3,V2)*V2-dot(V3,V1)*V1;
    m3_t2=sqrt((V3(1)^2)+(V3(2)^2)+(V3(3)^2)+(V3(4)^2)+(V3(5)^2)+(V3(6)^2)+(V3(7)^2));
    lambda_3(i)=(1/t_step)*log(m3_t2);
    lam3=lam3+lambda_3(i);             %summing instantaneous lyapunov exponent for each iteration
    lambda3(i)=lam3/length(lambda_3);  %computing Finite time lyapunov exponent for each iteration
    
    %% Normalizing dx,dy and dz for 3rd-DOF 
    V3= V3/norm(V3);
    d3x1=V3(1); d3x2=V3(2); d3x3=V3(3);d3x4=V3(4); d3x5=V3(5); d3x6=V3(6);d3x7=V3(7);
    %% Calculating Magnitude of dx,dy,dz at t'', and finite and Instantaneous lambda1 for 4th-DOF
    V4=[d4x1,d4x2,d4x3,d4x4,d4x5,d4x6,d4x7];
    V4=V4-dot(V4,V3)*V3-dot(V4,V2)*V2-dot(V4,V1)*V1;
    m4_t2=sqrt((V4(1)^2)+(V4(2)^2)+(V4(3)^2)+(V4(4)^2)+(V4(5)^2)+(V4(6)^2)+(V4(7)^2));
    lambda_4(i)=(1/t_step)*log(m4_t2);
    lam4=lam4+lambda_4(i);             %summing instantaneous lyapunov exponent for each iteration
    lambda4(i)=lam4/length(lambda_4);  %computing Finite time lyapunov exponent for each iteration
    
    %% Normalizing dx,dy and dz for 4th-DOF 
    V4= V4/norm(V4);
    d4x1=V4(1); d4x2=V4(2); d4x3=V4(3);d4x4=V4(4); d4x5=V4(5); d4x6=V4(6); d4x7=V4(7);
    
    %% Calculating Magnitude of dx,dy,dz at t'', and finite and Instantaneous lambda1 for 5th-DOF
    V5=[d5x1,d5x2,d5x3,d5x4,d5x5,d5x6,d5x7];
    V5=V5-dot(V5,V4)*V4-dot(V5,V3)*V3-dot(V5,V2)*V2-dot(V5,V1)*V1;
    m5_t2=sqrt((V5(1)^2)+(V5(2)^2)+(V5(3)^2)+(V5(4)^2)+(V5(5)^2)+(V5(6)^2)+(V5(7)^2));
    lambda_5(i)=(1/t_step)*log(m5_t2);
    lam5=lam5+lambda_5(i);             %summing instantaneous lyapunov exponent for each iteration
    lambda5(i)=lam5/length(lambda_5);  %computing Finite time lyapunov exponent for each iteration
    
    %% Normalizing dx,dy and dz for 5th-DOF 
    V5= V5/norm(V5);
    d5x1=V5(1); d5x2=V5(2); d5x3=V5(3);d5x4=V5(4); d5x5=V5(5); d5x6=V5(6);d5x7=V5(7);
    
    %% Calculating Magnitude of dx,dy,dz at t'', and finite and Instantaneous lambda1 for 6th-DOF
    V6=[d6x1,d6x2,d6x3,d6x4,d6x5,d6x6,d6x7];
    V6=V6-dot(V6,V5)*V5-dot(V6,V4)*V4-dot(V6,V3)*V3-dot(V6,V2)*V2-dot(V6,V1)*V1;
    m6_t2=sqrt((V6(1)^2)+(V6(2)^2)+(V6(3)^2)+(V6(4)^2)+(V6(5)^2)+(V6(6)^2)+(V6(7)^2));
    lambda_6(i)=(1/t_step)*log(m6_t2);
    lam6=lam6+lambda_6(i);             %summing instantaneous lyapunov exponent for each iteration
    lambda6(i)=lam6/length(lambda_6);  %computing Finite time lyapunov exponent for each iteration
    
    %% Normalizing dx,dy and dz for 6th-DOF 
    V6= V6/norm(V6);
    d6x1=V6(1); d6x2=V6(2); d6x3=V6(3);d6x4=V6(4); d6x5=V6(5); d6x6=V6(6); d6x7=V6(7);
    
     %% Calculating Magnitude of dx,dy,dz at t'', and finite and Instantaneous lambda1 for 7th-DOF
    V7=[d7x1,d7x2,d7x3,d7x4,d7x5,d7x6,d7x7];
    V7=V7-dot(V7,V6)*V6-dot(V7,V5)*V5-dot(V7,V4)*V4-dot(V7,V3)*V3-dot(V7,V2)*V2-dot(V7,V1)*V1;
    m7_t2=sqrt((V7(1)^2)+(V7(2)^2)+(V7(3)^2)+(V7(4)^2)+(V7(5)^2)+(V7(6)^2)+(V7(7)^2));
    lambda_7(i)=(1/t_step)*log(m7_t2);
    lam7=lam7+lambda_7(i);             %summing instantaneous lyapunov exponent for each iteration
    lambda7(i)=lam7/length(lambda_7);  %computing Finite time lyapunov exponent for each iteration
    
    %% Normalizing dx,dy and dz for 6th-DOF 
    V7= V7/norm(V7);
    d7x1=V7(1); d7x2=V7(2); d7x3=V7(3);d7x4=V7(4); d7x5=V7(5); d7x6=V7(6);d7x7=V7(7);
%     keyboard
    
    
t1=t1+t_step;
end
toc
la1=lambda1(end);
la2=lambda2(end);
la3=lambda3(end);
la4=lambda4(end);
la5=lambda5(end);
la6=lambda6(end);
la7=lambda7(end);
