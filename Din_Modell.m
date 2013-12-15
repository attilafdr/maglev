clear all; clc; close all;

% Adatok

m = 0.032;        % Golyotomeg Nagyo: 107g; Kicsi: 32g
l = 0.07;       % Tekercshossz
a = 0.03^2;     % Felulet
n = 1200;       % Menetszam
mu= 4*pi*10^-5;% permeabilitas; mu_r=100
R = 17.6;       % Tekercs ellenallas
Q = 0.001917;      % Légrés_induktivitás, nagyobb golyó: 0.00458 golyó: ; kisebb: 0.001917
g = 9.80665;

% Szamitott adatok

L0 = mu*a*n^2/l;    % Tekercs induktivitas
%L0 = 0.5; %probaertek

% L(x)=L0 + Q/y
% L0: konstans indukcio a golyo nelkul
% Q: legreshez tartozo induktivitas

% u(t) = R*i+di/di*(L*i)
% u(t) = R*i+(-Q/y^2*dy/dt*i+(Q/y+L0)*di/dt

% di/dt=(u-R*i+Q*v*i/y^2)/(Q/y+L0)

% Allapotvaltozok

% x1 = x;
% x2 = v;
% x3 = I;

% Rendszeregyenletek:

% v = dx/dt
% m*a=m*g-Q*i^2/(2*y^2)
% u(t)=R*i+(-Q/y^2*dy/dt*i+(Q/y+L0)*di/dt
echo on
% Dinamikus modell:

% --  --   --                                  -- --         --
% | x1'|   | x2                                 | |     0     |
% | x2'| = | g-Q*x3^2/(2*m*x1^2)                |+|     0     |*u
% | x3'|   | (u-R*x3+Q*x2*x3/x1^2)/(Q/x1+L0)    | | x1/Q+L0*x1|
% --  --   --                                  -- --         --

% y = x1
echo off

% x0: Munkapont: egyensulyi pont
% x0 = [x01 x02 x03]

x01 = 0.005;
x02 = 0;
x03 = x01*sqrt(2*m*g/Q); 

u0  = R*x03;

A31=((u0-R*x03)*(Q+L0*x01)-L0*x01*(u0-R*x03))/(Q+L0*x01^2)+Q*x01^2*x02*x03*(Q+L0*x01)-(Q^2*x01^2*x02*x03*(2+3*L0*x01))/(Q*x01^2+L0*x01^2)

echo on
A = [0                  1                               0;
    Q*x03^2/(m*x01^3)   0                               -Q*x03/(m*x01^2);
    0                   Q*x01*x03/x01^2/(Q+L0*x01)      (Q*x02-R*x01^2)/x01/(Q+L0*x01)];

B = [0;
     0;
     x01/(Q+L0*x01)];
 
C = [1 0 0];

D=0;
echo off

W=ss(A,B,C,D);
P=eig(A)

% Irányíthatóság ellenõrzése

Mc=ctrb(A,B);
if rank(Mc)>=length(P) 
    display('A szakasz iranyithato, az iranyithatosagi matrix rangja: '), display(rank(Mc))
else
    error('A szakasz nem iranyithato!')
end

szor=1;
plc=[-2.9 -2.9 -2.9]; %-0.5 -18 -1.5
obs=[-320 -320 -320 -320];
K = acker(A,B,plc)

eig(A-B*K)

t = 0:0.01:2;
u = 0*t;
x0 = [0.01 0 0];
%lsim(A-B*K,B,C,D,u,t,x0);
%bode(ss(A-B*K,B,C,D))

%initial(ss(A-B*K,[],C,[]),[0.01 0 0])

% Alapjelkövetés:

% A*Nx + B*Nu = 0
% C*Nx = 1

N=inv([A B; C D])*[zeros(size(B));1];
Nx=N(1:length(A))
Nu=N(length(A)+1)

% Allapotbecslo

% dx = Fx+Gy+Hu
% dxh = Fxh = (A - GC)xh

% Megfigyelhetõség ellenõrzése

Mo=ctrb(A,C');
if rank(Mo)>=length(P) 
    display('A szakasz megfigyelheto, a megfigyelhetosegi matrix rangja: '), display(rank(Mc))
else
    error('A szakasz nem iranyithato!')
end

G=acker(A',C',[-50 -50 -50])';
H=B;
F=A-G*C;

% Terhelesbecslovel:

% dx = Ax + B(u+d)

G=acker([A B; 0 0 0 0]',[C 0]',obs)'
H=[B;0];
%H=[[B;0]-G*C*B]
F=[A B;0 0 0 0]-G*[C 0];




Wc=ss(A-B*K,B,C,D);
[Wn,Z,P]=damp(ss(A-B*K,B,C,D))

