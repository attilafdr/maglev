clear all; clc; close all;

% Adatok

m = 0.032;        % Golyotomeg Nagyo: 107g; Kicsi: 32g
l = 0.07;       % Tekercshossz
a = 0.03^2;     % Felulet
n = 1200;       % Menetszam
mu= 4*pi*10^-5;% permeabilitas; mu_r=100
R = 30;       % Tekercs ellenallas
Q = 0.001917;      % Légrés_induktivitás, nagyobb golyó: 0.00458 golyó: ; kisebb: 0.001917
g = 0.980665;
Ts = 0.001;     % Sample time

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
Wd=c2d(W, Ts, 'zoh')
[Phi, Gamma, C, D] = ssdata(Wd)
P=eig(Wd)

% Irányíthatóság ellenõrzése

Mc=ctrb(Phi,Gamma);
if rank(Mc)>=length(P) 
    display('A szakasz iranyithato, az iranyithatosagi matrix rangja: '), display(rank(Mc))
else
    error('A szakasz nem iranyithato!')
end

szor=1;
%plc=[-2.9 -2.9 -2.9];
plc=exp(Ts.*[-10 -10 -10]);
%obs=[-320 -320 -320 -320];
obs=exp(Ts.*[-1800 -1800 -1800 -1800]);
K = acker(Phi,Gamma,plc)

eig(Phi-Gamma*K)
Wcl = ss(Phi-Gamma*K,Gamma,C,D,Ts);
t = 0:Ts:2;
u = 0*t;
x0 = [0.01 0 0];
lsim(Wcl,u,t,x0);
%bode(ss(A-B*K,B,C,D))

%initial(ss(A-B*K,[],C,[]),[0.01 0 0])

% Alapjelkövetés:

% A*Nx + B*Nu = 0
% C*Nx = 1

N=inv([Phi Gamma; C D])*[zeros(size(Gamma));1];
Nx=N(1:length(Phi))
Nu=N(length(Phi)+1)

% Allapotbecslo

% dx = Fx+Gy+Hu
% dxh = Fxh = (A - GC)xh

% Megfigyelhetõség ellenõrzése

Mo=ctrb(Phi,C');
if rank(Mo)>=length(P) 
    display('A szakasz megfigyelheto, a megfigyelhetosegi matrix rangja: '), display(rank(Mc))
else
    error('A szakasz nem iranyithato!')
end

% G=acker(Phi',C',[-50 -50 -50])';
% H=Gamma;
% F=Phi-G*C;


% Terhelesbecslovel:

% dx = Ax + B(u+d)

G=acker([Phi Gamma; 0 0 0 0]',[C 0]',obs)'
H=[Gamma;0];
%H=[[B;0]-G*C*B]
F=[Phi Gamma;0 0 0 0]-G*[C 0];

eig(F)


