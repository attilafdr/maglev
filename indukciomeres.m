% Inductance measurement
% Dspace: -10..10V --dspace--> -1..1 --*100--> -100..100
clc, clear all;
Rr = 220.3; % Ohm
R = Rr + 17.6; % Ohm
Phi0 = 0.383 % rad
A = 4 * Rr / R; % V
du = 0.015; % V
Phi1 = (du / A) + Phi0 % rad
w = 10*2*pi; % rad/sec
X0 = R * tan(Phi0);
X1 = R * tan(Phi1);
L0 = X0 / w
L1 = X1 / w
dL = L1-L0

% measurement no.2
Rr = 6; % Ohm
R = Rr + 17.6; % Ohm
Phi0 = 0.12305 % rad
A = 0.6 * Rr / R; % V
du = 0.006; % V
Phi1 = (du / A) + Phi0 % rad
w = 10*2*pi; % rad/sec
X0 = R * tan(Phi0);
X1 = R * tan(Phi1);
L0 = X0 / w
L1 = X1 / w
dL = L1-L0