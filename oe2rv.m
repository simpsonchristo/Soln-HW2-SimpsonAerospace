function [X] = oe2rv(oe, f)
%TWOBODYEPHEMERIS Determine position and velocity at time t by solving
%Kepler's Equation
%Solve Kepler's Equation for eccentric anomaly and determine true anomaly
%to update the initial state to some state at time t. This propagation can
%only be used in the two-body case. E0 and f are solved using the orbital
%elements at the initial time, t0. 
%Earth orbit
mu = 3.9860044e+14; %m^3/s^2

a    = oe(1);%semimajor axis, meters
e    = oe(2);%eccentricity
i    = oe(3);%inclination, degrees
Si = sind(i); Ci = cosd(i);
raan = oe(4);%ascending node, degrees
Sraan = sind(raan); Craan = cosd(raan);
w    = oe(5);%aop, degrees
Sw = sind(w); Cw = cosd(w);
M    = oe(6);%mean anomaly, degrees
Sm = sind(M); Cm = cosd(M);

Q = [((Cw*Craan)-(Sw*Sraan*Ci)) ((-Sw*Craan)-(Cw*Sraan*Ci));...
    ((Cw*Sraan)+(Sw*Craan*Ci)) ((-Sw*Sraan)+(Cw*Craan*Ci));...
    Sw*Si Cw*Si];

p = a*(1-e^2);%semilatus rectum
h = sqrt(mu*p);%mag of angular momentum
r = p/(1 + (e*cosd(f)));%mag of position vector, meters

Vr = (h*e/p)*sind(f);%radial velocity, meters/sec
Vtheta = h/r;%angular velocity, meters/sec

Xdotstar = (Vr*cosd(f)) - (Vtheta*sind(f));
Ydotstar = (Vr*sind(f)) + (Vtheta*cosd(f));

X = Q*[r*cosd(f) Xdotstar; r*sind(f) Ydotstar];

r = transpose(X(:,1));
v = transpose(X(:,2));
X = [r v];
end