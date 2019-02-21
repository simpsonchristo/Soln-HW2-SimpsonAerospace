function [d2Xdt2] = j2potential(t, X)
%J2POTENTIAL Return 2nd derivative of motion for orbit
%Considering only J2 perturbations and point mass force use the kinematics
%and return the acceleration.
% INPUTS:
%   X - Given state in J2000 frame at time t. [meters, meters/seconds]
%   t - State will be predicted for each point. [seconds]
% OUTPUTS:
%   d2Xdt2 - acceleration from point mass potential and J2 [meters/sec^2]
mu = 3.9860044e+14; %m^3/s^2, Earth gravitational parameter
we = (2*pi()/86164);%rad/sec, Earth avg rotational rate
re = 6378137;%meters, spherical Earth radius
ag = 0;%Greenwich angle
J2 = 0.001082636;
C2_0 = -0.0004841695;%normalized value of harmonic coefficient
r = norm(X(1:3));

%eom
fspherical = -(mu/(r^3))*X(1:3);
Tecf2eci = ecf2ecisimple(t, ag);
x = transpose(Tecf2eci)*X(1:3);
latlonalt = ecf2spherical(x);
Tsph2ecf = spherical2ecf(latlonalt);

fnsr = 3*mu*(re^2/r^4)*J2*((3*(sin(latlonalt(1))*sin(latlonalt(1)))-1)/2);%radial direction
fnsphi = -mu*(re^2/r^4)*J2*3*sin(latlonalt(1))*cos(latlonalt(1));%latitude
fnslam = 0;%longitude

d2Xdt2(1:3,1) = X(4:6);
d2Xdt2(4:6,1) = fspherical + Tecf2eci*Tsph2ecf*[fnsr; fnsphi; fnslam];

end