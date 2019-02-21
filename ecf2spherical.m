function [latlonalt] = ecf2spherical(xt)
%ECF2SPHERICAL Convert ITRF to geocentric coordinates
%Convert position vector in earth fixed frame to latitude, longitude, and
%altitude based coordinates or a spherical system. 
% INPUTS:
%   X0 - Given state in ITRF frame at time t. [meters, meters/seconds]
% OUTPUTS:
%   latlonalt - Geocentric latitude, longitude, and altitude. [rad, rad,
%   meters]

re = 6378137;%meters, spherical Earth radius
r = sqrt(xt(1)^2 + xt(2)^2 + xt(3)^2);
phi = asin(xt(3)/r);
lambda = atan2(xt(2)/(r*cos(phi)),xt(1)/(r*cos(phi)));
h = r - re;
latlonalt(1,:) = [phi lambda h];


end
