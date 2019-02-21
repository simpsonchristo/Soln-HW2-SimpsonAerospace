function [Xt, oe, tn] = J2perturbedephemeris(X0, t)
%J2PERTURBEDEPHEMERIS Produce ephemeris using EOM with J2 perturbation
%Using point mass and J2 perturbation, propagate from initial state through
%given time to produce the ephemeris.
% INPUTS:
%   X0 - Given state in J2000 frame at time t. [meters, meters/seconds]
%   t  - State will be predicted for each point. [seconds]
% OUTPUTS:

mu = 3.9860044e+14; %m^3/s^2, Earth gravitational parameter
% we = (2*pi()/86164);%rad/sec, Earth avg rotational rate
% re = 6378137;%meters, spherical Earth radius
% ag = 0;%Greenwich angle
% J2 = 0.001082636;
% C2_0 = -0.0004841695;%normalized value of harmonic coefficient

%eom
[tn, Xt] = ode45(@j2potential,[t(1) t(length(t))], X0);

for i=1:length(Xt)
    oe(i,:) = rv2oe(Xt(i,:), mu);
end
end