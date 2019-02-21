function [Xt, xt, latlonalt] = twobodyephemeris(X0, t)
%TWOBODYEPHEMERIS Produce two-body ephemeris table from time and state
%Using two-body point force, propagate from initial state through given
%time to produce ephemeris. 
% INPUTS:
%   X0 - Given state in inertial frame at time t. [meters, meters/seconds]
%   t  - State will be predicted for each point. [seconds]

mu = 3.9860044e+14; %m^3/s^2, Earth gravitational parameter
we = (2*pi()/86164);%rad/sec, Earth avg rotational rate
re = 6378137;%meters, spherical Earth radius
ag = 0;%Greenwich angle

[oe, T, rp, ra, f0] = rv2oe(X0, mu);
a    = oe(1);%semimajor axis, meters
e    = oe(2);%eccentricity
i    = oe(3);%inclination, degrees
Si = sind(i); Ci = cosd(i);
raan = oe(4);%ascending node, degrees
Sraan = sind(raan); Craan = cosd(raan);
w    = oe(5);%aop, degrees
Sw = sind(w); Cw = cosd(w);
M0    = oe(6);%mean anomaly, degrees
Sm0 = sind(M0); Cm0 = cosd(M0);



n = sqrt(mu/a^3);%mean motion
% fprintf('ECI\n')
% fprintf('t(min)   X(m)   Y(m)     Z(m)\n')
for i=1:length(t)
    Et = keplerseqn(e,n,t(i), deg2rad(M0));

    Et = rad2deg(Et);
    Mt = Et - (e*sind(Et));
    ft = acosd((cosd(Et) - e)/(1 - (e*cosd(Et))));
    oe(6) = Mt;
    Xt(i,:) = oe2rv(oe,ft);
%     fprintf('%.0f   %.0f   %.0f   %.0f\n',(t(i)/60), Xt(i,1), Xt(i,2), Xt(i,3));
end
%geocentric coordinate conversion
% fprintf('ECEF\n')
% fprintf('t(min)   x(m)   y(m)     z(m)\n')
for i=1:length(t)
    ag = we*t(i);
    rotz = [cos(ag) -sin(ag) 0; sin(ag) cos(ag) 0; 0 0 1];
    xt(i,:) = Xt(i,1:3)*rotz;
%     fprintf('%.0f   %.0f   %.0f   %.0f\n',(t(i)/60), xt(i,1), xt(i,2), xt(i,3));
end
% fprintf('Geocentric lat, lon, and height\n')
% fprintf('t(min)   phi(deg)   lambda(deg)     h(m)\n')
for i=1:length(t)
    r = sqrt(xt(i,1)^2 + xt(i,2)^2 + xt(i,3)^2);
    phi(i) = asin(xt(i,3)/r);
    lambda(i) = atan2(xt(i,2)/(r*cos(phi(i))),xt(i,1)/(r*cos(phi(i))));
    h(i) = r - re;
    latlonalt(i,:) = [phi(i) lambda(i) h(i)];
%     fprintf('%.0f       %.3f     %.3f       %.2f\n',(t(i)/60), rad2deg(phi(i)), rad2deg(lambda(i)), h(i));
end
