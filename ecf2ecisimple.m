function [T] = ecf2ecisimple(t, ag0)
%ECF2ECISIMPLE Convert from ecf to eci using Greenwich angle
%Simplified conversion by rotating solely Z axis where Z and z axis are
%coincident because of the angular velocity vector direction is constant.

mu = 3.9860044e+14; %m^3/s^2, Earth gravitational parameter
we = (2*pi()/86164);%rad/sec, Earth avg rotational rate
re = 6378137;%meters, spherical Earth radius

ag = we*t +ag0;
T = [cos(ag) -sin(ag) 0; sin(ag) cos(ag) 0; 0 0 1];

end