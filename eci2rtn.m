function [T] = eci2rtn(X)
%ECI2RTN Return transformation matrix for ECI to RTN
%Provided some reference orbit, rstar and vstar, the conversion matrix
%between the radial, transverse, and normal (radial, along-track,
%cross-track) and earth centered inertial frame is returned.
rstar = X(1,1:3);%position vector, meters
vstar = X(1,4:6);%velocity vector, meters/second

ur = rstar/norm(rstar);
un = cross(rstar,vstar)/norm(cross(rstar,vstar));
ut = cross(un,ur);

T = [ur(1) ur(2) ur(3); un(1) un(2) un(3); ut(1) ut(2) ut(3)];

end