function [T] = spherical2ecf(X)
%SPHERICAL2ECF Return transformation matrix for spherical to ECF
%Provided some position vector in spherical coord., the conversion matrix
%between the spherical (radial, latitude, and longitude) and earth centered
%fixed frame is returned.
rstar = X(1,1:3);%radial, latitude, and longitude, meters, rad, and rad

Cphi = cos(rstar(1));%cos(phi) or cosine of latitude
Sphi = sin(rstar(1));%sin(phi) or sine of latitude

Clon = cos(rstar(2));%cos(lambda) or cosine of longitude
Slon = sin(rstar(2));%sin(lambda) or sine of longitude

T = [Cphi*Clon, -Sphi*Clon, -Slon; Cphi*Slon, -Sphi*Slon, Clon; Sphi, Cphi, 0];

end