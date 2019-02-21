function [v1, v2] = lambertprob(r1, t1, r2, t2)
%LAMBERTPROB Solve Lambert's problem to determine velocity at t1 and t2
%Provided the position vector at two unique times determine the velocity at
%each point from the time-of-flight and change in position vectors.
% INPUTS:
%   r1 - Given position in ITRF frame at time t1. [meters]
%   t1 - Time when position vector r1 is known. [days] JD in UT1
%   r2 - Given position in ITRF frame at time t2. [meters]
%   t2 - Time when position vector r2 is known. [days] JD in UT1
% OUTPUTS:
%   v1 - Velocity vector in ITRF at time t1. [meters/sec]
%   v2 - Velocity vector in ITRF at time t2. [meters/sec]

u12 = cross(r1,r2);
if u12(3)<0
    theta12 = (2*pi) - acos(dot(r1,r2)/(norm(r1)*norm(r2)));
else
    theta12 = acos(dot(r1,r2)/(norm(r1)*norm(r2)));
end

A = sin(theta12)*sqrt((norm(r1)*norm(r2))/(1-cos(theta12)));
