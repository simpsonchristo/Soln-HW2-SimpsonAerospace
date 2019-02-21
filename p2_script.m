clear;
clc;
%p2_script
mu = 3.9860044e+14; %m^3/s^2, Earth gravitational parameter
we = (2*pi()/86164);%rad/sec, Earth avg rotational rate
re = 6378137;%meters, spherical Earth radius
ag = 0;%Greenwich angle
J2 = 0.001082636;
C2_0 = -0.0004841695;%normalized value of harmonic coefficient
% %Julian date from USNO
% api = 'https://api.usno.navy.mil/jdconverter?';
% url_0_18081997 = [api 'ID=crsimpson&date=08/18/1997&time=00:00:00.000000'];
% url_3_18081997 = [api 'ID=crsimpson&date=08/18/1997&time=00:00:03.000000'];
% url_0_19081997 = [api 'ID=crsimpson&date=08/19/1997&time=00:00:00.000000'];
% url_3_19081997 = [api 'ID=crsimpson&date=08/19/1997&time=00:00:03.000000'];
% web_t0_18 = webread(url_0_18081997);%JD in UT1
% t0_18 = web_t0_18.data.jd;
% web_t3_18 = webread(url_3_18081997);%JD in UT1
% t3_18 = web_t3_18.data.jd;
% web_t0_19 = webread(url_0_19081997);%JD in UT1
% t0_19 = web_t0_19.data.jd;
% web_t3_19 = webread(url_3_19081997);%JD in UT1
% t3_19 = web_t3_19.data.jd;

%gps reported positions (xSec_DDMMYYYY)
x0_18081997 = [3325396.441 5472597.483 -2057129.050];
x3_18081997 = [3309747.175 5485240.159 -2048664.333];
x0_19081997 = [4389882.255 -4444406.953 -2508462.520];
x3_19081997 = [4402505.030 -4428002.728 -2515303.456];

%determine vector normal to orbital plane
u_18 = cross(x0_18081997/norm(x0_18081997),x3_18081997/norm(x3_18081997));%specific angular momentum unit vector
check_18 = dot(u_18, x0_18081997);
u_19 = cross(x0_19081997/norm(x0_19081997),x3_19081997/norm(x3_19081997));%specific angular momentum unit vector
check_19 = dot(u_19, x0_19081997);
if (abs(check_18) > 1e-9)
    error('n_18 is not orthogonal to plane.')
elseif (abs(check_19) > 1e-9)
    error('n_19 is not orthogonal to plane.')
end

X = [1 0 0]; Y = [0 1 0]; Z = [0 0 1];
%determine inc and raan of orbital plane 18 Aug 1997
inc_18 = acos(dot(u_18/norm(u_18),Z)/(norm(u_18)*norm(Z)));%rad, inclination
linnodes = cross(u_18/norm(u_18),Z);
raan_18 = acos(dot(linnodes,X)/(norm(linnodes)*norm(X)));%rad, right ascension of the ascending node

%determine inc and raan of orbital plane 19 Aug 1997
inc_19 = acos(dot(u_19/norm(u_19),Z)/(norm(u_19)*norm(Z)));%rad, inclination
linnodes = cross(u_19/norm(u_19),Z);
raan_19 = acos(dot(linnodes,X)/(norm(linnodes)*norm(X)));%rad, right ascension of the ascending node

%computed dAN/dt (deg/day)
estdANdt = (raan_19-raan_18);
% tf = t3_18-t0_18;%day, time to transverse arc
tf = 3/(24*60*60);%day, time to transverse arc
m = 0;%num of orbits traversed
r1 = x0_18081997;%/1000;%km, position vec 1
r2 = x3_18081997;%/1000;%km, position vec 2
[v1,v2, extremal_distances,exitflag] = lambert(r1, r2, tf, m, mu);
if exitflag==-1
    error('The given problem has no solution and cannot be solved.')
elseif exitflag==-2
    error('Both algorithms failed to find a solution. Check your inputs.');
end
[oe] = rv2oe([r1 v1], mu);
n = sqrt(mu/oe(1,1)^3);%mean motion
dANdt = (-3/2)*J2*(n/((1-oe(1,2)^2)^2))*((re/oe(1,1))^2)*cosd(oe(1,3));
dANdt = dANdt*(24*60*60);

