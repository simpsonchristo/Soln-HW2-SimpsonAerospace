function [oe, T, rp, ra, f] = rv2oe(rv, mu)
%RV2OE Convert between orbital elements and pos. and velocity
%Given a known gravitational parameter for a planet, an artificial
%satellite initial state is provided in a nonrotation coordinate system.

rholder = rv(1,1:3);
vholder = rv(1,4:6);
clear rv;
rv(1,:) = rholder;
rv(2,:) = vholder;
%magnitudes
r0 = sqrt(dot(rv(1,:),rv(1,:)));
v0 = sqrt(dot(rv(2,:),rv(2,:)));
%angular momentum, h
h = cross(rv(1,:),rv(2,:));
hmag = sqrt(dot(h,h));
%inclination, i
i = acos(h(1,3)/sqrt(dot(h,h)));

%raan
hxy = sqrt(h(1,1)^2 + h(1,2)^2);
sinraan = h(1,1)/hxy;
cosraan = -h(1,2)/hxy;
raan = atan2(sinraan,cosraan);

%specific energy per unit mass
energy = ((v0^2)/2) - (mu/r0);

%semimajor axis, a
a = -mu/(2*energy);

%eccentricity, e
e = sqrt(1 + (2*energy*dot(h,h))/(mu^2));
evec = cross(rv(2,:),h/mu) - rv(1,:)/r0;

%semiminor axis, b
b = a*sqrt(1-e^2);

%semi-latus rectum, p
p = a*(1-e^2);

%true anomaly, f
if dot(rv(1,:),rv(2,:))<0
    f = (2*pi()) - acos(dot(evec,rv(1,:))/(e*r0));
else
    f = acos(dot(evec,rv(1,:))/(e*r0));
end

%perifocus
sinwf = rv(1,3)/(r0*sin(i));
coswf = ((rv(1,1)/r0)*cos(raan))+((rv(1,2)/r0)*sin(raan));
wf = atan2(sinwf,coswf);
w = wf-f;

%eccentric anomaly
cosE0 = ((r0/a)*cos(f)) + e;

sinE0 = ((r0/b)*sin(f));
E0 = atan2(sinE0,cosE0);

%mean anomaly
M0 = E0 - e*sin(E0);

%period
T = 2*pi()*sqrt(a^3/mu);
rp = a*(1-e);
ra = a*(1+e);

i    = rad2deg(i);
raan = rad2deg(raan);
w    = rad2deg(w);
M0   = rad2deg(M0);
f    = rad2deg(f);

oe = [a e i raan w M0];
end