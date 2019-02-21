clear;
clc;
%example 2.2.5 (p. 37)
t = 30*60:2*60:34*60;%seconds
%Et = 159.628138;%eccentric anomaly at time t, degrees
X0 = [5492000.34 3984001.40 2955.81; -3931.046491 5498.676921 3665.980697];
twobodyephemeris(X0,t);
