clear;
clc;
%p1
t = 3000;%seconds
X0 = [7088580.789 -64.326 920.514; -10.20544809 -522.85385193 7482.075141];
twobodyephemeris(X0,t);