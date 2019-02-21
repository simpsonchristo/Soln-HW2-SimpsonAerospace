clear;
clc;
%example 2.2.6.2 (p. 43)
%illustrate nature of orbit error
t = 0:60:11000;%seconds
X0 = [5492000.34 3984001.40 2955.81 -3931.046491 5498.676921 3665.980697];
[Xt] = twobodyephemeris(X0,t);
%transform to rtn
T = eci2rtn(X0);
RTN = T*transpose(X0(1,1:3));
%add 1 meter error in radial
RTN(1,1) = RTN(1,1)+1;
X0(1,1:3) = transpose(T)*RTN;
[Xerr] = twobodyephemeris(X0,t);

for i=1:length(Xt)
    T = eci2rtn(Xt(i,:));
    RTN(:,i) = T*transpose(Xt(i,1:3));
    T = eci2rtn(Xerr(i,:));
    RTNerr(:,i) = T*transpose(Xerr(i,1:3));
end

h = figure('Name','RTN error');
xlabel('Error (m)');
ylabel('Time (sec)');
plot(t,(RTN(1,:)-RTNerr(1,:)),'k-','DisplayName','R');
hold on
plot(t,(RTN(2,:)-RTNerr(2,:)),'b--','DisplayName','T');
plot(t,(RTN(3,:)-RTNerr(3,:)),'g-o','DisplayName','N');
legend()
saveas(h,'ex2262_rtnerr.png');