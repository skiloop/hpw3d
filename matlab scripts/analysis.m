
%load errDemo;
amp=1e-10;
tw=2e-9;
t0 = 4.0 * tw;
C = 2.99792458E8;
Pi = 3.14159265358979;
mu=4.0*pi*1.0e-7;
dx=50e-3;
dy=50e-3;
dz=50e-3;

dt = dx/2/C; %0.99 / (C * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));
len=1:500;
dist=12*dx;
delay=dist/C;
t=len*dt;
% numeric solve
figure('Name','Original');
yreal=fun(t-delay,dist,amp,t0,tw,mu,C);
plot(t/1e-9,yreal);
%yreal=yreal*0.1317/max(yreal);
% % 
% figure('Name','C float Compare');
% plot(t,yreal,'r',t,opml,'b');
% legend('numeric','pml');
% %
% figure('Name','C++ double Compare');
% plot(t,yreal,'r',t,npml,'b');
% legend('numeric','pml');