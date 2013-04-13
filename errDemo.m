%% load data
% load new dat
cd /home/skiloop/workspace/hpw3d/new
loadAllDat
cd /home/skiloop/matlab/hpw3d
save new
clear
% load old dat
cd /home/skiloop/workspace/hpw3d/old
loadAllDat
cd /home/skiloop/matlab/hpw3d
save old
clear
% load npml and opml
load old opml
load new npml
%% compare
% compute analysis results
amp=1000;
tw=53e-12;
t0 = 4.0 * tw;
C = 2.99792458E8;
Pi = 3.14159265358979;
mu=4.0*Pi*1.0e-7;
dx=1e-3;
dy=1e-3;
dz=1e-3;

dt = 0.99 / (C * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));
len=1:500;
dist=12E-3;
delay=dist/C;close
t=len*dt;
% numeric solve
yreal=fun(t-delay,dist,amp,t0,tw,mu,C);
%yreal=yreal*0.1317/max(yreal);

%
figure('Name','ErrorComparision');
plot(t,yreal,'r',t,npml,'b+',t,opml,'k.-');
legend('numeric','npml','opml');