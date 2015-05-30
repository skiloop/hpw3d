
%% load data
clear;
load lar/co.dat
cpml=co;
clear co;
load lar_no/co.dat;
ncpml=co;

%% prepare data
len=min([length(cpml),length(ncpml)]);
if len>1000
    len=1000;
end
cpml=cpml(1:len,:);
ncpml=ncpml(1:len,:);

ez_z_pml=cpml(:,1);
ez_z_no_pml=ncpml(:,1);

err=abs(ez_z_pml-ez_z_no_pml)/max(ez_z_pml);

dt=9.09091e-14;
t=dt*(1:len);

dbError=20*log10(err);

%% make figure
% comparision of Ex in time
figure('NumberTitle','OFF','Name','Time Macthing');
plot(t/1e-9,ez_z_pml,t/1e-9,ez_z_no_pml,'--','LineWidth',2);
xlabel('time (ns)');
ylabel('Ez');
grid on;
title('time domain compare');
legend('with cpml','reference results');
print -depsc -tiff -r300 cpml3dTimeCompare_r
print -dtiff -r300 cpml3dTimeCompare_r
%% Relative Error
figure('NumberTitle','OFF','Name','Relative Error');
semilogy(t/1e-9,err,'LineWidth',2);
xlabel('time (ns)');
ylabel('Relative Error');
grid on;
title('Relative Error');
print -depsc -tiff -r300 cpml3dRelativeError_r
print -dtiff -r300 cpml3dRelativeError_r
% Error in DB
figure('NumberTitle','OFF','Name','DB Error');
plot(t/1e-9,dbError,'LineWidth',2);
xlabel('time (ns)');
ylabel('Error (DB)');
grid on;
title('DB Error');
print -depsc -tiff -r300 cpml3dDBError_r
print -dtiff -r300 cpml3dDBError_r


