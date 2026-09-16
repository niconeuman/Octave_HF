%This script times the different angular momentum shell quartet calculations

nClasses = 14; %classes of shell quartet integrals
nPrimMax = 4000;
timeTable = zeros(nPrimMax/4,nClasses);

for k = 1:nPrimMax/4
    nPrim = 4*(k-1)+1;
    [t_ssss,t_psss,t_ppss,t_ppps,t_pppp,t_dsss,t_dpss,t_dpps,t_dppp,t_ddps,t_ddpp,t_ddds,t_dddp,t_dddd] = timeShellQuartets(nPrim);
    timeTable(k,:) = [t_ssss,t_psss,t_ppss,t_ppps,t_pppp,t_dsss,t_dpss,t_dpps,t_dppp,t_ddps,t_ddpp,t_ddds,t_dddp,t_dddd];

end

profile clear;
profile on;
[t_ssss,t_psss,t_ppss,t_ppps,t_pppp,t_dsss,t_dpss,t_dpps,t_dppp,t_ddps,t_ddpp,t_ddds,t_dddp,t_dddd] = timeShellQuartets(nPrim);
profile off;

plot(4*(1:nPrimMax/4)',timeTable);
legend('t_{ssss} (s)','t_{pppp} (s)','t_{dpss} (s)','t_{dddd} (s)');
xlabel('nPrim');
ylabel('t (s)');
