function [t_ssss,t_psss,t_ppss,t_ppps,t_pppp,t_dsss,t_dpss,t_dpps,t_dppp,t_ddps,t_ddpp,t_ddds,t_dddp,t_dddd] = timeShellQuartets(nPrim)

t_ssss = 0;
t_psss = 0;
t_ppss = 0;
t_ppps = 0;
t_pppp = 0;
t_dsss = 0;
t_dpss = 0;
t_dpps = 0;
t_dppp = 0;
t_ddps = 0;
t_ddpp = 0;
t_ddds = 0;
t_dddp = 0;
t_dddd = 0;

RPAValues = rand(nPrim,3);
RPBValues = rand(nPrim,3);
RQCValues = rand(nPrim,3);
RQDValues = rand(nPrim,3);
RWPValues = rand(nPrim,3);
RWQValues = rand(nPrim,3);
pValues = rand(nPrim,1);
qValues = rand(nPrim,1);
ppqValues = rand(nPrim,1);

order = 0;
gSSSSNValues = rand(nPrim,order+1);
tic;
gout = OSssss(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
t_ssss = t_ssss + toc;

order = 3;
gSSSSNValues = rand(nPrim,order+1);
tic;
gout = OSdpss(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
t_dpss = t_dpss + toc;

order = 4;
gSSSSNValues = rand(nPrim,order+1);
tic;
gout = OSpppp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
t_pppp = t_pppp + toc;

order = 8;
gSSSSNValues = rand(nPrim,order+1);
tic;
gout = OSdddd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
t_dddd = t_dddd + toc;









end
