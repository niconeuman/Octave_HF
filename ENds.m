function Vout = ENds(RPAValues,RPBValues,RPCValues,pValues,VssNValues)
Vs_s_0 = VssNValues(:,1);
Vs_s_1 = VssNValues(:,2);
Vs_s_2 = VssNValues(:,3);
PAx = RPAValues(:,1);
PAy = RPAValues(:,2);
PAz = RPAValues(:,3);
PCx = RPCValues(:,1);
PCy = RPCValues(:,2);
PCz = RPCValues(:,3);
PBx = RPBValues(:,1);
PBy = RPBValues(:,2);
PBz = RPBValues(:,3);
oo2p = 0.5./pValues;
ABx = PBx(1)-PAx(1);
ABy = PBy(1)-PAy(1);
ABz = PBz(1)-PAz(1);
%[p|VeN|s](0)
Vpx_s_0 = PAx.*Vs_s_0 - PCx.*Vs_s_1;
Vpy_s_0 = PAy.*Vs_s_0 - PCy.*Vs_s_1;
Vpz_s_0 = PAz.*Vs_s_0 - PCz.*Vs_s_1;
%[p|VeN|s](1)
Vpx_s_1 = PAx.*Vs_s_1 - PCx.*Vs_s_2;
Vpy_s_1 = PAy.*Vs_s_1 - PCy.*Vs_s_2;
Vpz_s_1 = PAz.*Vs_s_1 - PCz.*Vs_s_2;
%Contracted (p|VeN|s)
Vpx_s = sum(Vpx_s_0);
Vpy_s = sum(Vpy_s_0);
Vpz_s = sum(Vpz_s_0);
%[d|VeN|s](0)
Vdxx_s_0 = PAx.*Vpx_s_0 - PCx.*Vpx_s_1 + 1.*oo2p.*(Vs_s_0 - Vs_s_1);
Vdxy_s_0 = PAx.*Vpy_s_0 - PCx.*Vpy_s_1;
Vdxz_s_0 = PAx.*Vpz_s_0 - PCx.*Vpz_s_1;
Vdyy_s_0 = PAy.*Vpy_s_0 - PCy.*Vpy_s_1 + 1.*oo2p.*(Vs_s_0 - Vs_s_1);
Vdyz_s_0 = PAy.*Vpz_s_0 - PCy.*Vpz_s_1;
Vdzz_s_0 = PAz.*Vpz_s_0 - PCz.*Vpz_s_1 + 1.*oo2p.*(Vs_s_0 - Vs_s_1);
%Contracted (d|VeN|s)
Vdxx_s = sum(Vdxx_s_0);
Vdxy_s = sum(Vdxy_s_0);
Vdxz_s = sum(Vdxz_s_0);
Vdyy_s = sum(Vdyy_s_0);
Vdyz_s = sum(Vdyz_s_0);
Vdzz_s = sum(Vdzz_s_0);
V = zeros(6,1);
V(1,1) = 1/1.73205080756888*Vdxx_s;
V(2,1) = Vdxy_s;
V(3,1) = Vdxz_s;
V(4,1) = 1/1.73205080756888*Vdyy_s;
V(5,1) = Vdyz_s;
V(6,1) = 1/1.73205080756888*Vdzz_s;
Vout = V;end
