function gout = OSdsss(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues)
s_s_s_s_0 = gSSSSNValues(:,1);
s_s_s_s_1 = gSSSSNValues(:,2);
s_s_s_s_2 = gSSSSNValues(:,3);
PAx = RPAValues(:,1);
PAy = RPAValues(:,2);
PAz = RPAValues(:,3);
PBx = RPBValues(:,1);
PBy = RPBValues(:,2);
PBz = RPBValues(:,3);
QCx = RQCValues(:,1);
QCy = RQCValues(:,2);
QCz = RQCValues(:,3);
QDx = RQDValues(:,1);
QDy = RQDValues(:,2);
QDz = RQDValues(:,3);
WPx = RWPValues(:,1);
WPy = RWPValues(:,2);
WPz = RWPValues(:,3);
WQx = RWQValues(:,1);
WQy = RWQValues(:,2);
WQz = RWQValues(:,3);
ABx = PBx(1)-PAx(1);
ABy = PBy(1)-PAy(1);
ABz = PBz(1)-PAz(1);
CDx = QDx(1)-QCx(1);
CDy = QDy(1)-QCy(1);
CDz = QDz(1)-QCz(1);
oo2p = 0.5./pValues;
oo2q = 0.5./qValues;
qoppq = qValues./ppqValues;
oo2pq = 0.5./ppqValues;
%[ps|ss](0)
px_s_s_s_0 = PAx.*s_s_s_s_0 + WPx.*s_s_s_s_1;
py_s_s_s_0 = PAy.*s_s_s_s_0 + WPy.*s_s_s_s_1;
pz_s_s_s_0 = PAz.*s_s_s_s_0 + WPz.*s_s_s_s_1;
%[ps|ss](1)
px_s_s_s_1 = PAx.*s_s_s_s_1 + WPx.*s_s_s_s_2;
py_s_s_s_1 = PAy.*s_s_s_s_1 + WPy.*s_s_s_s_2;
pz_s_s_s_1 = PAz.*s_s_s_s_1 + WPz.*s_s_s_s_2;
%Contracted (ps|ss)
px_s_s_s = sum(px_s_s_s_0);
py_s_s_s = sum(py_s_s_s_0);
pz_s_s_s = sum(pz_s_s_s_0);
%[ds|ss](0)
dxx_s_s_s_0 = PAx.*px_s_s_s_0 + WPx.*px_s_s_s_1 + 1.*oo2p.*(s_s_s_s_0 - qoppq.*s_s_s_s_1);
dxy_s_s_s_0 = PAx.*py_s_s_s_0 + WPx.*py_s_s_s_1;
dxz_s_s_s_0 = PAx.*pz_s_s_s_0 + WPx.*pz_s_s_s_1;
dyy_s_s_s_0 = PAy.*py_s_s_s_0 + WPy.*py_s_s_s_1 + 1.*oo2p.*(s_s_s_s_0 - qoppq.*s_s_s_s_1);
dyz_s_s_s_0 = PAy.*pz_s_s_s_0 + WPy.*pz_s_s_s_1;
dzz_s_s_s_0 = PAz.*pz_s_s_s_0 + WPz.*pz_s_s_s_1 + 1.*oo2p.*(s_s_s_s_0 - qoppq.*s_s_s_s_1);
%Contracted (ds|ss)
dxx_s_s_s = sum(dxx_s_s_s_0);
dxy_s_s_s = sum(dxy_s_s_s_0);
dxz_s_s_s = sum(dxz_s_s_s_0);
dyy_s_s_s = sum(dyy_s_s_s_0);
dyz_s_s_s = sum(dyz_s_s_s_0);
dzz_s_s_s = sum(dzz_s_s_s_0);
gout = zeros(6,1,1,1);
gout(1,1,1,1) = dxx_s_s_s;
gout(2,1,1,1) = dxy_s_s_s;
gout(3,1,1,1) = dxz_s_s_s;
gout(4,1,1,1) = dyy_s_s_s;
gout(5,1,1,1) = dyz_s_s_s;
gout(6,1,1,1) = dzz_s_s_s;
end