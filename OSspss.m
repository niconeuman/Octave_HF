function gout = OSspss(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues)
s_s_s_s_0 = gSSSSNValues(:,1);
s_s_s_s_1 = gSSSSNValues(:,2);
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
%Contracted (ps|ss)
px_s_s_s = sum(px_s_s_s_0);
py_s_s_s = sum(py_s_s_s_0);
pz_s_s_s = sum(pz_s_s_s_0);
%(sp|ss) from horizontal recursion on the bra 
s_px_s_s = ABx*s_s_s_s + px_s_s_s;
s_py_s_s = ABy*s_s_s_s + py_s_s_s;
s_pz_s_s = ABz*s_s_s_s + pz_s_s_s;
gout = zeros(1,3,1,1);
gout(1,1,1,1) = s_px_s_s;
gout(1,2,1,1) = s_py_s_s;
gout(1,3,1,1) = s_pz_s_s;
end