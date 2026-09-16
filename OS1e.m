function [Sab,Tab] = OS1e(KabValues,RPAValues,RPBValues,aValues,bValues,pValues,WeightValues,R12sq,La,Lb)

%This function calculates one electron integrals (overlap and kinetic energy)
%R12sq is a function only on the atoms, not the contracted functions
L1 = La;
L2 = Lb;

PAx = RPAValues(:,1);
PAy = RPAValues(:,2);
PAz = RPAValues(:,3);

PBx = RPBValues(:,1);
PBy = RPBValues(:,2);
PBz = RPBValues(:,3);

oo2p = 0.5./pValues;

if (L1 == 0 && L2 == 0) %[s|s]
    Sab = WeightValues.*KabValues.*(pi./pValues).^1.5;
    ab = aValues.*bValues;
    abop = ab./pValues;
    Tab = Sab.*(abop.*(3-2*abop.*R12sq));

    Sab = sum(Sab);
    Tab = sum(Tab);

elseif (L1 == 1 && L2 == 0) %[p|s]

    Ss_s = WeightValues.*KabValues.*(pi./pValues).^1.5;
    abop = aValues.*bValues./pValues;
    Ts_s= Ss_s.*(abop.*(3-2*abop.*R12sq));

    Spx_s = PAx.*Ss_s;
    Spy_s = PAy.*Ss_s;
    Spz_s = PAz.*Ss_s;

    Tpx_s = PAx.*Ts_s+2*abop.*(Spx_s);
    Tpy_s = PAy.*Ts_s+2*abop.*(Spy_s);
    Tpz_s = PAz.*Ts_s+2*abop.*(Spz_s);

    Sab = zeros(3,1);
    Sab(1) = sum(Spx_s);
    Sab(2) = sum(Spy_s);
    Sab(3) = sum(Spz_s);

    Tab = zeros(3,1);
    Tab(1) = sum(Tpx_s);
    Tab(2) = sum(Tpy_s);
    Tab(3) = sum(Tpz_s);

elseif (L1 == 2 && L2 == 0) %[d|s]

    Ss_s = WeightValues.*KabValues.*(pi./pValues).^1.5;
    abop = aValues.*bValues./pValues;
    oo2p = 0.5./pValues;
    Ts_s= Ss_s.*(abop.*(3-2*abop.*R12sq));

    Spx_s = PAx.*Ss_s;
    Spy_s = PAy.*Ss_s;
    Spz_s = PAz.*Ss_s;

    Tpx_s = PAx.*Ts_s+2*abop.*(Spx_s);
    Tpy_s = PAy.*Ts_s+2*abop.*(Spy_s);
    Tpz_s = PAz.*Ts_s+2*abop.*(Spz_s);

    Sdxx_s = PAx.*Spx_s+1*oo2p.*Ss_s;
    Sdxy_s = PAy.*Spx_s;
    Sdxz_s = PAz.*Spx_s;
    Sdyy_s = PAy.*Spy_s+1*oo2p.*Ss_s;
    Sdyz_s = PAz.*Spy_s;
    Sdzz_s = PAz.*Spz_s+1*oo2p.*Ss_s;

    Tdxx_s = PAx.*Tpx_s+2*abop.*(Sdxx_s+1*oo2p.*Ss_s)+1*oo2p.*Ts_s;
    Tdxy_s = PAy.*Tpx_s+2*abop.*(Sdxy_s);
    Tdxz_s = PAz.*Tpx_s+2*abop.*(Sdxz_s);
    Tdyy_s = PAy.*Tpy_s+2*abop.*(Sdyy_s+1*oo2p.*Ss_s)+1*oo2p.*Ts_s;
    Tdyz_s = PAz.*Tpy_s+2*abop.*(Sdyz_s);
    Tdzz_s = PAz.*Tpz_s+2*abop.*(Sdzz_s+1*oo2p.*Ss_s)+1*oo2p.*Ts_s;

    Sab = zeros(6,1);
    Sab(1) = sum(Sdxx_s);
    Sab(2) = sum(Sdxy_s);
    Sab(3) = sum(Sdxz_s);
    Sab(4) = sum(Sdyy_s);
    Sab(5) = sum(Sdyz_s);
    Sab(6) = sum(Sdzz_s);

    Tab = zeros(6,1);
    Tab(1) = sum(Tdxx_s);
    Tab(2) = sum(Tdxy_s);
    Tab(3) = sum(Tdxz_s);
    Tab(4) = sum(Tdyy_s);
    Tab(5) = sum(Tdyz_s);
    Tab(6) = sum(Tdzz_s);

elseif (L1 == 0 && L2 == 1) %[s|p]

    Ss_s = WeightValues.*KabValues.*(pi./pValues).^1.5;
    abop = aValues.*bValues./pValues;
    Ts_s= Ss_s.*(abop.*(3-2*abop.*R12sq));
%Ss_p
    Ss_px = PBx.*Ss_s;
    Ss_py = PBy.*Ss_s;
    Ss_pz = PBz.*Ss_s;
%Ts_p
    Ts_px = PBx.*Ts_s+2*abop.*(Ss_px);
    Ts_py = PBy.*Ts_s+2*abop.*(Ss_py);
    Ts_pz = PBz.*Ts_s+2*abop.*(Ss_pz);

    Sab = zeros(1,3);
    Sab(1) = sum(Ss_px);
    Sab(2) = sum(Ss_py);
    Sab(3) = sum(Ss_pz);

    Tab = zeros(1,3);
    Tab(1) = sum(Ts_px);
    Tab(2) = sum(Ts_py);
    Tab(3) = sum(Ts_pz);

elseif (L1 == 0 && L2 == 2) %[s|d]

    Ss_s = WeightValues.*KabValues.*(pi./pValues).^1.5;
    abop = aValues.*bValues./pValues;
    oo2p = 0.5./pValues;
    Ts_s= Ss_s.*(abop.*(3-2*abop.*R12sq));

    Ss_px = PBx.*Ss_s;
    Ss_py = PBy.*Ss_s;
    Ss_pz = PBz.*Ss_s;

    Ts_px = PBx.*Ts_s+2*abop.*(Ss_px);
    Ts_py = PBy.*Ts_s+2*abop.*(Ss_py);
    Ts_pz = PBz.*Ts_s+2*abop.*(Ss_pz);

    Ss_dxx = PBx.*Ss_px+1*oo2p.*Ss_s;
    Ss_dxy = PBy.*Ss_px;
    Ss_dxz = PBz.*Ss_px;
    Ss_dyy = PBy.*Ss_py+1*oo2p.*Ss_s;
    Ss_dyz = PBz.*Ss_py;
    Ss_dzz = PBz.*Ss_pz+1*oo2p.*Ss_s;

    Ts_dxx = PBx.*Ts_px+2*abop.*(Ss_dxx+1*oo2p.*Ss_s)+1*oo2p.*Ts_s;
    Ts_dxy = PBy.*Ts_px+2*abop.*(Ss_dxy);
    Ts_dxz = PBz.*Ts_px+2*abop.*(Ss_dxz);
    Ts_dyy = PBy.*Ts_py+2*abop.*(Ss_dyy+1*oo2p.*Ss_s)+1*oo2p.*Ts_s;
    Ts_dyz = PBz.*Ts_py+2*abop.*(Ss_dyz);
    Ts_dzz = PBz.*Ts_pz+2*abop.*(Ss_dzz+1*oo2p.*Ss_s)+1*oo2p.*Ts_s;

    Sab = zeros(1,6);
    Sab(1) = sum(Ss_dxx);
    Sab(2) = sum(Ss_dxy);
    Sab(3) = sum(Ss_dxz);
    Sab(4) = sum(Ss_dyy);
    Sab(5) = sum(Ss_dyz);
    Sab(6) = sum(Ss_dzz);

    Tab = zeros(1,6);
    Tab(1) = sum(Ts_dxx);
    Tab(2) = sum(Ts_dxy);
    Tab(3) = sum(Ts_dxz);
    Tab(4) = sum(Ts_dyy);
    Tab(5) = sum(Ts_dyz);
    Tab(6) = sum(Ts_dzz);

elseif (L1 == 1 && L2 == 1) %[p|p]

    Ss_s = WeightValues.*KabValues.*(pi./pValues).^1.5;
    abop = aValues.*bValues./pValues;
    oo2p = 0.5./pValues;
    Ts_s= Ss_s.*(abop.*(3-2*abop.*R12sq));

    Spx_s = PAx.*Ss_s;
    Spy_s = PAy.*Ss_s;
    Spz_s = PAz.*Ss_s;

    Tpx_s = PAx.*Ts_s+2*abop.*(Spx_s);
    Tpy_s = PAy.*Ts_s+2*abop.*(Spy_s);
    Tpz_s = PAz.*Ts_s+2*abop.*(Spz_s);

    Spx_px = PBx.*Spx_s+1*oo2p.*Ss_s;
    Spx_py = PBy.*Spx_s;
    Spx_pz = PBz.*Spx_s;

    Spy_px = PBx.*Spy_s;
    Spy_py = PBy.*Spy_s+1*oo2p.*Ss_s;
    Spy_pz = PBz.*Spy_s;

    Spz_px = PBx.*Spz_s;
    Spz_py = PBy.*Spz_s;
    Spz_pz = PBz.*Spz_s+1*oo2p.*Ss_s;

    Tpx_px = PBx.*Tpx_s+1*oo2p.*Ts_s+2*abop.*(Spx_px);
    Tpx_py = PBy.*Tpx_s             +2*abop.*(Spx_py);
    Tpx_pz = PBz.*Tpx_s             +2*abop.*(Spx_pz);

    Tpy_px = PBx.*Tpy_s             +2*abop.*(Spy_px);
    Tpy_py = PBy.*Tpy_s+1*oo2p.*Ts_s+2*abop.*(Spy_py);
    Tpy_pz = PBz.*Tpy_s             +2*abop.*(Spy_pz);

    Tpz_px = PBx.*Tpz_s             +2*abop.*(Spz_px);
    Tpz_py = PBy.*Tpz_s             +2*abop.*(Spz_py);
    Tpz_pz = PBz.*Tpz_s+1*oo2p.*Ts_s+2*abop.*(Spz_pz);

    Sab = zeros(3,3);
    Sab(1,1) = sum(Spx_px);
    Sab(1,2) = sum(Spx_py);
    Sab(1,3) = sum(Spx_pz);

    Sab(2,1) = sum(Spy_px);
    Sab(2,2) = sum(Spy_py);
    Sab(2,3) = sum(Spy_pz);

    Sab(3,1) = sum(Spz_px);
    Sab(3,2) = sum(Spz_py);
    Sab(3,3) = sum(Spz_pz);

    Tab = zeros(3,3);
    Tab(1,1) = sum(Tpx_px);
    Tab(1,2) = sum(Tpx_py);
    Tab(1,3) = sum(Tpx_pz);

    Tab(2,1) = sum(Tpy_px);
    Tab(2,2) = sum(Tpy_py);
    Tab(2,3) = sum(Tpy_pz);

    Tab(3,1) = sum(Tpz_px);
    Tab(3,2) = sum(Tpz_py);
    Tab(3,3) = sum(Tpz_pz);
elseif (L1 == 1 && L2 == 2) %[p|d]

    Ss_s = WeightValues.*KabValues.*(pi./pValues).^1.5;
    abop = aValues.*bValues./pValues;
    oo2p = 0.5./pValues;
    Ts_s= Ss_s.*(abop.*(3-2*abop.*R12sq));
%Sp_s
    Spx_s = PAx.*Ss_s;
    Spy_s = PAy.*Ss_s;
    Spz_s = PAz.*Ss_s;
%Tp_s
    Tpx_s = PAx.*Ts_s+2*abop.*(Spx_s);
    Tpy_s = PAy.*Ts_s+2*abop.*(Spy_s);
    Tpz_s = PAz.*Ts_s+2*abop.*(Spz_s);
%Ss_p
    Ss_px = PBx.*Ss_s;
    Ss_py = PBy.*Ss_s;
    Ss_pz = PBz.*Ss_s;
%Ts_p
    Ts_px = PBx.*Ts_s+2*abop.*(Ss_px);
    Ts_py = PBy.*Ts_s+2*abop.*(Ss_py);
    Ts_pz = PBz.*Ts_s+2*abop.*(Ss_pz);
%Sp_p
    Spx_px = PBx.*Spx_s+1*oo2p.*Ss_s;
    Spx_py = PBy.*Spx_s;
    Spx_pz = PBz.*Spx_s;

    Spy_px = PBx.*Spy_s;
    Spy_py = PBy.*Spy_s+1*oo2p.*Ss_s;
    Spy_pz = PBz.*Spy_s;

    Spz_px = PBx.*Spz_s;
    Spz_py = PBy.*Spz_s;
    Spz_pz = PBz.*Spz_s+1*oo2p.*Ss_s;
%Tp_p
    Tpx_px = PBx.*Tpx_s+1*oo2p.*Ts_s+2*abop.*(Spx_px);
    Tpx_py = PBy.*Tpx_s             +2*abop.*(Spx_py);
    Tpx_pz = PBz.*Tpx_s             +2*abop.*(Spx_pz);

    Tpy_px = PBx.*Tpy_s             +2*abop.*(Spy_px);
    Tpy_py = PBy.*Tpy_s+1*oo2p.*Ts_s+2*abop.*(Spy_py);
    Tpy_pz = PBz.*Tpy_s             +2*abop.*(Spy_pz);

    Tpz_px = PBx.*Tpz_s             +2*abop.*(Spz_px);
    Tpz_py = PBy.*Tpz_s             +2*abop.*(Spz_py);
    Tpz_pz = PBz.*Tpz_s+1*oo2p.*Ts_s+2*abop.*(Spz_pz);

%Sp_d
    Spx_dxx = PBx.*Spx_px + 1.*oo2p.*(Spx_s) + 1.*oo2p.*(Ss_px); %I think this could be wrong, could there be a 2 ?
    Spx_dxy = PBx.*Spx_py + 1.*oo2p.*(Ss_py);
    Spx_dxz = PBx.*Spx_pz + 1.*oo2p.*(Ss_pz);
    Spx_dyy = PBy.*Spx_py + 1.*oo2p.*(Spx_s);
    Spx_dyz = PBy.*Spx_pz;
    Spx_dzz = PBz.*Spx_pz + 1.*oo2p.*(Spx_s);
    Spy_dxx = PBx.*Spy_px + 1.*oo2p.*(Spy_s);
    Spy_dxy = PBx.*Spy_py;
    Spy_dxz = PBx.*Spy_pz;
    Spy_dyy = PBy.*Spy_py + 1.*oo2p.*(Spy_s) + 1.*oo2p.*(Ss_py);
    Spy_dyz = PBy.*Spy_pz + 1.*oo2p.*(Ss_pz);
    Spy_dzz = PBz.*Spy_pz + 1.*oo2p.*(Spy_s);
    Spz_dxx = PBx.*Spz_px + 1.*oo2p.*(Spz_s);
    Spz_dxy = PBx.*Spz_py;
    Spz_dxz = PBx.*Spz_pz;
    Spz_dyy = PBy.*Spz_py + 1.*oo2p.*(Spz_s);
    Spz_dyz = PBy.*Spz_pz;
    Spz_dzz = PBz.*Spz_pz + 1.*oo2p.*(Spz_s) + 1.*oo2p.*(Ss_pz);
%Tp_d
    Tpx_dxx = PBx.*Tpx_px + 1.*oo2p.*(Tpx_s) + 1.*oo2p.*(Ts_px) + 2.*abop.*Spx_dxx -2.*abop.*1.*oo2p.*(Spx_s);
    Tpx_dxy = PBx.*Tpx_py + 1.*oo2p.*(Ts_py) + 2.*abop.*Spx_dxy;
    Tpx_dxz = PBx.*Tpx_pz + 1.*oo2p.*(Ts_pz) + 2.*abop.*Spx_dxz;
    Tpx_dyy = PBy.*Tpx_py + 1.*oo2p.*(Tpx_s) + 2.*abop.*Spx_dyy -2.*abop.*1.*oo2p.*(Spx_s);
    Tpx_dyz = PBy.*Tpx_pz + 2.*abop.*Spx_dyz;
    Tpx_dzz = PBz.*Tpx_pz + 1.*oo2p.*(Tpx_s) + 2.*abop.*Spx_dzz -2.*abop.*1.*oo2p.*(Spx_s);
    Tpy_dxx = PBx.*Tpy_px + 1.*oo2p.*(Tpy_s) + 2.*abop.*Spy_dxx -2.*abop.*1.*oo2p.*(Spy_s);
    Tpy_dxy = PBx.*Tpy_py + 2.*abop.*Spy_dxy;
    Tpy_dxz = PBx.*Tpy_pz + 2.*abop.*Spy_dxz;
    Tpy_dyy = PBy.*Tpy_py + 1.*oo2p.*(Tpy_s) + 1.*oo2p.*(Ts_py) + 2.*abop.*Spy_dyy -2.*abop.*1.*oo2p.*(Spy_s);
    Tpy_dyz = PBy.*Tpy_pz + 1.*oo2p.*(Ts_pz) + 2.*abop.*Spy_dyz;
    Tpy_dzz = PBz.*Tpy_pz + 1.*oo2p.*(Tpy_s) + 2.*abop.*Spy_dzz -2.*abop.*1.*oo2p.*(Spy_s);
    Tpz_dxx = PBx.*Tpz_px + 1.*oo2p.*(Tpz_s) + 2.*abop.*Spz_dxx -2.*abop.*1.*oo2p.*(Spz_s);
    Tpz_dxy = PBx.*Tpz_py + 2.*abop.*Spz_dxy;
    Tpz_dxz = PBx.*Tpz_pz + 2.*abop.*Spz_dxz;
    Tpz_dyy = PBy.*Tpz_py + 1.*oo2p.*(Tpz_s) + 2.*abop.*Spz_dyy -2.*abop.*1.*oo2p.*(Spz_s);
    Tpz_dyz = PBy.*Tpz_pz + 2.*abop.*Spz_dyz;
    Tpz_dzz = PBz.*Tpz_pz + 1.*oo2p.*(Tpz_s) + 1.*oo2p.*(Ts_pz) + 2.*abop.*Spz_dzz -2.*abop.*1.*oo2p.*(Spz_s);

    S = zeros(3,6);
    S(1,1) = sum(Spx_dxx);
    S(1,2) = sum(Spx_dxy);
    S(1,3) = sum(Spx_dxz);
    S(1,4) = sum(Spx_dyy);
    S(1,5) = sum(Spx_dyz);
    S(1,6) = sum(Spx_dzz);
    S(2,1) = sum(Spy_dxx);
    S(2,2) = sum(Spy_dxy);
    S(2,3) = sum(Spy_dxz);
    S(2,4) = sum(Spy_dyy);
    S(2,5) = sum(Spy_dyz);
    S(2,6) = sum(Spy_dzz);
    S(3,1) = sum(Spz_dxx);
    S(3,2) = sum(Spz_dxy);
    S(3,3) = sum(Spz_dxz);
    S(3,4) = sum(Spz_dyy);
    S(3,5) = sum(Spz_dyz);
    S(3,6) = sum(Spz_dzz);

    T = zeros(3,6);
    T(1,1) = sum(Tpx_dxx);
    T(1,2) = sum(Tpx_dxy);
    T(1,3) = sum(Tpx_dxz);
    T(1,4) = sum(Tpx_dyy);
    T(1,5) = sum(Tpx_dyz);
    T(1,6) = sum(Tpx_dzz);
    T(2,1) = sum(Tpy_dxx);
    T(2,2) = sum(Tpy_dxy);
    T(2,3) = sum(Tpy_dxz);
    T(2,4) = sum(Tpy_dyy);
    T(2,5) = sum(Tpy_dyz);
    T(2,6) = sum(Tpy_dzz);
    T(3,1) = sum(Tpz_dxx);
    T(3,2) = sum(Tpz_dxy);
    T(3,3) = sum(Tpz_dxz);
    T(3,4) = sum(Tpz_dyy);
    T(3,5) = sum(Tpz_dyz);
    T(3,6) = sum(Tpz_dzz);

    Sab = S;
    Tab = T;
elseif (L1 == 2 && L2 == 1) %[d|p]

    Ss_s = WeightValues.*KabValues.*(pi./pValues).^1.5;
    abop = aValues.*bValues./pValues;
    oo2p = 0.5./pValues;
    Ts_s= Ss_s.*(abop.*(3-2*abop.*R12sq));
%Sp_s
    Spx_s = PAx.*Ss_s;
    Spy_s = PAy.*Ss_s;
    Spz_s = PAz.*Ss_s;
%Tp_s
    Tpx_s = PAx.*Ts_s+2*abop.*(Spx_s);
    Tpy_s = PAy.*Ts_s+2*abop.*(Spy_s);
    Tpz_s = PAz.*Ts_s+2*abop.*(Spz_s);

%Sd_s
    Sdxx_s = PAx.*Spx_s+1*oo2p.*Ss_s;
    Sdxy_s = PAy.*Spx_s;
    Sdxz_s = PAz.*Spx_s;
    Sdyy_s = PAy.*Spy_s+1*oo2p.*Ss_s;
    Sdyz_s = PAz.*Spy_s;
    Sdzz_s = PAz.*Spz_s+1*oo2p.*Ss_s;
%Td_s
    Tdxx_s = PAx.*Tpx_s+2*abop.*(Sdxx_s+1*oo2p.*Ss_s)+1*oo2p.*Ts_s;
    Tdxy_s = PAy.*Tpx_s+2*abop.*(Sdxy_s);
    Tdxz_s = PAz.*Tpx_s+2*abop.*(Sdxz_s);
    Tdyy_s = PAy.*Tpy_s+2*abop.*(Sdyy_s+1*oo2p.*Ss_s)+1*oo2p.*Ts_s;
    Tdyz_s = PAz.*Tpy_s+2*abop.*(Sdyz_s);
    Tdzz_s = PAz.*Tpz_s+2*abop.*(Sdzz_s+1*oo2p.*Ss_s)+1*oo2p.*Ts_s;

%Sd_p
    Sdxx_px = PBx.*Sdxx_s + 2.*oo2p.*(Spx_s);
    Sdxx_py = PBy.*Sdxx_s;
    Sdxx_pz = PBz.*Sdxx_s;

    Sdxy_px = PBx.*Sdxy_s + 1.*oo2p.*(Spy_s);
    Sdxy_py = PBy.*Sdxy_s + 1.*oo2p.*(Spx_s);
    Sdxy_pz = PBz.*Sdxy_s;

    Sdxz_px = PBx.*Sdxz_s + 1.*oo2p.*(Spz_s);
    Sdxz_py = PBy.*Sdxz_s;
    Sdxz_pz = PBz.*Sdxz_s + 1.*oo2p.*(Spx_s);

    Sdyy_px = PBx.*Sdyy_s;
    Sdyy_py = PBy.*Sdyy_s + 2.*oo2p.*(Spy_s);
    Sdyy_pz = PBz.*Sdyy_s;

    Sdyz_px = PBx.*Sdyz_s;
    Sdyz_py = PBy.*Sdyz_s + 1.*oo2p.*(Spz_s);
    Sdyz_pz = PBz.*Sdyz_s + 1.*oo2p.*(Spy_s);

    Sdzz_px = PBx.*Sdzz_s;
    Sdzz_py = PBy.*Sdzz_s;
    Sdzz_pz = PBz.*Sdzz_s + 2.*oo2p.*(Spz_s);
%Td_p
    Tdxx_px = PBx.*Tdxx_s + 2.*oo2p.*(Tpx_s) + 2.*abop.*Sdxx_px;
    Tdxx_py = PBy.*Tdxx_s + 2.*abop.*Sdxx_py;
    Tdxx_pz = PBz.*Tdxx_s + 2.*abop.*Sdxx_pz;

    Tdxy_px = PBx.*Tdxy_s + 1.*oo2p.*(Tpy_s) + 2.*abop.*Sdxy_px;
    Tdxy_py = PBy.*Tdxy_s + 1.*oo2p.*(Tpx_s) + 2.*abop.*Sdxy_py;
    Tdxy_pz = PBz.*Tdxy_s + 2.*abop.*Sdxy_pz;

    Tdxz_px = PBx.*Tdxz_s + 1.*oo2p.*(Tpz_s) + 2.*abop.*Sdxz_px;
    Tdxz_py = PBy.*Tdxz_s + 2.*abop.*Sdxz_py;
    Tdxz_pz = PBz.*Tdxz_s + 1.*oo2p.*(Tpx_s) + 2.*abop.*Sdxz_pz;

    Tdyy_px = PBx.*Tdyy_s + 2.*abop.*Sdyy_px;
    Tdyy_py = PBy.*Tdyy_s + 2.*oo2p.*(Tpy_s) + 2.*abop.*Sdyy_py;
    Tdyy_pz = PBz.*Tdyy_s + 2.*abop.*Sdyy_pz;

    Tdyz_px = PBx.*Tdyz_s + 2.*abop.*Sdyz_px;
    Tdyz_py = PBy.*Tdyz_s + 1.*oo2p.*(Tpz_s) + 2.*abop.*Sdyz_py;
    Tdyz_pz = PBz.*Tdyz_s + 1.*oo2p.*(Tpy_s) + 2.*abop.*Sdyz_pz;

    Tdzz_px = PBx.*Tdzz_s + 2.*abop.*Sdzz_px;
    Tdzz_py = PBy.*Tdzz_s + 2.*abop.*Sdzz_py;
    Tdzz_pz = PBz.*Tdzz_s + 2.*oo2p.*(Tpz_s) + 2.*abop.*Sdzz_pz;

    Sab = zeros(6,3);
    Sab(1,1) = sum(Sdxx_px);
    Sab(1,2) = sum(Sdxx_py);
    Sab(1,3) = sum(Sdxx_pz);

    Sab(2,1) = sum(Sdxy_px);
    Sab(2,2) = sum(Sdxy_py);
    Sab(2,3) = sum(Sdxy_pz);

    Sab(3,1) = sum(Sdxz_px);
    Sab(3,2) = sum(Sdxz_py);
    Sab(3,3) = sum(Sdxz_pz);

    Sab(4,1) = sum(Sdyy_px);
    Sab(4,2) = sum(Sdyy_py);
    Sab(4,3) = sum(Sdyy_pz);

    Sab(5,1) = sum(Sdyz_px);
    Sab(5,2) = sum(Sdyz_py);
    Sab(5,3) = sum(Sdyz_pz);

    Sab(6,1) = sum(Sdzz_px);
    Sab(6,2) = sum(Sdzz_py);
    Sab(6,3) = sum(Sdzz_pz);

    Tab = zeros(6,3);
    Tab(1,1) = sum(Tdxx_px);
    Tab(1,2) = sum(Tdxx_py);
    Tab(1,3) = sum(Tdxx_pz);

    Tab(2,1) = sum(Tdxy_px);
    Tab(2,2) = sum(Tdxy_py);
    Tab(2,3) = sum(Tdxy_pz);

    Tab(3,1) = sum(Tdxz_px);
    Tab(3,2) = sum(Tdxz_py);
    Tab(3,3) = sum(Tdxz_pz);

    Tab(4,1) = sum(Tdyy_px);
    Tab(4,2) = sum(Tdyy_py);
    Tab(4,3) = sum(Tdyy_pz);

    Tab(5,1) = sum(Tdyz_px);
    Tab(5,2) = sum(Tdyz_py);
    Tab(5,3) = sum(Tdyz_pz);

    Tab(6,1) = sum(Tdzz_px);
    Tab(6,2) = sum(Tdzz_py);
    Tab(6,3) = sum(Tdzz_pz);

elseif (L1 == 2 && L2 == 2) %[d|d]

    Ss_s = WeightValues.*KabValues.*(pi./pValues).^1.5;
    abop = aValues.*bValues./pValues;
    oo2p = 0.5./pValues;
    Ts_s= Ss_s.*(abop.*(3-2*abop.*R12sq));

    Spx_s = PAx.*Ss_s;
    Spy_s = PAy.*Ss_s;
    Spz_s = PAz.*Ss_s;

    Tpx_s = PAx.*Ts_s+2*abop.*(Spx_s);
    Tpy_s = PAy.*Ts_s+2*abop.*(Spy_s);
    Tpz_s = PAz.*Ts_s+2*abop.*(Spz_s);
    % %Ss_p
    %     Ss_px = PBx.*Ss_s;
    %     Ss_py = PBy.*Ss_s;
    %     Ss_pz = PBz.*Ss_s;
    % %Ts_p
    %     Ts_px = PBx.*Ts_s+2*abop.*(Ss_px);
    %     Ts_py = PBy.*Ts_s+2*abop.*(Ss_py);
    %     Ts_pz = PBz.*Ts_s+2*abop.*(Ss_pz);
%Sp_p
    Spx_px = PBx.*Spx_s+1*oo2p.*Ss_s;
    Spx_py = PBy.*Spx_s;
    Spx_pz = PBz.*Spx_s;

    Spy_px = PBx.*Spy_s;
    Spy_py = PBy.*Spy_s+1*oo2p.*Ss_s;
    Spy_pz = PBz.*Spy_s;

    Spz_px = PBx.*Spz_s;
    Spz_py = PBy.*Spz_s;
    Spz_pz = PBz.*Spz_s+1*oo2p.*Ss_s;
%Tp_p
    Tpx_px = PBx.*Tpx_s+1*oo2p.*Ts_s+2*abop.*(Spx_px);
    Tpx_py = PBy.*Tpx_s             +2*abop.*(Spx_py);
    Tpx_pz = PBz.*Tpx_s             +2*abop.*(Spx_pz);

    Tpy_px = PBx.*Tpy_s             +2*abop.*(Spy_px);
    Tpy_py = PBy.*Tpy_s+1*oo2p.*Ts_s+2*abop.*(Spy_py);
    Tpy_pz = PBz.*Tpy_s             +2*abop.*(Spy_pz);

    Tpz_px = PBx.*Tpz_s             +2*abop.*(Spz_px);
    Tpz_py = PBy.*Tpz_s             +2*abop.*(Spz_py);
    Tpz_pz = PBz.*Tpz_s+1*oo2p.*Ts_s+2*abop.*(Spz_pz);
%Sd_s
    Sdxx_s = PAx.*Spx_s+1*oo2p.*Ss_s;
    Sdxy_s = PAy.*Spx_s;
    Sdxz_s = PAz.*Spx_s;
    Sdyy_s = PAy.*Spy_s+1*oo2p.*Ss_s;
    Sdyz_s = PAz.*Spy_s;
    Sdzz_s = PAz.*Spz_s+1*oo2p.*Ss_s;
%Td_s
    Tdxx_s = PAx.*Tpx_s+2*abop.*(Sdxx_s+1*oo2p.*Ss_s)+1*oo2p.*Ts_s;
    Tdxy_s = PAy.*Tpx_s+2*abop.*(Sdxy_s);
    Tdxz_s = PAz.*Tpx_s+2*abop.*(Sdxz_s);
    Tdyy_s = PAy.*Tpy_s+2*abop.*(Sdyy_s+1*oo2p.*Ss_s)+1*oo2p.*Ts_s;
    Tdyz_s = PAz.*Tpy_s+2*abop.*(Sdyz_s);
    Tdzz_s = PAz.*Tpz_s+2*abop.*(Sdzz_s+1*oo2p.*Ss_s)+1*oo2p.*Ts_s;

%Sd_p
    Sdxx_px = PBx.*Sdxx_s + 2.*oo2p.*(Spx_s);
    Sdxx_py = PBy.*Sdxx_s;
    Sdxx_pz = PBz.*Sdxx_s;

    Sdxy_px = PBx.*Sdxy_s + 1.*oo2p.*(Spy_s);
    Sdxy_py = PBy.*Sdxy_s + 1.*oo2p.*(Spx_s);
    Sdxy_pz = PBz.*Sdxy_s;

    Sdxz_px = PBx.*Sdxz_s + 1.*oo2p.*(Spz_s);
    Sdxz_py = PBy.*Sdxz_s;
    Sdxz_pz = PBz.*Sdxz_s + 1.*oo2p.*(Spx_s);

    Sdyy_px = PBx.*Sdyy_s;
    Sdyy_py = PBy.*Sdyy_s + 2.*oo2p.*(Spy_s);
    Sdyy_pz = PBz.*Sdyy_s;

    Sdyz_px = PBx.*Sdyz_s;
    Sdyz_py = PBy.*Sdyz_s + 1.*oo2p.*(Spz_s);
    Sdyz_pz = PBz.*Sdyz_s + 1.*oo2p.*(Spy_s);

    Sdzz_px = PBx.*Sdzz_s;
    Sdzz_py = PBy.*Sdzz_s;
    Sdzz_pz = PBz.*Sdzz_s + 2.*oo2p.*(Spz_s);
%Td_p
    Tdxx_px = PBx.*Tdxx_s + 2.*oo2p.*(Tpx_s) + 2.*abop.*Sdxx_px;
    Tdxx_py = PBy.*Tdxx_s + 2.*abop.*Sdxx_py;
    Tdxx_pz = PBz.*Tdxx_s + 2.*abop.*Sdxx_pz;

    Tdxy_px = PBx.*Tdxy_s + 1.*oo2p.*(Tpy_s) + 2.*abop.*Sdxy_px;
    Tdxy_py = PBy.*Tdxy_s + 1.*oo2p.*(Tpx_s) + 2.*abop.*Sdxy_py;
    Tdxy_pz = PBz.*Tdxy_s + 2.*abop.*Sdxy_pz;

    Tdxz_px = PBx.*Tdxz_s + 1.*oo2p.*(Tpz_s) + 2.*abop.*Sdxz_px;
    Tdxz_py = PBy.*Tdxz_s + 2.*abop.*Sdxz_py;
    Tdxz_pz = PBz.*Tdxz_s + 1.*oo2p.*(Tpx_s) + 2.*abop.*Sdxz_pz;

    Tdyy_px = PBx.*Tdyy_s + 2.*abop.*Sdyy_px;
    Tdyy_py = PBy.*Tdyy_s + 2.*oo2p.*(Tpy_s) + 2.*abop.*Sdyy_py;
    Tdyy_pz = PBz.*Tdyy_s + 2.*abop.*Sdyy_pz;

    Tdyz_px = PBx.*Tdyz_s + 2.*abop.*Sdyz_px;
    Tdyz_py = PBy.*Tdyz_s + 1.*oo2p.*(Tpz_s) + 2.*abop.*Sdyz_py;
    Tdyz_pz = PBz.*Tdyz_s + 1.*oo2p.*(Tpy_s) + 2.*abop.*Sdyz_pz;

    Tdzz_px = PBx.*Tdzz_s + 2.*abop.*Sdzz_px;
    Tdzz_py = PBy.*Tdzz_s + 2.*abop.*Sdzz_py;
    Tdzz_pz = PBz.*Tdzz_s + 2.*oo2p.*(Tpz_s) + 2.*abop.*Sdzz_pz;

%Sd_d
    Sdxx_dxx = PBx.*Sdxx_px + 1*oo2p.*(Sdxx_s) + 2*oo2p.*(Spx_px);
    Sdxx_dxy = PBx.*Sdxx_py + 2*oo2p.*(Spx_py);
    Sdxx_dxz = PBx.*Sdxx_pz + 2*oo2p.*(Spx_pz);
    Sdxx_dyy = PBy.*Sdxx_py + 1*oo2p.*(Sdxx_s);
    Sdxx_dyz = PBy.*Sdxx_pz;
    Sdxx_dzz = PBz.*Sdxx_pz + 1*oo2p.*(Sdxx_s);
    Sdxy_dxx = PBx.*Sdxy_px + 1*oo2p.*(Sdxy_s) + 1*oo2p.*(Spy_px);
    Sdxy_dxy = PBx.*Sdxy_py + 1*oo2p.*(Spy_py);
    Sdxy_dxz = PBx.*Sdxy_pz + 1*oo2p.*(Spy_pz);
    Sdxy_dyy = PBy.*Sdxy_py + 1*oo2p.*(Sdxy_s) + 1*oo2p.*(Spx_py);
    Sdxy_dyz = PBy.*Sdxy_pz + 1*oo2p.*(Spx_pz);
    Sdxy_dzz = PBz.*Sdxy_pz + 1*oo2p.*(Sdxy_s);
    Sdxz_dxx = PBx.*Sdxz_px + 1*oo2p.*(Sdxz_s) + 1*oo2p.*(Spz_px);
    Sdxz_dxy = PBx.*Sdxz_py + 1*oo2p.*(Spz_py);
    Sdxz_dxz = PBx.*Sdxz_pz + 1*oo2p.*(Spz_pz);
    Sdxz_dyy = PBy.*Sdxz_py + 1*oo2p.*(Sdxz_s);
    Sdxz_dyz = PBy.*Sdxz_pz;
    Sdxz_dzz = PBz.*Sdxz_pz + 1*oo2p.*(Sdxz_s) + 1*oo2p.*(Spx_pz);
    Sdyy_dxx = PBx.*Sdyy_px + 1*oo2p.*(Sdyy_s);
    Sdyy_dxy = PBx.*Sdyy_py;
    Sdyy_dxz = PBx.*Sdyy_pz;
    Sdyy_dyy = PBy.*Sdyy_py + 1*oo2p.*(Sdyy_s) + 2*oo2p.*(Spy_py);
    Sdyy_dyz = PBy.*Sdyy_pz + 2*oo2p.*(Spy_pz);
    Sdyy_dzz = PBz.*Sdyy_pz + 1*oo2p.*(Sdyy_s);
    Sdyz_dxx = PBx.*Sdyz_px + 1*oo2p.*(Sdyz_s);
    Sdyz_dxy = PBx.*Sdyz_py;
    Sdyz_dxz = PBx.*Sdyz_pz;
    Sdyz_dyy = PBy.*Sdyz_py + 1*oo2p.*(Sdyz_s) + 1*oo2p.*(Spz_py);
    Sdyz_dyz = PBy.*Sdyz_pz + 1*oo2p.*(Spz_pz);
    Sdyz_dzz = PBz.*Sdyz_pz + 1*oo2p.*(Sdyz_s) + 1*oo2p.*(Spy_pz);
    Sdzz_dxx = PBx.*Sdzz_px + 1*oo2p.*(Sdzz_s);
    Sdzz_dxy = PBx.*Sdzz_py;
    Sdzz_dxz = PBx.*Sdzz_pz;
    Sdzz_dyy = PBy.*Sdzz_py + 1*oo2p.*(Sdzz_s);
    Sdzz_dyz = PBy.*Sdzz_pz;
    Sdzz_dzz = PBz.*Sdzz_pz + 1*oo2p.*(Sdzz_s) + 2*oo2p.*(Spz_pz);

%Td_d
    Tdxx_dxx = PBx.*Tdxx_px + 1*oo2p.*(Tdxx_s) + 2*oo2p.*(Tpx_px) + 2*abop.*Sdxx_dxx -2*abop.*1*oo2p.*(Sdxx_s);
    Tdxx_dxy = PBx.*Tdxx_py + 2*oo2p.*(Tpx_py) + 2*abop.*Sdxx_dxy;
    Tdxx_dxz = PBx.*Tdxx_pz + 2*oo2p.*(Tpx_pz) + 2*abop.*Sdxx_dxz;
    Tdxx_dyy = PBy.*Tdxx_py + 1*oo2p.*(Tdxx_s) + 2*abop.*Sdxx_dyy -2*abop.*1*oo2p.*(Sdxx_s);
    Tdxx_dyz = PBy.*Tdxx_pz + 2*abop.*Sdxx_dyz;
    Tdxx_dzz = PBz.*Tdxx_pz + 1*oo2p.*(Tdxx_s) + 2*abop.*Sdxx_dzz -2*abop.*1*oo2p.*(Sdxx_s);
    Tdxy_dxx = PBx.*Tdxy_px + 1*oo2p.*(Tdxy_s) + 1*oo2p.*(Tpy_px) + 2*abop.*Sdxy_dxx -2*abop.*1*oo2p.*(Sdxy_s);
    Tdxy_dxy = PBx.*Tdxy_py + 1*oo2p.*(Tpy_py) + 2*abop.*Sdxy_dxy;
    Tdxy_dxz = PBx.*Tdxy_pz + 1*oo2p.*(Tpy_pz) + 2*abop.*Sdxy_dxz;
    Tdxy_dyy = PBy.*Tdxy_py + 1*oo2p.*(Tdxy_s) + 1*oo2p.*(Tpx_py) + 2*abop.*Sdxy_dyy -2*abop.*1*oo2p.*(Sdxy_s);
    Tdxy_dyz = PBy.*Tdxy_pz + 1*oo2p.*(Tpx_pz) + 2*abop.*Sdxy_dyz;
    Tdxy_dzz = PBz.*Tdxy_pz + 1*oo2p.*(Tdxy_s) + 2*abop.*Sdxy_dzz -2*abop.*1*oo2p.*(Sdxy_s);
    Tdxz_dxx = PBx.*Tdxz_px + 1*oo2p.*(Tdxz_s) + 1*oo2p.*(Tpz_px) + 2*abop.*Sdxz_dxx -2*abop.*1*oo2p.*(Sdxz_s);
    Tdxz_dxy = PBx.*Tdxz_py + 1*oo2p.*(Tpz_py) + 2*abop.*Sdxz_dxy;
    Tdxz_dxz = PBx.*Tdxz_pz + 1*oo2p.*(Tpz_pz) + 2*abop.*Sdxz_dxz;
    Tdxz_dyy = PBy.*Tdxz_py + 1*oo2p.*(Tdxz_s) + 2*abop.*Sdxz_dyy -2*abop.*1*oo2p.*(Sdxz_s);
    Tdxz_dyz = PBy.*Tdxz_pz + 2*abop.*Sdxz_dyz;
    Tdxz_dzz = PBz.*Tdxz_pz + 1*oo2p.*(Tdxz_s) + 1*oo2p.*(Tpx_pz) + 2*abop.*Sdxz_dzz -2*abop.*1*oo2p.*(Sdxz_s);
    Tdyy_dxx = PBx.*Tdyy_px + 1*oo2p.*(Tdyy_s) + 2*abop.*Sdyy_dxx -2*abop.*1*oo2p.*(Sdyy_s);
    Tdyy_dxy = PBx.*Tdyy_py + 2*abop.*Sdyy_dxy;
    Tdyy_dxz = PBx.*Tdyy_pz + 2*abop.*Sdyy_dxz;
    Tdyy_dyy = PBy.*Tdyy_py + 1*oo2p.*(Tdyy_s) + 2*oo2p.*(Tpy_py) + 2*abop.*Sdyy_dyy -2*abop.*1*oo2p.*(Sdyy_s);
    Tdyy_dyz = PBy.*Tdyy_pz + 2*oo2p.*(Tpy_pz) + 2*abop.*Sdyy_dyz;
    Tdyy_dzz = PBz.*Tdyy_pz + 1*oo2p.*(Tdyy_s) + 2*abop.*Sdyy_dzz -2*abop.*1*oo2p.*(Sdyy_s);
    Tdyz_dxx = PBx.*Tdyz_px + 1*oo2p.*(Tdyz_s) + 2*abop.*Sdyz_dxx -2*abop.*1*oo2p.*(Sdyz_s);
    Tdyz_dxy = PBx.*Tdyz_py + 2*abop.*Sdyz_dxy;
    Tdyz_dxz = PBx.*Tdyz_pz + 2*abop.*Sdyz_dxz;
    Tdyz_dyy = PBy.*Tdyz_py + 1*oo2p.*(Tdyz_s) + 1*oo2p.*(Tpz_py) + 2*abop.*Sdyz_dyy -2*abop.*1*oo2p.*(Sdyz_s);
    Tdyz_dyz = PBy.*Tdyz_pz + 1*oo2p.*(Tpz_pz) + 2*abop.*Sdyz_dyz;
    Tdyz_dzz = PBz.*Tdyz_pz + 1*oo2p.*(Tdyz_s) + 1*oo2p.*(Tpy_pz) + 2*abop.*Sdyz_dzz -2*abop.*1*oo2p.*(Sdyz_s);
    Tdzz_dxx = PBx.*Tdzz_px + 1*oo2p.*(Tdzz_s) + 2*abop.*Sdzz_dxx -2*abop.*1*oo2p.*(Sdzz_s);
    Tdzz_dxy = PBx.*Tdzz_py + 2*abop.*Sdzz_dxy;
    Tdzz_dxz = PBx.*Tdzz_pz + 2*abop.*Sdzz_dxz;
    Tdzz_dyy = PBy.*Tdzz_py + 1*oo2p.*(Tdzz_s) + 2*abop.*Sdzz_dyy -2*abop.*1*oo2p.*(Sdzz_s);
    Tdzz_dyz = PBy.*Tdzz_pz + 2*abop.*Sdzz_dyz;
    Tdzz_dzz = PBz.*Tdzz_pz + 1*oo2p.*(Tdzz_s) + 2*oo2p.*(Tpz_pz) + 2*abop.*Sdzz_dzz -2*abop.*1*oo2p.*(Sdzz_s);

    S = zeros(6,6);

    S(1,1) = sum(Sdxx_dxx);
    S(1,2) = sum(Sdxx_dxy);
    S(1,3) = sum(Sdxx_dxz);
    S(1,4) = sum(Sdxx_dyy);
    S(1,5) = sum(Sdxx_dyz);
    S(1,6) = sum(Sdxx_dzz);
    S(2,1) = sum(Sdxy_dxx);
    S(2,2) = sum(Sdxy_dxy);
    S(2,3) = sum(Sdxy_dxz);
    S(2,4) = sum(Sdxy_dyy);
    S(2,5) = sum(Sdxy_dyz);
    S(2,6) = sum(Sdxy_dzz);
    S(3,1) = sum(Sdxz_dxx);
    S(3,2) = sum(Sdxz_dxy);
    S(3,3) = sum(Sdxz_dxz);
    S(3,4) = sum(Sdxz_dyy);
    S(3,5) = sum(Sdxz_dyz);
    S(3,6) = sum(Sdxz_dzz);
    S(4,1) = sum(Sdyy_dxx);
    S(4,2) = sum(Sdyy_dxy);
    S(4,3) = sum(Sdyy_dxz);
    S(4,4) = sum(Sdyy_dyy);
    S(4,5) = sum(Sdyy_dyz);
    S(4,6) = sum(Sdyy_dzz);
    S(5,1) = sum(Sdyz_dxx);
    S(5,2) = sum(Sdyz_dxy);
    S(5,3) = sum(Sdyz_dxz);
    S(5,4) = sum(Sdyz_dyy);
    S(5,5) = sum(Sdyz_dyz);
    S(5,6) = sum(Sdyz_dzz);
    S(6,1) = sum(Sdzz_dxx);
    S(6,2) = sum(Sdzz_dxy);
    S(6,3) = sum(Sdzz_dxz);
    S(6,4) = sum(Sdzz_dyy);
    S(6,5) = sum(Sdzz_dyz);
    S(6,6) = sum(Sdzz_dzz);

    T = zeros(6,6);

    T(1,1) = sum(Tdxx_dxx);
    T(1,2) = sum(Tdxx_dxy);
    T(1,3) = sum(Tdxx_dxz);
    T(1,4) = sum(Tdxx_dyy);
    T(1,5) = sum(Tdxx_dyz);
    T(1,6) = sum(Tdxx_dzz);
    T(2,1) = sum(Tdxy_dxx);
    T(2,2) = sum(Tdxy_dxy);
    T(2,3) = sum(Tdxy_dxz);
    T(2,4) = sum(Tdxy_dyy);
    T(2,5) = sum(Tdxy_dyz);
    T(2,6) = sum(Tdxy_dzz);
    T(3,1) = sum(Tdxz_dxx);
    T(3,2) = sum(Tdxz_dxy);
    T(3,3) = sum(Tdxz_dxz);
    T(3,4) = sum(Tdxz_dyy);
    T(3,5) = sum(Tdxz_dyz);
    T(3,6) = sum(Tdxz_dzz);
    T(4,1) = sum(Tdyy_dxx);
    T(4,2) = sum(Tdyy_dxy);
    T(4,3) = sum(Tdyy_dxz);
    T(4,4) = sum(Tdyy_dyy);
    T(4,5) = sum(Tdyy_dyz);
    T(4,6) = sum(Tdyy_dzz);
    T(5,1) = sum(Tdyz_dxx);
    T(5,2) = sum(Tdyz_dxy);
    T(5,3) = sum(Tdyz_dxz);
    T(5,4) = sum(Tdyz_dyy);
    T(5,5) = sum(Tdyz_dyz);
    T(5,6) = sum(Tdyz_dzz);
    T(6,1) = sum(Tdzz_dxx);
    T(6,2) = sum(Tdzz_dxy);
    T(6,3) = sum(Tdzz_dxz);
    T(6,4) = sum(Tdzz_dyy);
    T(6,5) = sum(Tdzz_dyz);
    T(6,6) = sum(Tdzz_dzz);
    Sab = S;
    Tab = T;
end

end
