function gabcd = shellOS(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d)

% Nints = basis_a.n*basis_b.n*basis_c.n*basis_d.n;
% boysValues = zeros(Nints,(L1+L2+L3+L4+1));
% gSSSSNValues = zeros(Nints,(L1+L2+L3+L4+1));
% xValues = zeros(Nints,1);
% indexValues = zeros(Nints,1);
% xIndexValues = zeros(Nints,1);
% DxValues = zeros(Nints,1);
% KabValues = zeros(Nints,1);
% KcdValues = zeros(Nints,1);
% RPAValues = zeros(Nints,3);
% RPBValues = zeros(Nints,3);
% RQCValues = zeros(Nints,3);
% RQDValues = zeros(Nints,3);
% RWPValues = zeros(Nints,3);
% RWQValues = zeros(Nints,3);
% pValues = zeros(Nints,1);
% qValues = zeros(Nints,1);
% alphaValues = zeros(Nints,1);
% RPQ2Values = zeros(Nints,1);
% PrefactorValues = zeros(Nints,1);
% WeightValues = zeros(Nints,1);
%
% xstep = 0.1; %Nlast/Npoints

if ( L1 == 0 && L2 == 0 && (L3 > 0 || L4 > 0))

gswap = shellOS(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table,pair_data2,c,d,a,b);
gabcd = permute(gswap,[3 4 1 2]);

elseif ( L2 > L1)

gswap = shellOS(basis_b,basis_a,basis_c,basis_d,L2,L1,L3,L4,Boys_Table,pair_data2,b,a,c,d);
gabcd = permute(gswap,[2 1 3 4]);

elseif ( L4 > L3)

gswap = shellOS(basis_a,basis_b,basis_c,basis_d,L1,L2,L4,L3,Boys_Table,pair_data2,a,b,d,c);
gabcd = permute(gswap,[1 2 4 3]);

elseif ( L1 == 1 && L2 == 0 && L3 == 0 && L4 == 0) %px_s_s_s_0

%[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors2(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table);
[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gpsss = OSpsss(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
gabcd = gpsss;

elseif ( L1 == 0 && L2 == 1 && L3 == 0 && L4 == 0) %s_px_s_s_0

%[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors2(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table);
[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gspss = OSspss(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
gabcd = gspss;

elseif ( L1 == 0 && L2 == 0 && L3 == 1 && L4 == 0) %s_s_p_s_0

gssps_swap = shellOS(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table,pair_data2,c,d,a,b);
gabcd = permute(gssps_swap,[3 4 1 2]);

elseif ( L1 == 0 && L2 == 0 && L3 == 0 && L4 == 1) %s_s_p_s_0

gsssp_swap = shellOS(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table,pair_data2,c,d,a,b);
gabcd = permute(gsssp_swap,[3 4 1 2]);

elseif ( L1 == 1 && L2 == 0 && L3 == 1 && L4 == 0) %px_s_px_s_0

%[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors2(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table);
[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);

gpsps = OSpsps(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
gabcd = gpsps;

elseif ( L1 == 1 && L2 == 0 && L3 == 1 && L4 == 1)

    gpspp_swap = shellOS(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table,pair_data2,c,d,a,b);
    gabcd = permute(gpspp_swap,[3 4 1 2]);

elseif ( L1 == 0 && L2 == 1 && L3 == 1 && L4 == 1)

    gsppp_swap = shellOS(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table,pair_data2,c,d,a,b);
    gabcd = permute(gsppp_swap,[3 4 1 2]);

elseif( L1 == 1 && L2 == 1 && L3 == 0 && L4 == 0)

%[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors2(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table);
[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gppss = OSppss(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
gabcd = gppss;

elseif( L1 == 0 && L2 == 0 && L3 == 1 && L4 == 1)

%[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors2(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table);
[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gsspp = OSsspp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
gabcd = gsspp;

elseif( L1 == 0 && L2 == 1 && L3 == 1 && L4 == 0)

%[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors2(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table);
[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gspps = OSspps(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
gabcd = gspps;

elseif( L1 == 0 && L2 == 1 && L3 == 0 && L4 == 1)

%[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors2(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table);
[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gspsp = OSspsp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
gabcd = gspsp;

elseif( L1 == 1 && L2 == 0 && L3 == 0 && L4 == 1)

%[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors2(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table);
[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gpssp = OSpssp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
gabcd = gpssp;

elseif( L1 == 1 && L2 == 1 && L3 == 1 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gppps = OSppps(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
gabcd = gppps;

elseif( L1 == 1 && L2 == 1 && L3 == 0 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gppsp = OSppsp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
gabcd = gppsp;

elseif( L1 == 1 && L2 == 1 && L3 == 1 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gpppp = OSpppp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
gabcd = gpppp;






%dsXX (9 integrals)--------------------------------------------------------------------------

elseif( L1 == 2 && L2 == 0 && L3 == 0 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdsss(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 0 && L3 == 0 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdssp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 0 && L3 == 0 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdssd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 0 && L3 == 1 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdsps(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 0 && L3 == 1 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdspp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 0 && L3 == 1 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdspd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 0 && L3 == 2 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdsds(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 0 && L3 == 2 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdsdp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 0 && L3 == 2 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdsdd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

%dpXX (9 integrals)----------------------------------------------------------------------

elseif( L1 == 2 && L2 == 1 && L3 == 0 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdpss(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 1 && L3 == 0 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdpsp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 1 && L3 == 0 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdpsd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 1 && L3 == 1 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdpps(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 1 && L3 == 1 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdppp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 1 && L3 == 1 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdppd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 1 && L3 == 2 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdpds(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 1 && L3 == 2 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdpdp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 1 && L3 == 2 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdpdd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

%sdXX (9 integrals)--------------------------------------------------------------------------------------------------
elseif( L1 == 0 && L2 == 2 && L3 == 0 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSsdss(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 0 && L2 == 2 && L3 == 0 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSsdsp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 0 && L2 == 2 && L3 == 0 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSsdsd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 0 && L2 == 2 && L3 == 1 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSpdps(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 0 && L2 == 2 && L3 == 1 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSsdpp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 0 && L2 == 2 && L3 == 1 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSsdpd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 0 && L2 == 2 && L3 == 2 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSsdds(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 0 && L2 == 2 && L3 == 2 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSsddp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 0 && L2 == 2 && L3 == 2 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSsddd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

%pdXX (9 integrals)--------------------------------------------------------------------------------------------------
elseif( L1 == 1 && L2 == 2 && L3 == 0 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSpdss(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 1 && L2 == 2 && L3 == 0 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSpdsp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 1 && L2 == 2 && L3 == 0 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSpdsd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 1 && L2 == 2 && L3 == 1 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSpdps(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 1 && L2 == 2 && L3 == 1 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSpdpp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 1 && L2 == 2 && L3 == 1 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSpdpd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 1 && L2 == 2 && L3 == 2 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSpdds(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 1 && L2 == 2 && L3 == 2 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSpddp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 1 && L2 == 2 && L3 == 2 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSpddd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

%ddXX (9 integrals)---------------------------------------------------------------------------------------------------------------------
elseif( L1 == 2 && L2 == 2 && L3 == 0 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSddss(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 2 && L3 == 0 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSddsp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 2 && L3 == 0 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSddsd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 2 && L3 == 1 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSddps(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 2 && L3 == 1 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSddpp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 2 && L3 == 1 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSddpd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 2 && L3 == 2 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSddds(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 2 && L3 == 2 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdddp(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

elseif( L1 == 2 && L2 == 2 && L3 == 2 && L4 == 2)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d);
gabcd = OSdddd(RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);

else

gabcd_swap = shellOS(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table,pair_data2,c,d,a,b);
gabcd = permute(gabcd_swap,[3 4 1 2]);

end  %if

end
