function [RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues,boysTableTime,timeReadPairData] = primitiveFactors3(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d)

%Nints = basis_a.n*basis_b.n*basis_c.n*basis_d.n;
%boysValues = zeros(Nints,(L1+L2+L3+L4+1));
%gSSSSNValues = zeros(Nints,(L1+L2+L3+L4+1));
%xValues = zeros(Nints,1);
%indexValues = zeros(Nints,1);
%xIndexValues = zeros(Nints,1);
%DxValues = zeros(Nints,1);
%KabValues = zeros(Nints,1);
%KcdValues = zeros(Nints,1);
%RPAValues = zeros(Nints,3);
%RPBValues = zeros(Nints,3);
%RQCValues = zeros(Nints,3);
%RQDValues = zeros(Nints,3);
%RWPValues = zeros(Nints,3);
%RWQValues = zeros(Nints,3);
%pValues = zeros(Nints,1);
%qValues = zeros(Nints,1);
%alphaValues = zeros(Nints,1);
%RPQ2Values = zeros(Nints,1);
%PrefactorValues = zeros(Nints,1);
%WeightValues = zeros(Nints,1);


%Nints = Nintsab*Nintscd;

%gSSSSNValues = zeros(Nints,1);
%xValues = zeros(Nints,1);
%indexValues = zeros(Nints,1);
%xIndexValues = zeros(Nints,1);
%DxValues = zeros(Nints,1);
%{
KabValues = zeros(Nintsab,1);
KcdValues = zeros(Nintscd,1);

pValues = zeros(Nintsab,1);
qValues = zeros(Nintscd,1);

WeightValuesab = zeros(Nintsab,1);
WeightValuescd = zeros(Nintscd,1);

PxValues = zeros(Nintsab,1);
PyValues = zeros(Nintsab,1);
PzValues = zeros(Nintsab,1);

QxValues = zeros(Nintscd,1);
QyValues = zeros(Nintscd,1);
QzValues = zeros(Nintscd,1);

RPAValues = zeros(Nintsab,3);
RPBValues = zeros(Nintsab,3);
RQCValues = zeros(Nintscd,3);
RQDValues = zeros(Nintscd,3);
%}
%RWPValues = zeros(Nints,3);
%RWQValues = zeros(Nints,3);


%alphaValues = zeros(Nints,1);
%RPQ2Values = zeros(Nints,1);
%PrefactorValues = zeros(Nints,1);



%RAB = [basis_a.g(1).x0-basis_b.g(1).x0;basis_a.g(1).y0-basis_b.g(1).y0;basis_a.g(1).z0-basis_b.g(1).z0];
%RCD = [basis_c.g(1).x0-basis_d.g(1).x0;basis_c.g(1).y0-basis_d.g(1).y0;basis_c.g(1).z0-basis_d.g(1).z0];

xstep = 0.1; %Nlast/Npoints
%t = 1;
%{
for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
    g1 = basis_a.g(na);
    aa = g1.alpha;
    c1 = basis_a.c(na);
    N1 = g1.N;
    c1N1 = c1*N1;
    for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
        g2 = basis_b.g(nb);
        ab = g2.alpha;
        c2 = basis_b.c(nb);
        N2 = g2.N;
        c1N1c2N2 = c1N1*c2*N2;

        p = aa + ab;
        Px = (aa*g1.x0 + ab*g2.x0)/p;
        Py = (aa*g1.y0 + ab*g2.y0)/p;
        Pz = (aa*g1.z0 + ab*g2.z0)/p;

        RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
        RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
        RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector

        rhoAB = aa*ab/p;
        Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
        pValues(t) = p;
        PxValues(t) = Px;
        PyValues(t) = Py;
        PzValues(t) = Pz;
        RPAValues(t,:) = [RPA(1),RPA(2),RPA(3)];
        RPBValues(t,:) = [RPB(1),RPB(2),RPB(3)];
        KabValues(t) = Kab;
        WeightValuesab(t) = c1N1*c2*N2;
        t = t + 1;
    end
end

s = 1;
for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function
        g3 = basis_c.g(nc);
        ac = g3.alpha;
        c3 = basis_c.c(nc);
        N3 = g3.N;
        c3N3 = c3*N3;

        for nd = 1:basis_d.n; %loops over number of primitives in 4th contracted basis functions.
            g4 = basis_d.g(nd);
            ad = g4.alpha;
            c4 = basis_d.c(nd);
            N4 = g4.N;

            q = ac+ad;
            Qx = (ac*g3.x0 + ad*g4.x0)/q;
            Qy = (ac*g3.y0 + ad*g4.y0)/q;
            Qz = (ac*g3.z0 + ad*g4.z0)/q;

            RQC = [Qx-g3.x0;Qy-g3.y0;Qz-g3.z0];
            RQD = [Qx-g4.x0;Qy-g4.y0;Qz-g4.z0];
            RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];

            rhoCD = ac*ad/q;
            Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));

            %RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
            %RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
            %alpha = q*p/(q+p);

            %x = alpha*RPQ2;

            qValues(s) = q;
            QxValues(s) = Qx;
            QyValues(s) = Qy;
            QzValues(s) = Qz;

            RQCValues(s,:) = [RQC(1),RQC(2),RQC(3)];
            RQDValues(s,:) = [RQD(1),RQD(2),RQD(3)];
            KcdValues(s) = Kcd;
            WeightValuescd(s) = c3N3*c4*N4;

            s = s + 1;

        end
end
%}
%{
disp('pValues calculated in-situ');
disp(pValues);
disp('PxValues calculated in-situ');
disp(PxValues);
disp('RPAValues calculated in-situ');
disp(RPAValues);
disp('qValues calculated in-situ');
disp(qValues);
disp('QxValues calculated in-situ');
disp(QxValues);
disp('RQCValues calculated in-situ');
disp(RQCValues);
%}
%pair_data2 is a cell structure with the following fields.
%pair_data2{a,b} = [pValues PxValues PyValues PzValues KabValues WeightValuesab RPAValues RPBValues]
prevtime = time();
pair_data2_ab = pair_data2{a,b};
pair_data2_cd = pair_data2{c,d};

pValues  = pair_data2_ab(:,1);
PxValues = pair_data2_ab(:,2);
PyValues = pair_data2_ab(:,3);
PzValues = pair_data2_ab(:,4);
KabValues = pair_data2_ab(:,5);
WeightValuesab = pair_data2_ab(:,6);
RPAValues = pair_data2_ab(:,7:9);
RPBValues = pair_data2_ab(:,10:12);
pPxValues = pair_data2_ab(:,13);
pPyValues = pair_data2_ab(:,14);
pPzValues = pair_data2_ab(:,15);

qValues  = pair_data2_cd(:,1);
QxValues = pair_data2_cd(:,2);
QyValues = pair_data2_cd(:,3);
QzValues = pair_data2_cd(:,4);
KcdValues = pair_data2_cd(:,5);
WeightValuescd = pair_data2_cd(:,6);
RQCValues = pair_data2_cd(:,7:9);
RQDValues = pair_data2_cd(:,10:12);
qQxValues = pair_data2_cd(:,13);
qQyValues = pair_data2_cd(:,14);
qQzValues = pair_data2_cd(:,15);

timeReadPairData = time()-prevtime;
%{
disp('pValues from pair_data2');
disp(pValues);
disp('PxValues from pair_data2');
disp(PxValues);
disp('RPAValues from pair_data2');
disp(RPAValues);
disp('qValues from pair_data2');
disp(qValues);
disp('QxValues from pair_data2');
disp(QxValues);
disp('RQCValues from pair_data2');
disp(RQCValues);
%}

%The size of RPAValues, etc was (Nintsab,3), it has to be expanded to be (Nints,3)
Nintsab = length(pValues);
Nintscd = length(qValues);

onesab = ones(Nintsab,1);
onescd = ones(Nintscd,1);

RPAValues = kron(RPAValues,onescd); %Have to check that these Kronecker products give me the correct ordering
RPBValues = kron(RPBValues,onescd);
RQCValues = kron(onesab,RQCValues);
RQDValues = kron(onesab,RQDValues);

%These are the values that depend on quartet quantities
RPQx = PxValues'-QxValues; %' %This will give a matrix
RPQy = PyValues'-QyValues; %' %This will give a matrix
RPQz = PzValues'-QzValues; %' %This will give a matrix

RPQ2 = RPQx.*RPQx+RPQy.*RPQy+RPQz.*RPQz;

ppqValues = qValues + pValues'; %' This should be a matrix
%disp(ppqValues);
pqValues = qValues*pValues'; %'
alpha = pqValues./ppqValues;

x = alpha.*RPQ2;

%For now, size(pValues) = Nintsab;
% WxValues = ((qValues.*QxValues)+(pValues.*PxValues)')./ppqValues; %(Ncd,Nab) matrix '
% WyValues = ((qValues.*QyValues)+(pValues.*PyValues)')./ppqValues; %(Ncd,Nab) matrix '
% WzValues = ((qValues.*QzValues)+(pValues.*PzValues)')./ppqValues; %(Ncd,Nab) matrix '

WxValues = ((qQxValues)+(pPxValues)')./ppqValues; %(Ncd,Nab) matrix '
WyValues = ((qQyValues)+(pPyValues)')./ppqValues; %(Ncd,Nab) matrix '
WzValues = ((qQzValues)+(pPzValues)')./ppqValues; %(Ncd,Nab) matrix '

%Now that I have used pValues and qValues, I can expand them to the full size.
pValues = kron(pValues,onescd);
qValues = kron(onesab,qValues);

%disp('WxValues');
%disp(WxValues);
%disp('ones(Nintscd,1)*PxValues');
%disp(ones(Nintscd,1)*PxValues'); %'

RWPValues_x = WxValues-onescd*PxValues'; %'
RWPValues_y = WyValues-onescd*PyValues'; %'
RWPValues_z = WzValues-onescd*PzValues'; %'

onesabT = onesab';
RWQValues_x = WxValues-QxValues*onesabT; %
RWQValues_y = WyValues-QyValues*onesabT; %
RWQValues_z = WzValues-QzValues*onesabT; %


RWPValues = [RWPValues_x(:),RWPValues_y(:),RWPValues_z(:)]; %Now this should be a (Nab*Ncd,3) matrix
RWQValues = [RWQValues_x(:),RWQValues_y(:),RWQValues_z(:)];

%disp('RWPValues_x');
%disp(RWPValues_x);
%disp('RWPValues');
%disp(RWPValues);


KabcdValues = KcdValues*KabValues'; %'

PrefactorValues = KabcdValues.*2*17.493418327624862./(pqValues.*sqrt(ppqValues));
PrefactorValues = PrefactorValues(:);

WeightValues = WeightValuescd*WeightValuesab'; %'
WeightValues = WeightValues(:);

xValues = x(:);
ppqValues = ppqValues(:);

below35 = (xValues < 35);
xValuesLow = xValues(below35);
xValuesHigh = xValues(~below35);

indexValues = floor(xValuesLow/xstep)+1;
xIndexValues = (indexValues-1)*xstep;
DxValues = (xValuesLow-xIndexValues); %Difference which enters the Taylor expansion

Dx2Values = DxValues.*DxValues;
Dx3Values = Dx2Values.*DxValues;
Dx4Values = Dx2Values.*Dx2Values;
Dx5Values = Dx3Values.*Dx2Values;

Dx2_o_2_Values = .5*Dx2Values;
Dx3_o_6_Values = .166666666666667*Dx3Values;
Dx4_o_24_Values = 4.166666666666666e-2*Dx4Values;
Dx5_o_120_Values = 8.333333333333334e-3*Dx5Values;

%ppqValues = pValues + qValues;
%PrefactorValues = KabValues.*KcdValues.*2*17.493418327624862./(pValues.*qValues.*sqrt(ppqValues));

prevtime = time();
% boysValuesLow = zeros(length(xValuesLow),L1+L2+L3+L4+1);
% for order = 0:(L1+L2+L3+L4)
%     boysValuesLow(:,order+1) = Boys_Table(indexValues,order+1)-Boys_Table(indexValues,order+2).*DxValues+...
%                                Boys_Table(indexValues,order+3).*Dx2_o_2_Values-Boys_Table(indexValues,order+4).*Dx3_o_6_Values+...
%                                Boys_Table(indexValues,order+5).*Dx4_o_24_Values-Boys_Table(indexValues,order+6).*Dx5_o_120_Values;
% end

MaxOrder = L1+L2+L3+L4;
switch (MaxOrder)
case (0)
  boysValuesLow = Boys_Table(indexValues,1)-Boys_Table(indexValues,2).*DxValues+...
                             Boys_Table(indexValues,3).*Dx2_o_2_Values-Boys_Table(indexValues,4).*Dx3_o_6_Values+...
                             Boys_Table(indexValues,5).*Dx4_o_24_Values-Boys_Table(indexValues,6).*Dx5_o_120_Values;
case (1)
    bt1 = Boys_Table(indexValues,1);
    bt2 = Boys_Table(indexValues,2);
    bt3 = Boys_Table(indexValues,3);
    bt4 = Boys_Table(indexValues,4);
    bt5 = Boys_Table(indexValues,5);
    bt6 = Boys_Table(indexValues,6);
    bt7 = Boys_Table(indexValues,7);
    boysValuesLow0 = bt1-bt2.*DxValues+...
                             bt3.*Dx2_o_2_Values-bt4.*Dx3_o_6_Values+...
                             bt5.*Dx4_o_24_Values-bt6.*Dx5_o_120_Values;
    boysValuesLow1 = bt2-bt3.*DxValues+...
                              bt4.*Dx2_o_2_Values-bt5.*Dx3_o_6_Values+...
                              bt6.*Dx4_o_24_Values-bt7.*Dx5_o_120_Values;
    boysValuesLow = [boysValuesLow0,boysValuesLow1];
case (2)
      bt1 = Boys_Table(indexValues,1);
      bt2 = Boys_Table(indexValues,2);
      bt3 = Boys_Table(indexValues,3);
      bt4 = Boys_Table(indexValues,4);
      bt5 = Boys_Table(indexValues,5);
      bt6 = Boys_Table(indexValues,6);
      bt7 = Boys_Table(indexValues,7);
      bt8 = Boys_Table(indexValues,8);
      boysValuesLow0 = bt1-bt2.*DxValues+...
                         bt3.*Dx2_o_2_Values-bt4.*Dx3_o_6_Values+...
                         bt5.*Dx4_o_24_Values-bt6.*Dx5_o_120_Values;
      boysValuesLow1 = bt2-bt3.*DxValues+...
                          bt4.*Dx2_o_2_Values-bt5.*Dx3_o_6_Values+...
                          bt6.*Dx4_o_24_Values-bt7.*Dx5_o_120_Values;
      boysValuesLow2 = bt3-bt4.*DxValues+...
                          bt5.*Dx2_o_2_Values-bt6.*Dx3_o_6_Values+...
                          bt7.*Dx4_o_24_Values-bt8.*Dx5_o_120_Values;
      boysValuesLow = [boysValuesLow0,boysValuesLow1,boysValuesLow2];
case (3)
      bt1 = Boys_Table(indexValues,1);
      bt2 = Boys_Table(indexValues,2);
      bt3 = Boys_Table(indexValues,3);
      bt4 = Boys_Table(indexValues,4);
      bt5 = Boys_Table(indexValues,5);
      bt6 = Boys_Table(indexValues,6);
      bt7 = Boys_Table(indexValues,7);
      bt8 = Boys_Table(indexValues,8);
      bt9 = Boys_Table(indexValues,9);
      boysValuesLow0 = bt1-bt2.*DxValues+...
                         bt3.*Dx2_o_2_Values-bt4.*Dx3_o_6_Values+...
                         bt5.*Dx4_o_24_Values-bt6.*Dx5_o_120_Values;
      boysValuesLow1 = bt2-bt3.*DxValues+...
                          bt4.*Dx2_o_2_Values-bt5.*Dx3_o_6_Values+...
                          bt6.*Dx4_o_24_Values-bt7.*Dx5_o_120_Values;
      boysValuesLow2 = bt3-bt4.*DxValues+...
                          bt5.*Dx2_o_2_Values-bt6.*Dx3_o_6_Values+...
                          bt7.*Dx4_o_24_Values-bt8.*Dx5_o_120_Values;
      boysValuesLow3 = bt4-bt5.*DxValues+...
                          bt6.*Dx2_o_2_Values-bt7.*Dx3_o_6_Values+...
                          bt8.*Dx4_o_24_Values-bt9.*Dx5_o_120_Values;
      boysValuesLow = [boysValuesLow0,boysValuesLow1,boysValuesLow2,boysValuesLow3];
case (4)

      boysValuesLow0 = Boys_Table(indexValues,1)-Boys_Table(indexValues,2).*DxValues+...
                               Boys_Table(indexValues,3).*Dx2_o_2_Values-Boys_Table(indexValues,4).*Dx3_o_6_Values+...
                               Boys_Table(indexValues,5).*Dx4_o_24_Values-Boys_Table(indexValues,6).*Dx5_o_120_Values;
      boysValuesLow1 = Boys_Table(indexValues,2)-Boys_Table(indexValues,3).*DxValues+...
                                Boys_Table(indexValues,4).*Dx2_o_2_Values-Boys_Table(indexValues,5).*Dx3_o_6_Values+...
                                Boys_Table(indexValues,6).*Dx4_o_24_Values-Boys_Table(indexValues,7).*Dx5_o_120_Values;
      boysValuesLow2 = Boys_Table(indexValues,3)-Boys_Table(indexValues,4).*DxValues+...
                                Boys_Table(indexValues,5).*Dx2_o_2_Values-Boys_Table(indexValues,6).*Dx3_o_6_Values+...
                                Boys_Table(indexValues,7).*Dx4_o_24_Values-Boys_Table(indexValues,8).*Dx5_o_120_Values;
      boysValuesLow3 = Boys_Table(indexValues,4)-Boys_Table(indexValues,5).*DxValues+...
                                Boys_Table(indexValues,6).*Dx2_o_2_Values-Boys_Table(indexValues,7).*Dx3_o_6_Values+...
                                Boys_Table(indexValues,8).*Dx4_o_24_Values-Boys_Table(indexValues,9).*Dx5_o_120_Values;
      boysValuesLow4 = Boys_Table(indexValues,5)-Boys_Table(indexValues,6).*DxValues+...
                                Boys_Table(indexValues,7).*Dx2_o_2_Values-Boys_Table(indexValues,8).*Dx3_o_6_Values+...
                                Boys_Table(indexValues,9).*Dx4_o_24_Values-Boys_Table(indexValues,10).*Dx5_o_120_Values;
      boysValuesLow = [boysValuesLow0,boysValuesLow1,boysValuesLow2,boysValuesLow3,boysValuesLow4];
otherwise
  boysValuesLow = zeros(length(xValuesLow),L1+L2+L3+L4+1);
  for order = 0:(L1+L2+L3+L4)
      boysValuesLow(:,order+1) = Boys_Table(indexValues,order+1)-Boys_Table(indexValues,order+2).*DxValues+...
                                 Boys_Table(indexValues,order+3).*Dx2_o_2_Values-Boys_Table(indexValues,order+4).*Dx3_o_6_Values+...
                                 Boys_Table(indexValues,order+5).*Dx4_o_24_Values-Boys_Table(indexValues,order+6).*Dx5_o_120_Values;
  endfor
endswitch


boysTableTime = time()-prevtime;
if ~isempty(xValuesHigh)
    sqrt_pi = 1.772453850905516;
    boysValuesHigh_0 = 0.5*sqrt_pi./sqrt(xValuesHigh); %this calculates the order 0 boys function for x > 35

    %I precompute the quantity (prod(n+1:2*n))/2^(2*n+1)*sqrt_pi
    Prefactor = [443.113462726379e-003
                  664.670194089569e-003
                  1.66167548522392e+000
                  5.81586419828372e+000
                  26.1713888922768e+000
                  143.942638907522e+000
                  935.627152898894e+000
                  7.01720364674171e+003
                  59.6462309973045e+003
                  566.639194474393e+003
                  5.94971154198112e+006
                  68.4216827327829e+006
                  855.271034159787e+006
                  11.5461589611571e+009
                  167.419304936778e+009
                  2.59499922652006e+012
                  42.8174872375810e+012
                  749.306026657668e+012
                  13.8621614931669e+015
                  270.312149116754e+015];
    % boysValuesHigh_1_N = zeros(length(xValuesHigh),L1+L2+L3+L4);

    % for order = 1:(L1+L2+L3+L4)
    %     boysValuesHigh_1_N(:,order) = Prefactor(order)./xValuesHigh.^(order+0.5);
    % end

    % disp(boysValuesHigh_1_N);
    orders = 1:(L1+L2+L3+L4);
    boysValuesHigh_1_N = Prefactor(orders)'./(xValuesHigh.^(orders+0.5));
    % disp('Now using broadcasting')
    % disp(boysValuesHigh_1_N);
    boysValuesHigh_0_N = [boysValuesHigh_0,boysValuesHigh_1_N];

    boysValues(below35==1,:) = boysValuesLow;
    boysValues(below35==0,:) = boysValuesHigh_0_N;
    %disp(boysValues);

    %boysValues = [boysValues;boysValuesHigh_0_N];
else
    boysValues = boysValuesLow;
end

% for order = 0:(L1+L2+L3+L4)
%     gSSSSNValues(:,order+1) = WeightValues.*PrefactorValues.*boysValues(:,order+1);
% end
% disp(gSSSSNValues);
gSSSSNValues = (WeightValues.*PrefactorValues).*boysValues; %I'm using broadcasting
% disp('Now using broadcasting')
%disp(gSSSSNValues);
%Outputs: RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,ppqValues,gSSSSNValues

end
