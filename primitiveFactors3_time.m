function [RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues,boysTableTime,timeReadPairData] = primitiveFactors3_time(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,pair_data2,a,b,c,d)

prevtime = time();
pair_data2_ab = pair_data2{a,b};
pair_data2_cd = pair_data2{c,d};

timeReadPairData = time()-prevtime;

%Call oct-function
[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,ooppqValues,xValues,PrefactorValues,WeightValues] = loopPrimitives(pair_data2_ab,pair_data2_cd);

%Now the Boys function algorithm starts
xstep = 0.1;
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
