clear, close all
tic
%Boys_Table 400,4000,order 0...10

Boys_Table = loadBoysTable2;

%basissetStr = 'def2-svp';
basissetStr = 'STO-3G';
Z = [6 1 1 1 1];
Charge = 0;
spin = 0;
AL = [-0.000000000000   0.000000000000   0.000000000000
       1.183772681898  -1.183771681898  -1.183771681898
       1.183771681898   1.183772681898   1.183771681898
      -1.183771681898   1.183771681898  -1.183771681898
      -1.183771681898  -1.183771681898   1.183771681898];

unrestricted = 1;

%This function creates a cell array of column vectors. Each vector contains the contracted shell information
[basis2, atomBasisList,Nel] = Build_Basis_2(Z,AL,basissetStr);
Nel = Nel-Charge;

%This function creates a matrix listing the contracted basis function indices
[Shell_List,Sph_Shell_List] = Build_Shell_List(basis2);

%This are thresholds for truncating the following pair_data2 matrices.
primitiveThr = 1e-14; %Threshold for discarding primitives
overlapThr = 1e-14; %Threshold for combining primitives

%This function calculates pair_data2, a 2-D cell array of matrices which allows to retrieve
%information on shell pairs, necessary to compute 1- and 2- electron integrals
[pair_data2,RABdata,pOverlap] = basis_products_2(basis2,primitiveThr,overlapThr);

%This function has been programmed in may 2019. Seems to work so far
[S,T,Ssph,Tsph] = Build_OS1e(basis2,pair_data2,Shell_List,Sph_Shell_List);

%Function for calculating electron-nuclear attraction Ven
VenOctave = Build_OS_EN(basis2,pair_data2,Shell_List,AL,Z,Boys_Table);
[Ven,VenSph,tInitVectorsAccu,tLoopVectorsAccu] = Build_OS_EN_time(basis2,pair_data2,Shell_List,Sph_Shell_List,AL,Z,Boys_Table);

%Now I need to generate the Hcore matrix
Hcore0 = T+Ven;

gabcd = 0;
if (rem (Nel,2) == 0) && (unrestricted == 0)
[Min_Energy,E,ncycle,D,Dinit,Coulomb,Exchange,epsilon,F,Fprime,G] = rSCF(Hcore0,T,Ven,S,Nel,basis2,pair_data2,Shell_List,Boys_Table,gabcd);
%[Min_Energy,E,ncycle,D,Dinit,Coulomb,Exchange,epsilon,F,Fprime,G] = rSCF(HcoreRef,T,Ven,Sref,Nel,basis2,pair_data2,Shell_List,Boys_Table,gabcd);
  Da = D/2;
  Db = Da;
end

if (rem (Nel, 2) == 1) || (unrestricted == 1)
  if rem(Nel,2) == 0
    Nela = Nel/2;
    Nelb = Nel/2;
  else
    Nela = (Nel+2*spin)/2;
    Nelb = Nel-Nela;
  end
[Min_Energy,E,ncycle,Da,Db,Dinita,Dinitb,Coulomb,Exchangea,Exchangeb,epsilona,epsilonb,Fa,Fprimea,Fb,Fprimeb] = uSCF(Hcore0,T,Ven,S,Nela,Nelb,basis2,pair_data2,Shell_List,Boys_Table,gabcd,Z,atomBasisList);
  D = Da + Db;
end


TotalEnergy = E + Nuclear_Repulsion(Z,AL);

[MullikenCharge,MullikenSpin] = mullikenAnalysis(S,Da,Db,Z,atomBasisList,Shell_List);

%{
basis = basis2;
a = 8;
b = 3;
c = 11;
d = 10;

tmp = primitiveFactorsSSSS_3(basis{a},basis{b},basis{c},basis{d},Boys_Table,pair_data2,a,b,c,d);
[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors3(basis{a},basis{b},basis{c},basis{d},basis{a}.L,basis{b}.L,basis{c}.L,basis{d}.L,Boys_Table,pair_data2,a,b,c,d);

figure
xx = [basis{a}.g(1).x0;basis{b}.g(1).x0;basis{c}.g(1).x0;basis{d}.g(1).x0];
yy = [basis{a}.g(1).y0;basis{b}.g(1).y0;basis{c}.g(1).y0;basis{d}.g(1).y0];
zz = [basis{a}.g(1).z0;basis{b}.g(1).z0;basis{c}.g(1).z0;basis{d}.g(1).z0];

xx = [xx; RPAValues(:,1); RWPValues(:,1)];
yy = [yy; RPAValues(:,2); RWPValues(:,2)];
zz = [zz; RPAValues(:,3); RWPValues(:,3)];

size_nuc = 400;
size_RPA = 100;
size_RWP = 100;
Sizes = [size_nuc*ones(4,1);size_RPA*ones(size(RPAValues(:,1)));size_RWP*ones(size(RWPValues(:,1)))];

Pvalues_x = RPAValues(:,1)+basis{a}.g(1).x0;
Pvalues_y = RPAValues(:,2)+basis{a}.g(1).y0;
Pvalues_z = RPAValues(:,3)+basis{a}.g(1).z0;

Qvalues_x = RQCValues(:,1)+basis{c}.g(1).x0;
Qvalues_y = RQCValues(:,2)+basis{c}.g(1).y0;
Qvalues_z = RQCValues(:,3)+basis{c}.g(1).z0;

Wvalues_x = RWPValues(:,1)+Pvalues_x;
Wvalues_y = RWPValues(:,2)+Pvalues_y;
Wvalues_z = RWPValues(:,3)+Pvalues_z;

%scatter3(xx,yy,zz,Sizes);
scatter3(basis{a}.g(1).x0,basis{a}.g(1).y0,basis{a}.g(1).z0,600*ones(size(basis{a}.g(1).z0)),'ok','MarkerFaceColor','k'); hold on
scatter3(basis{b}.g(1).x0,basis{b}.g(1).y0,basis{b}.g(1).z0,600*ones(size(basis{a}.g(1).z0)),'ok','MarkerFaceColor','k'); hold on
scatter3(basis{c}.g(1).x0,basis{c}.g(1).y0,basis{c}.g(1).z0,600*ones(size(basis{a}.g(1).z0)),'ok','MarkerFaceColor','k'); hold on
scatter3(basis{d}.g(1).x0,basis{d}.g(1).y0,basis{d}.g(1).z0,600*ones(size(basis{a}.g(1).z0)),'ok','MarkerFaceColor','k'); hold on
%scatter3(RPAValues(:,1),RPAValues(:,2),RPAValues(:,3),200*ones(size(RPAValues(:,3))),'or','MarkerFaceColor','r'); hold on
%scatter3(RWPValues(:,1),RWPValues(:,2),RWPValues(:,3),200*ones(size(RPAValues(:,3))),'ok','MarkerFaceColor','k');
scatter3(Pvalues_x,Pvalues_y,Pvalues_z,300*ones(size(Pvalues_x)),'or','MarkerFaceColor','r'); hold on
scatter3(Qvalues_x,Qvalues_y,Qvalues_z,300*ones(size(Qvalues_x)),'ob','MarkerFaceColor','b'); hold on
scatter3(Wvalues_x,Wvalues_y,Wvalues_z,200*ones(size(Wvalues_x)),'om','MarkerFaceColor','m');
axis([-3 3 -3 3 -5 2]);
figure
plot(pValues,'.k','MarkerFaceColor','k');


%}
