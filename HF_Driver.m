clear, close all
tic
%Boys_Table 400,4000,order 0...10

Boys_Table = loadBoysTable2;

nz = cell(1,12);
nz{2} = [1;0;0;1;0;1];
nz{3} = [2;1;1;0;0;0;2;1;0;2];
nz{4} = [3;2;2;1;1;1;0;0;0;0;3;2;1;0;3];
nz{5} = [4;3;3;2;2;2;1;1;1;1;0;0;0;0;0;4;3;2;1;0;4];
nz{6} = [5;4;4;3;3;3;2;2;2;2;1;1;1;1;1;0;0;0;0;0;0;5;4;3;2;1;0;5];
nz{7} = [6;5;5;4;4;4;3;3;3;3;2;2;2;2;2;1;1;1;1;1;1;0;0;0;0;0;0;0;6;5;4;3;2;1;0;6];
nz{8} = [7;6;6;5;5;5;4;4;4;4;3;3;3;3;3;2;2;2;2;2;2;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;7;6;5;4;3;2;1;0;7];
nz{9} = [8;7;7;6;6;6;5;5;5;5;4;4;4;4;4;3;3;3;3;3;3;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;8;7;6;5;4;3;2;1;0;8];
nz{10} = [9;8;8;7;7;7;6;6;6;6;5;5;5;5;5;4;4;4;4;4;4;3;3;3;3;3;3;3;2;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;9;8;7;6;5;4;3;2;1;0;9];
nz{11} = [10;9;9;8;8;8;7;7;7;7;6;6;6;6;6;5;5;5;5;5;5;4;4;4;4;4;4;4;3;3;3;3;3;3;3;3;2;2;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;10;9;8;7;6;5;4;3;2;1;0;10];
nz{12} = [11;10;10;9;9;9;8;8;8;8;7;7;7;7;7;6;6;6;6;6;6;5;5;5;5;5;5;5;4;4;4;4;4;4;4;4;3;3;3;3;3;3;3;3;3;2;2;2;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;11;10;9;8;7;6;5;4;3;2;1;0;11];

unique = [2, 4,5, 7,8,9, 11,12,13,14, 16,17,18,19,20, 22,23,24,25,26,27,  29,30,31,32,33,34,35, 37,38,39,40,41,42,43,44, 46,47,48,49,50,51,52,53,54, 56,57,58,59,60,61,62,63,64,65, 67,68,69,70,71,72,73,74,75,76,77, 79,80,81,82,83,84,85,86,87,88,89,90]; %up to [ms|ss]



%   Z = [2 6 1 1 1 1]; %H, He
%  Z = [Z Z];
% %Z = [Z Z Z];
%
%acetilene
% Z = [6 6 1 1];
%% AL = [0.05 0 -0.5
%%       -0.05 0 0.5
%%       0.05 0 -1.4
%%       -0.05 0 1.4]*1.8897;
%Charge = 0;
%AL = [0.094486    0.000000   -0.944863
%-0.094486    0.000000    0.944863
% 0.094486    0.000000   -2.645617
%-0.094486    0.000000    2.645617];



%H2O
% Z = [8 1 1];
% AL = [0.00 0.000 0.116
%       0.000 0.751 -0.465
%       0.000 -0.751 -0.465]*1.8897;
%CO
% Z = [6 8];
% r = 1.10*1.8897;
% theta = 30*pi/180;
% phi = 45*pi/180;
% AL = [0.00 0.00 0
%       sin(theta)*cos(phi)*r sin(theta)*sin(phi)*r cos(theta)*r];

%Methane
 Z = [6 1 1 1 1];
 Charge = 0;
 spin = 0;
 AL = [-0.000000000000   0.000000000000   0.000000000000
      1.183772681898  -1.183771681898  -1.183771681898
      1.183771681898   1.183772681898   1.183771681898
      -1.183771681898   1.183771681898  -1.183771681898
      -1.183771681898  -1.183771681898   1.183771681898];

% Z = [6 1 1];
% Charge = 0;
% spin = 0;
% AL = [-0.000000000000   0.000000000000   0.000000000000
%      1.183772681898  -1.183771681898  -1.183771681898
%      1.183771681898   1.183772681898   1.183771681898];
% Z = [6];
% Charge = 0;
% spin = 0;
% AL = [-0.000000000000   0.000000000000   0.000000000000];
unrestricted = 1;
 % Z = [Z Z];
 % AL = [AL;AL+3*ones(size(AL))];

%Fictituous Molecule with
%      Z = [6 1 1 1 1];
%      Charge = 0;
%      AL = [-0.000000000000   0.000000000000   0.000000000000
%           2.183772681898  -2.183771681898  -1.143771681898
%           1.183771681898   3.183772681898   1.169771681898
%           -1.444771681898   1.883771681898  -1.03771681898
%           -1.123771681898  -1.783771681898   1.483771681898];
% %
% % %1 Bohr = 0.5291772108 A;
% %
% AL = [-0.000000000000   0.000000000000   0.000000000000
%     1.183772  -1.183772  -1.183772
%     1.183772   1.183772   1.183772
%     -1.183772   1.183772  -1.183772
%     -1.183772  -1.183772   1.1837712];
% AL = AL*1.5;

% Z = [6];
% Charge = 0;
% AL = [-0.000000000000   0.000000000000   0.000000000000];

%Fraction of molecule_004000
% Z = [8 7 6];
% Charge = 0;
% AL = [-0.2849591644	1.7190324354	0.0525888862
% 0.1014649249	0.4895501827	0.5551454214
% -0.8014553823	-0.4157525006	0.401989349]*1.8897;

%Methane + 2 HeH+
% Z = [6 1 1 1 1 2 1 2 1];
% Charge = 2;
% AL = [-0.000000000000   0.000000000000   0.000000000000
%    1.183771681898  -1.183771681898  -1.183771681898
%    1.183771681898   1.183771681898   1.183771681898
%    -1.183771681898   1.183771681898  -1.183771681898
%    -1.183771681898  -1.183771681898   1.183771681898
%    0   0   5
%    0   0.3 5.8
%    0   0   -5
%    0   -0.3 -5.8];

%CHe2H2
% Z = [6 2 2 1 1];
% AL = [-0.000000000000   0.000000000000   0.000000000000
%     1.183771681898  -1.183771681898  -1.183771681898
%     1.183771681898   1.183771681898   1.183771681898
%     -1.183771681898   1.183771681898  -1.183771681898
%     -1.183771681898  -1.183771681898   1.183771681898];
%
%  Z = [Z Z];
%  AL = [AL;AL+4*ones(size(AL))];

% AL =	[-4 -4 -4
%         -0.000000000000   0.000000000000   0.000000000000
%          1.183771681898  -1.183771681898  -1.183771681898
%          1.183771681898   1.183771681898   1.183771681898
%          -1.183771681898   1.183771681898  -1.183771681898
%          -1.183771681898  -1.183771681898   1.183771681898];
% % AL = AL*1.8897;
%  AL = [AL;AL+3*ones(size(AL))];
%AL = [AL;AL+2*ones(size(AL));AL-2*ones(size(AL))];

%  Z = [1 1];
%  AL = [0 0 0;0 0 0.7408478]*1.8897;
%HeH+ hexamer
% Z = [2 1];
% Z = [Z Z Z Z Z Z];
% AL = [0 0 0;0 0 1.4632];
% AL2 = [0.21 0.5 0;0 0 1.3632];
% AL = [AL;AL+2*ones(size(AL));AL-[3 0 0;3 0 0];AL+4*ones(size(AL));AL2;AL2-5*ones(size(AL2))];

% Ndistances = 1;
% rmsd_gabcd_k = zeros(1,Ndistances);
% for ka = 1:Ndistances


%Commented (11/05/2019)
% [basis,Nel] = Build_Basis(Z,AL,'6-31G');
% Nel = Nel-Charge;

%This function creates a cell array of column vectors. Each vector contains the contracted shell information
[basis2, atomBasisList,Nel] = Build_Basis_2(Z,AL,'def2-svp');
Nel = Nel-Charge;

%This function creates a matrix listing the contracted basis function indices
%(for example the indices in the S, T, VNuc, and F matrices) and the corresponding
%contracted shell indices (remember that a shell can have matrix dimensions larger than 1)
[Shell_List,Sph_Shell_List] = Build_Shell_List(basis2);

%This are thresholds for truncating the following pair_data2 matrices.
%They are now set on very low values so they don't bother
primitiveThr = 1e-14; %Threshold for discarding primitives
overlapThr = 1e-14; %Threshold for combining primitives

%This function calculates pair_data2, a 2-D cell array of matrices which allows to retrieve
%information on shell pairs, necessary to compute 1- and 2- electron integrals
%This function is right now the slowest part, because the shell pair formation uses loops
%I will not optimize it right now because ERI integrals will take much more.
[pair_data2,RABdata,pOverlap] = basis_products_2(basis2,primitiveThr,overlapThr);

%This function has been programmed in may 2019. Seems to work so far
[S,T,Ssph,Tsph] = Build_OS1e(basis2,pair_data2,Shell_List,Sph_Shell_List);

% cart2sphMat = [sqrt(3)/2 0 0 -sqrt(3)/2 0 0
%               -1/2 0 0 -1/2 0 1
%                        0 1 0          0 0 0
%                        0 0 1          0 0 0
%                        0 0 0          0 1 0];
% cart2sphMat = [sqrt(3)/2 0 0 -sqrt(3)/2 0 0
%              -1/2 0 0 -1/2 0 1
%                       0 sqrt(3) 0          0 0 0
%                       0 0 sqrt(3)          0 0 0
%                       0 0 0          0 sqrt(3) 0];
%
% Ncont = Shell_List(end,2);
% cart2sphBigMat = blkdiag(eye(9),cart2sphMat',eye(Ncont-15));
%
% Ssph = cart2sphBigMat'*S*cart2sphBigMat;

Sref = [1.0000000   0.2483624   0.0000000  -0.0000000   0.0000000   0.0630068   0.0630068   0.0630068   0.0630068
       0.2483624   1.0000000   0.0000000   0.0000000   0.0000000   0.4936348   0.4936348   0.4936348   0.4936348
       0.0000000   0.0000000   1.0000000   0.0000000   0.0000000   0.2707690   0.2707690  -0.2707690  -0.2707690
      -0.0000000   0.0000000   0.0000000   1.0000000   0.0000000  -0.2707690   0.2707690   0.2707690  -0.2707690
       0.0000000   0.0000000   0.0000000   0.0000000   1.0000000  -0.2707690   0.2707690  -0.2707690   0.2707690
       0.0630068   0.4936348   0.2707690  -0.2707690  -0.2707690   1.0000000   0.1714170   0.1714170   0.1714170
       0.0630068   0.4936348   0.2707690   0.2707690   0.2707690   0.1714170   1.0000000   0.1714170   0.1714170
       0.0630068   0.4936348  -0.2707690   0.2707690  -0.2707690   0.1714170   0.1714170   1.0000000   0.1714170
       0.0630068   0.4936348  -0.2707690  -0.2707690   0.2707690   0.1714170   0.1714170   0.1714170   1.0000000];
% disp('S difference is');
% disp(S-Sref);


%Function for calculating electron-nuclear attraction Ven
VenOctave = Build_OS_EN(basis2,pair_data2,Shell_List,AL,Z,Boys_Table);
[Ven,VenSph,tInitVectorsAccu,tLoopVectorsAccu] = Build_OS_EN_time(basis2,pair_data2,Shell_List,Sph_Shell_List,AL,Z,Boys_Table);

if(~isequal(Ven,VenOctave))
  disp('Something is wrong with Ven');
end

%Now I need to generate the Hcore matrix
Hcore0 = T+Ven;

% HcoreRef = [-19.7125033  -4.7626451   0.0000000   0.0000000   0.0000000  -1.2180554  -1.2180554  -1.2180554  -1.2180554
%       -4.7626451  -6.6122791   0.0000000   0.0000000   0.0000000  -2.8730395  -2.8730395  -2.8730395  -2.8730395
%         0.0000000   0.0000000  -5.5649805   0.0000000   0.0000000  -1.3591225  -1.3591225   1.3591225   1.3591225
%         0.0000000   0.0000000   0.0000000  -5.5649805   0.0000000   1.3591225  -1.3591225  -1.3591225   1.3591225
%         0.0000000   0.0000000   0.0000000   0.0000000  -5.5649805   1.3591225  -1.3591225   1.3591225  -1.3591225
%        -1.2180554  -2.8730395  -1.3591225   1.3591225   1.3591225  -4.2234865  -0.9394502  -0.9394502  -0.9394502
%        -1.2180554  -2.8730395  -1.3591225  -1.3591225  -1.3591225  -0.9394502  -4.2234865  -0.9394502  -0.9394502
%        -1.2180554  -2.8730395   1.3591225  -1.3591225   1.3591225  -0.9394502  -0.9394502  -4.2234865  -0.9394502
%        -1.2180554  -2.8730395   1.3591225   1.3591225  -1.3591225  -0.9394502  -0.9394502  -0.9394502  -4.2234865];

% disp('Hcore difference is');
% disp(Hcore0-HcoreRef);

%gabcd = Build_ERI_OS_2(basis2,Shell_List,Boys_Table,pair_data2);

% gsymtest1 = squeeze(all((gabcd-permute(gabcd,[2,1,3,4]))<1e-8));
% gsymtest2 = squeeze(all((gabcd-permute(gabcd,[1,2,4,3]))<1e-8));
% gsymtest3 = squeeze(all((gabcd-permute(gabcd,[2,1,4,3]))<1e-8));
% gsymtest4 = squeeze(all((gabcd-permute(gabcd,[3,4,2,1]))<1e-8));
% gsymtest5 = squeeze(all((gabcd-permute(gabcd,[4,3,1,2]))<1e-8));
% gsymtest6 = squeeze(all((gabcd-permute(gabcd,[4,3,2,1]))<1e-8));
% gsymtest7 = squeeze(all((gabcd-permute(gabcd,[3,4,1,2]))<1e-8));
% disp('gabcd permutational symmetry tests');
% disp(squeeze(all(gsymtest1)));
% disp(squeeze(all(gsymtest2)));
% disp(squeeze(all(gsymtest3)));
% disp(squeeze(all(gsymtest4)));
% disp(squeeze(all(gsymtest5)));
% disp(squeeze(all(gsymtest6)));
% disp(squeeze(all(gsymtest7)));

% gabcd_sherrill = load('eri_CH4_Sherrill.dat');
%
% diffERI = zeros(912,1);
% for t = 1:912
% diffERI(t) = (gabcd(gabcd_sherrill(t,1),gabcd_sherrill(t,2),gabcd_sherrill(t,3),gabcd_sherrill(t,4))-gabcd_sherrill(t,5));
% end
% figure; plot(diffERI);
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

%From D. Crawford data
% Dref = [1.0329323  -0.1010877  -0.0000000   0.0000000   0.0000000  -0.0469085  -0.0469085  -0.0469085  -0.0469085
%       -0.1010877   0.3960228   0.0000000  -0.0000000  -0.0000000   0.1131900   0.1131900   0.1131900   0.1131900
%       -0.0000000   0.0000000   0.3268651   0.0000000   0.0000000   0.1721029   0.1721029  -0.1721029  -0.1721029
%        0.0000000  -0.0000000   0.0000000   0.3268651  -0.0000000  -0.1721029   0.1721029   0.1721029  -0.1721029
%        0.0000000  -0.0000000   0.0000000  -0.0000000   0.3268651  -0.1721029   0.1721029  -0.1721029   0.1721029
%       -0.0469085   0.1131900   0.1721029  -0.1721029  -0.1721029   0.3045239  -0.0579428  -0.0579428  -0.0579428
%       -0.0469085   0.1131900   0.1721029   0.1721029   0.1721029  -0.0579428   0.3045239  -0.0579428  -0.0579428
%       -0.0469085   0.1131900  -0.1721029   0.1721029  -0.1721029  -0.0579428  -0.0579428   0.3045239  -0.0579428
%       -0.0469085   0.1131900  -0.1721029  -0.1721029   0.1721029  -0.0579428  -0.0579428  -0.0579428   0.3045239];
%
% Fref = [-8.9530154  -2.5417508  -0.0000000  -0.0000000  -0.0000000  -0.6617374  -0.6617374  -0.6617374  -0.6617374
%        -2.5417508  -0.8174464  -0.0000000  -0.0000000   0.0000000  -0.3883845  -0.3883845  -0.3883845  -0.3883845
%        -0.0000000  -0.0000000   0.4192075   0.0000000   0.0000000  -0.0541874  -0.0541874   0.0541874   0.0541874
%        -0.0000000  -0.0000000   0.0000000   0.4192075  -0.0000000   0.0541874  -0.0541874  -0.0541874   0.0541874
%        -0.0000000   0.0000000   0.0000000  -0.0000000   0.4192075   0.0541874  -0.0541874   0.0541874  -0.0541874
%        -0.6617374  -0.3883845  -0.0541874   0.0541874   0.0541874   0.0420011  -0.1549906  -0.1549906  -0.1549906
%        -0.6617374  -0.3883845  -0.0541874  -0.0541874  -0.0541874  -0.1549906   0.0420011  -0.1549906  -0.1549906
%        -0.6617374  -0.3883845   0.0541874  -0.0541874   0.0541874  -0.1549906  -0.1549906   0.0420011  -0.1549906
%        -0.6617374  -0.3883845   0.0541874   0.0541874  -0.0541874  -0.1549906  -0.1549906  -0.1549906   0.0420011];
%
% DensityMatrixDeviation = sum(sum(abs(D/2-Dref)))/length(D)^2;
% FockMatrixDeviation = sum(sum(abs(F-Fref)))/length(F)^2;
% FockMatrixDifference = F-Fref;
%
% disp('DensityMatrixDeviation is');
% disp(DensityMatrixDeviation);
% disp('FockMatrixDeviation is');
% disp(FockMatrixDeviation);
% disp('FockMatrixDifference is');
% disp(FockMatrixDifference);

TotalEnergy = E + Nuclear_Repulsion(Z,AL);

[MullikenCharge,MullikenSpin] = mullikenAnalysis(S,Da,Db,Z,atomBasisList,Shell_List);


%The function for building the one electron matrices S and T has been changed in may 2019 to use the OS1e function
%{
[S,T,Pcell,pcell,S00cell] = Build_One_Electron_3(basis,Shell_Doublets,NShell_Doublets);

%Nuclear_Attraction = Build_Nuclear_Attraction_2(basis,Shell_Doublets,NShell_Doublets,AL,Z,Boys_Table);
Nuclear_Attraction = Build_Nuclear_Attraction_3(basis,Shell_List,AL,Z,Boys_Table);

%[Shells,NShells] = Build_Shells(basis);
%gabcd = Build_ERI(basis,Shells,NShells,Boys_Table,pair_data);
%gabcd = Build_ERI_4(basis,Shell_List,Boys_Table,pair_data);
gabcd = Build_ERI_OS(basis,Shell_List,Boys_Table,pair_data2);
shellThr = 1e-7;
PairThr = 1e-8;
[gabcd2,UniquePairData] = Build_Flat_ERI(basis,Shell_List,Boys_Table,pair_data2,pOverlap,shellThr,PairThr);
gabcd3 = Build_Flat_ERI_2(basis,Shell_List,Boys_Table,pair_data2,pOverlap,shellThr,PairThr);
%[gabcd_fast,int_error] = Build_ERI_fast(basis,Shell_List,Boys_Table,pair_data,S,T,Pcell,pcell,S00cell,gabcd);

% t = 0;
% Ncont = Shell_List(end,1);
% gabcd_list = zeros(nb^4,5);
% for i = 1:Ncont
%     for j = 1:Ncont
%         for k = 1:Ncont
%             for l = 1:Ncont
%                 t = t + 1;
%                 gabcd_list(t,:) = [i j k l gabcd(i,j,k,l)];
%             end
%         end
%     end
% end
%gabcd_list = reshape(gabcd,[],1);
%gabcd_fast_list = reshape(gabcd_fast,[],1);
%int_error_list = gabcd_list-gabcd_fast_list;
%rmsd_gabcd = sqrt(sum((gabcd_list-gabcd_fast_list).^2))/length(gabcd_list);

%This calculates the rmsd of the integrals whose error lies below a certain
%threshold.
%This is because some of the integrals have to big errors that dominate the
%rmsd calculation.
%int_error_list_below_threshold = int_error_list(abs(int_error_list)<0.01);
%rmsd_below_threshold = sqrt(sum(int_error_list_below_threshold.^2))/length(int_error_list_below_threshold);
%rmsd_gabcd_k(ka) = rmsd_below_threshold;

%
% if rmsd_gabcd < 0.05
%     rmsd_gabcd_k(k) = rmsd_gabcd;
% else
%     rmsd_gabcd_k(k) = rmsd_gabcd_k(k-1)
% end
%gabcd = Build_Electron_Repulsion(basis);

%This calculates the number of integrals which are below a certain
%threshold
%threshold = 0.002;
%Int_below_threshold = sum(abs(int_error_list)<threshold);
%Percent_below_threshold = Int_below_threshold/length(int_error_list)*100;

H0 = T+Nuclear_Attraction;

[Min_Energy,E,ncycle,D,Dinit,G,epsilon,F,Fprime] = SCF(H0,T,Nuclear_Attraction,gabcd,S,Nel,Shell_List);
%Min_Energy = Min_Energy + Nuclear_Repulsion(Z,AL);
E = E + Nuclear_Repulsion(Z,AL);
toc

if ka<4
 figure;
 plot((1:ncycle)',E,'.k','MarkerFaceColor','k'); %'
 xlabel('Cycle');
 ylabel('Energy (Eh)');

Ncont = Shell_List(end,2);
%figure;
%plot((1:Ncont^4)',gabcd_list,'.k',(1:Ncont^4)',gabcd_fast_list,'-r');% hold on
%plot((1:Ncont^4)',int_error_list,'-b'); %hold on %'
%plot((1:Ncont^4)',gabcd_fast_list./gabcd_list,'-m'); %'
%legend(num2str(rmsd_gabcd_k(ka)));
end

% end
% figure
% plot((1:Ndistances)/2,rmsd_gabcd_k,'sqk','MarkerFaceColor','k');
gabcd_sherrill = load('')
% diff = zeros(912,1);
% for t = 1:912
% diff(t) = (gabcd(gabcd_sherrill(t,1),gabcd_sherrill(t,2),gabcd_sherrill(t,3),gabcd_sherrill(t,4))-gabcd_sherrill(t,5));
% end
% figure; plot(diff);


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
