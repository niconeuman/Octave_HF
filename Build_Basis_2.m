function [basis, atomBasisList,Nel] = Build_Basis_2(Z,AL,basisset)

%Build_Basis2 created a cell array of vectors


Natoms = size(Z,2);
atomBasisList = zeros(Natoms,2);
Nel = 0;
nb = 0;
basis = cell(1000,1);

%TODO Complete the periodicTable
elementTable = {'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar'};



if (strcmp(basisset,'def2-svp') || strcmp(basisset,'def2-SVP') || strcmp(basisset,'Def2-SVP'))
    fileName = 'def2_svp.nw';
% if (strcmp(basisset,'mini'))
%     fileName = 'mini.nw';
    for index = 1:Natoms
        x0 = AL(index,1);
        y0 = AL(index,2);
        z0 = AL(index,3);
        Nel = Nel + Z(index);
        %nb = nb+1;
        atomBasisList(index,1) = nb+1; %Here I tell which basis set is the first for atom(index)

        elementName = elementTable{Z(index)};

        [coeffMatrices,matrixTypes] = readBasisSetFile(fileName,elementName);
        nMatrices = size(coeffMatrices,1);
        for k = 1:nMatrices
            nb = nb + 1;
            coeffMat = coeffMatrices{k};
            typeMat = matrixTypes{k};
            if strcmp(typeMat,'S')
                L = 0;
            elseif strcmp(typeMat,'P')
                L = 1;
            elseif strcmp(typeMat,'D')
                L = 2;
            elseif strcmp(typeMat,'SP')
                L = [0;1];
            end

            nPrim = size(coeffMat,1);

            if length(L) == 1
                %Data which is common to all primitives
                basis{nb}(1,1) = nPrim; %Number of primitives
                basis{nb}(2,1) = L; %Angular momentum
                basis{nb}(3) = x0;
                basis{nb}(4) = y0;
                basis{nb}(5) = z0;

                %Again, because normalization coefficients are different depending on the individual ml_x, ml_y ml_z of the
                %function with L > 1, here I choose the cartesian indices that give me the lowest normalization coefficient
                %for the angular momentum shell
                %The ratio of the other coefficients to this minimum value is built into the code-generated OS expressions
                if L == 0
                    indices = [0,0,0];
                elseif L == 1
                    indices = [1,0,0]; %here all are the same
                elseif L == 2
                    indices = [2,0,0]; %here [2,0,0] would have a different N (normalization) than [1,1,0]
                elseif L == 3
                    indices = [3,0,0];
                elseif L == 4
                    indices = [4,0,0];
                elseif L == 5
                    indices = [5,0,0];
                end

                % globalOverlap = 0;
                % globalOverlapPref = (pi^1.5)/(2^(L))*double_factorial(2*indices(1,1)-1,2*indices(1,2)-1,2*indices(1,3)-1);
                %
                % for k = 1:nPrim
                %     for l = 1:nPrim
                %         globalOverlap = globalOverlap + coeffMat(k,2)*coeffMat(l,2)/(coeffMat(k,1)+coeffMat(k,2))^(L+1.5);
                %
                %
                %     end
                % end
                % globalNormalization = (globalOverlapPref*globalOverlap)^0.5;

                %This part is left for compatibility for now, because I abandoned the idea of putting the specific normalization
                %coefficients for d- f- etc functions in the basis set. Now these coefficients are built in the OS rr's
                DimL = (L+1)*(L+2)/2;
                basis{nb}(6:6+DimL-1) = ones(DimL,1);

                for k = 1:nPrim
                    alpha = coeffMat(k,1);

                    normPrefactor = 2^(L)*(2/pi)^.75*(alpha^(0.75+0.5*(L)));

                    double_factorial_product = 1/double_factorial(2*indices(1,1)-1,2*indices(1,2)-1,2*indices(1,3)-1)^0.5;

                    N = normPrefactor*double_factorial_product;

                    basis{nb}(6+DimL+(k-1)*3+0) = alpha;  %alpha coefficient
                    basis{nb}(6+DimL+(k-1)*3+1) = coeffMat(k,2); %/globalNormalization;  %c coefficient
                    basis{nb}(6+DimL+(k-1)*3+2) = N; %N; for debugging
                end



            else
                disp('this case, SP, SPD, etc, is not yet dealt with for this basis set');
            end

        end
        atomBasisList(index,2) = nb; %here I tell which is the last basis for atom(index)
    end



elseif (strcmp(basisset,'6-31G') || strcmp(basisset,'6-31g'))

    for index = 1:Natoms
        x0 = AL(index,1);
        y0 = AL(index,2);
        z0 = AL(index,3);
        switch Z(index)
            case 1 %Z(index) == 1, H atom, 6-31g, 2 CFG in total
                Nel = Nel+1; %number of electrons?
                S = [18.7311370, 0.03349460
                     2.8253937, 0.23472695
                     0.6401217, 0.81375733];
                nb = nb+1;
                atomBasisList(index,1) = nb;

                %Data which is common to all primitives
                basis{nb}(1,1) = 3; %Number of primitives
                basis{nb}(2,1) = 0; %Angular momentum
                basis{nb}(3) = x0;
                basis{nb}(4) = y0;
                basis{nb}(5) = z0;

                %s-functions have a normalization factor for each primitive,
                %p-functions have the same normalization factor for all 3
                %directions.
                %But starting on d-functions, different groups of cartesian indices
                %have different normalization factors. The ratios between those
                %normalization factors is common to all primitives, so a scalar
                %normalization can be applied to the primitives, and after
                %contraction, the vector of ratios can be element-multiplied to
                %the appropriate dimensions of the contracted shell matrices.
                %So for s-functions this vector of ratios is simply a scalar 1.
                basis{nb}(6) = 1;

                basis{nb}(7) = S(1,1); %alpha value
                basis{nb}(8) = S(1,2); %contraction coefficient c
                basis{nb}(9) = (2*S(1,1)/pi)^.75; %normalization for s-function

                basis{nb}(10) = S(2,1); %alpha value
                basis{nb}(11) = S(2,2); %contraction coefficient c
                basis{nb}(12) = (2*S(2,1)/pi)^.75; %normalization for s-function

                basis{nb}(13) = S(3,1); %alpha value
                basis{nb}(14) = S(3,2); %contraction coefficient c
                basis{nb}(15) = (2*S(3,1)/pi)^.75; %normalization for s-function

                nb = nb + 1;

                basis{nb}(1,1) = 1; %Number of primitives
                basis{nb}(2,1) = 0; %Angular momentum
                basis{nb}(3) = x0;
                basis{nb}(4) = y0;
                basis{nb}(5) = z0;

                basis{nb}(6) = 1;

                basis{nb}(7) = 0.1612778; %alpha value
                basis{nb}(8) = 1; %contraction coefficient c
                basis{nb}(9) = (2*0.1612778/pi)^.75; %normalization for s-function

            case 2 %He atom, 6-31g
            Nel = Nel+2; %number of electrons?
            S = [38.4216340, 0.0237660
                 5.7780300, 0.1546790
                 1.2417740, 0.8835300];
            nb = nb+1;
            atomBasisList(index,1) = nb;

            %Data which is common to all primitives
            basis{nb}(1,1) = 3; %Number of primitives
            basis{nb}(2,1) = 0; %Angular momentum
            basis{nb}(3) = x0;
            basis{nb}(4) = y0;
            basis{nb}(5) = z0;

            basis{nb}(6) = 1;

            basis{nb}(7) = S(1,1); %alpha value
            basis{nb}(8) = S(1,2); %contraction coefficient c
            basis{nb}(9) = (2*S(1,1)/pi)^.75; %normalization for s-function

            basis{nb}(10) = S(2,1); %alpha value
            basis{nb}(11) = S(2,2); %contraction coefficient c
            basis{nb}(12) = (2*S(2,1)/pi)^.75; %normalization for s-function

            basis{nb}(13) = S(3,1); %alpha value
            basis{nb}(14) = S(3,2); %contraction coefficient c
            basis{nb}(15) = (2*S(3,1)/pi)^.75; %normalization for s-function

            nb = nb + 1;

            basis{nb}(1,1) = 1; %Number of primitives
            basis{nb}(2,1) = 0; %Angular momentum
            basis{nb}(3) = x0;
            basis{nb}(4) = y0;
            basis{nb}(5) = z0;

            basis{nb}(6) = 1;

            basis{nb}(7) = 0.2979640; %alpha value
            basis{nb}(8) = 1; %contraction coefficient c
            basis{nb}(9) = (2*0.2979640/pi)^.75; %normalization for s-function

            case 6 %C atom, 6-31g, 9 CGF in total
                Nel = Nel+6; %number of electrons
                S = [3047.5249000 0.0018347
                    457.3695100 0.0140373
                    103.9486900 0.0688426
                    29.2101550 0.2321844
                    9.2866630 0.4679413
                    3.1639270 0.3623120];
                SP1 = [7.8682724 -0.1193324 0.0689991
                       1.8812885 -0.1608542 0.3164240
                       0.5442493 1.1434564 0.7443083]; %3rd column are c coefficients for p
                SP2 = [0.1687144 1.000 1.000];

                nb = nb+1;
                atomBasisList(index,1) = nb;
                %Data which is common to all primitives
                basis{nb}(1,1) = 6; %Number of primitives
                basis{nb}(2,1) = 0; %Angular momentum
                basis{nb}(3) = x0;
                basis{nb}(4) = y0;
                basis{nb}(5) = z0;

                %s-functions have a normalization factor for each primitive,
                %p-functions have the same normalization factor for all 3
                %directions.
                %But starting on d-functions, different groups of cartesian indices
                %have different normalization factors. The ratios between those
                %normalization factors is common to all primitives, so a scalar
                %normalization can be applied to the primitives, and after
                %contraction, the vector of ratios can be element-multiplied to
                %the appropriate dimensions of the contracted shell matrices.
                %So for s-functions this vector of ratios is simply a scalar 1.
                basis{nb}(6) = 1;

                basis{nb}(7) = S(1,1); %alpha value
                basis{nb}(8) = S(1,2); %contraction coefficient c
                basis{nb}(9) = (2*S(1,1)/pi)^.75; %normalization for s-function

                basis{nb}(10) = S(2,1); %alpha value
                basis{nb}(11) = S(2,2); %contraction coefficient c
                basis{nb}(12) = (2*S(2,1)/pi)^.75; %normalization for s-function

                basis{nb}(13) = S(3,1); %alpha value
                basis{nb}(14) = S(3,2); %contraction coefficient c
                basis{nb}(15) = (2*S(3,1)/pi)^.75; %normalization for s-function

                basis{nb}(16) = S(4,1); %alpha value
                basis{nb}(17) = S(4,2); %contraction coefficient c
                basis{nb}(18) = (2*S(4,1)/pi)^.75; %normalization for s-function

                basis{nb}(19) = S(5,1); %alpha value
                basis{nb}(20) = S(5,2); %contraction coefficient c
                basis{nb}(21) = (2*S(5,1)/pi)^.75; %normalization for s-function

                basis{nb}(22) = S(6,1); %alpha value
                basis{nb}(23) = S(6,2); %contraction coefficient c
                basis{nb}(24) = (2*S(6,1)/pi)^.75; %normalization for s-function

                nb = nb + 1;

                basis{nb}(1,1) = 3; %Number of primitives
                basis{nb}(2,1) = 0; %Angular momentum
                basis{nb}(3) = x0;
                basis{nb}(4) = y0;
                basis{nb}(5) = z0;

                basis{nb}(6) = 1;

                basis{nb}(7) = SP1(1,1); %alpha value
                basis{nb}(8) = SP1(1,2); %contraction coefficient c
                basis{nb}(9) = (2*SP1(1,1)/pi)^.75; %normalization for s-function

                basis{nb}(10) = SP1(2,1); %alpha value
                basis{nb}(11) = SP1(2,2); %contraction coefficient c
                basis{nb}(12) = (2*SP1(2,1)/pi)^.75; %normalization for s-function

                basis{nb}(13) = SP1(3,1); %alpha value
                basis{nb}(14) = SP1(3,2); %contraction coefficient c
                basis{nb}(15) = (2*SP1(3,1)/pi)^.75; %normalization for s-function

                nb = nb + 1;

                basis{nb}(1,1) = 3; %Number of primitives
                basis{nb}(2,1) = 1; %Angular momentum, p-function
                basis{nb}(3) = x0;
                basis{nb}(4) = y0;
                basis{nb}(5) = z0;

                basis{nb}(6) = 1;
                basis{nb}(7) = 1;
                basis{nb}(8) = 1;

                basis{nb}(9) = SP1(1,1); %alpha value
                basis{nb}(10) = SP1(1,3); %contraction coefficient c
                alpha = SP1(1,1);
                basis{nb}(11) = 2^(1)*(2/pi)^.75*alpha^(0.75+0.5*(1))/(double_factorial(1,0,0))^.5; %normalization for p-function

                basis{nb}(12) = SP1(2,1); %alpha value
                basis{nb}(13) = SP1(2,3); %contraction coefficient c
                alpha = SP1(2,1);
                basis{nb}(14) = 2^(1)*(2/pi)^.75*alpha^(0.75+0.5*(1))/(double_factorial(1,0,0))^.5; %normalization for p-function

                basis{nb}(15) = SP1(3,1); %alpha value
                basis{nb}(16) = SP1(3,3); %contraction coefficient c
                alpha = SP1(3,1);
                basis{nb}(17) = 2^(1)*(2/pi)^.75*alpha^(0.75+0.5*(1))/(double_factorial(1,0,0))^.5; %normalization for p-function

                nb = nb + 1;

                basis{nb}(1,1) = 1; %Number of primitives
                basis{nb}(2,1) = 0; %Angular momentum
                basis{nb}(3) = x0;
                basis{nb}(4) = y0;
                basis{nb}(5) = z0;

                basis{nb}(6) = 1;

                basis{nb}(7) = SP2(1,1); %alpha value
                basis{nb}(8) = SP2(1,2); %contraction coefficient c
                basis{nb}(9) = (2*SP2(1,1)/pi)^.75; %normalization for s-function

                nb = nb + 1;

                basis{nb}(1,1) = 1; %Number of primitives
                basis{nb}(2,1) = 1; %Angular momentum
                basis{nb}(3) = x0;
                basis{nb}(4) = y0;
                basis{nb}(5) = z0;

                basis{nb}(6) = 1;
                basis{nb}(7) = 1;
                basis{nb}(8) = 1;

                basis{nb}(9) = SP2(1,1); %alpha value
                basis{nb}(10) = SP2(1,3); %contraction coefficient c
                alpha = SP2(1,1);
                basis{nb}(11) = 2^(1)*(2/pi)^.75*alpha^(0.75+0.5*(1))/(double_factorial(1,0,0))^.5; %normalization for p-function
        end
        atomBasisList(index,2) = nb;
    end

elseif (strcmp(basisset,'STO-3G') || strcmp(basisset,'sto-3g'))

    for index = 1:Natoms
                x0 = AL(index,1);
                y0 = AL(index,2);
                z0 = AL(index,3);
                switch Z(index)
                    case 1 %Z(index) == 1, H atom, 6-31g, 2 CFG in total
                        Nel = Nel+1; %number of electrons?
                        S = [3.42525091,0.15432897
                             0.62391373,0.53532814
                             0.16885540,0.44463454];

                        nb = nb+1;
                        atomBasisList(index,1) = nb;

                        basis{nb}(1,1) = 3; %Number of primitives
                        basis{nb}(2,1) = 0; %Angular momentum
                        basis{nb}(3) = x0;
                        basis{nb}(4) = y0;
                        basis{nb}(5) = z0;

                        basis{nb}(6) = 1;

                        basis{nb}(7) = S(1,1); %alpha value
                        basis{nb}(8) = S(1,2); %contraction coefficient c
                        basis{nb}(9) = (2*S(1,1)/pi)^.75; %normalization for s-function

                        basis{nb}(10) = S(2,1); %alpha value
                        basis{nb}(11) = S(2,2); %contraction coefficient c
                        basis{nb}(12) = (2*S(2,1)/pi)^.75; %normalization for s-function

                        basis{nb}(13) = S(3,1); %alpha value
                        basis{nb}(14) = S(3,2); %contraction coefficient c
                        basis{nb}(15) = (2*S(3,1)/pi)^.75; %normalization for s-function

                    case 2 %He atom, 6-31g
                        Nel = Nel+2; %number of electrons?
                        S = [6.36242139,0.15432897
                             1.15892300,0.53532814
                             0.31364979,0.44463454];

                        nb = nb+1;
                        atomBasisList(index,1) = nb;
                        basis{nb}(1,1) = 3; %Number of primitives
                        basis{nb}(2,1) = 0; %Angular momentum
                        basis{nb}(3) = x0;
                        basis{nb}(4) = y0;
                        basis{nb}(5) = z0;

                        basis{nb}(6) = 1;

                        basis{nb}(7) = S(1,1); %alpha value
                        basis{nb}(8) = S(1,2); %contraction coefficient c
                        basis{nb}(9) = (2*S(1,1)/pi)^.75; %normalization for s-function

                        basis{nb}(10) = S(2,1); %alpha value
                        basis{nb}(11) = S(2,2); %contraction coefficient c
                        basis{nb}(12) = (2*S(2,1)/pi)^.75; %normalization for s-function

                        basis{nb}(13) = S(3,1); %alpha value
                        basis{nb}(14) = S(3,2); %contraction coefficient c
                        basis{nb}(15) = (2*S(3,1)/pi)^.75; %normalization for s-function

                    case 6 %C atom, STO-3G
                        Nel = Nel+6; %number of electrons?
                        S = [71.6168370, 0.15432897
                             13.0450960,0.53532814
                             3.5305122,0.44463454];
                        SP = [2.9412494,-0.09996723,0.15591627
                              0.6834831,0.39951283,0.60768372
                              0.2222899,0.70011547,0.39195739]; %3rd column are c coefficients for p

                        nb = nb+1;
                        atomBasisList(index,1) = nb;
                        basis{nb}(1,1) = 3; %Number of primitives
                        basis{nb}(2,1) = 0; %Angular momentum
                        basis{nb}(3) = x0;
                        basis{nb}(4) = y0;
                        basis{nb}(5) = z0;

                        basis{nb}(6) = 1;

                        basis{nb}(7) = S(1,1); %alpha value
                        basis{nb}(8) = S(1,2); %contraction coefficient c
                        basis{nb}(9) = (2*S(1,1)/pi)^.75; %normalization for s-function

                        basis{nb}(10) = S(2,1); %alpha value
                        basis{nb}(11) = S(2,2); %contraction coefficient c
                        basis{nb}(12) = (2*S(2,1)/pi)^.75; %normalization for s-function

                        basis{nb}(13) = S(3,1); %alpha value
                        basis{nb}(14) = S(3,2); %contraction coefficient c
                        basis{nb}(15) = (2*S(3,1)/pi)^.75; %normalization for s-function

                        nb = nb+1;

                        basis{nb}(1,1) = 3; %Number of primitives
                        basis{nb}(2,1) = 0; %Angular momentum
                        basis{nb}(3) = x0;
                        basis{nb}(4) = y0;
                        basis{nb}(5) = z0;

                        basis{nb}(6) = 1;

                        basis{nb}(7) = SP(1,1); %alpha value
                        basis{nb}(8) = SP(1,2); %contraction coefficient c
                        basis{nb}(9) = (2*SP(1,1)/pi)^.75; %normalization for s-function

                        basis{nb}(10) = SP(2,1); %alpha value
                        basis{nb}(11) = SP(2,2); %contraction coefficient c
                        basis{nb}(12) = (2*SP(2,1)/pi)^.75; %normalization for s-function

                        basis{nb}(13) = SP(3,1); %alpha value
                        basis{nb}(14) = SP(3,2); %contraction coefficient c
                        basis{nb}(15) = (2*SP(3,1)/pi)^.75; %normalization for s-function

                        nb = nb+1;

                        basis{nb}(1,1) = 3; %Number of primitives
                        basis{nb}(2,1) = 1; %Angular momentum
                        basis{nb}(3) = x0;
                        basis{nb}(4) = y0;
                        basis{nb}(5) = z0;

                        basis{nb}(6) = 1;
                        basis{nb}(7) = 1;
                        basis{nb}(8) = 1;

                        basis{nb}(9) = SP(1,1); %alpha value
                        basis{nb}(10) = SP(1,3); %contraction coefficient c
                        alpha = SP(1,1);
                        basis{nb}(11) = 2^(1)*(2/pi)^.75*alpha^(0.75+0.5*(1))/(double_factorial(1,0,0))^.5; %normalization for p-function

                        basis{nb}(12) = SP(2,1); %alpha value
                        basis{nb}(13) = SP(2,3); %contraction coefficient c
                        alpha = SP(2,1);
                        basis{nb}(14) = 2^(1)*(2/pi)^.75*alpha^(0.75+0.5*(1))/(double_factorial(1,0,0))^.5; %normalization for p-function

                        basis{nb}(15) = SP(3,1); %alpha value
                        basis{nb}(16) = SP(3,3); %contraction coefficient c
                        alpha = SP(3,1);
                        basis{nb}(17) = 2^(1)*(2/pi)^.75*alpha^(0.75+0.5*(1))/(double_factorial(1,0,0))^.5; %normalization for p-function

                    case 7 %N atom, STO-3G
                        Nel = Nel+7; %number of electrons?
                        S = [99.1061690,0.15432897
                             18.0523120,0.53532814
                             4.8856602,0.44463454];
                        SP = [3.7804559,-0.09996723,0.15591627
                              0.8784966,0.39951283,0.60768372
                              0.2857144,0.70011547,0.39195739]; %3rd column are c coefficients for p

                        nb = nb+1;
                        atomBasisList(index,1) = nb;
                        basis{nb}(1,1) = 3; %Number of primitives
                        basis{nb}(2,1) = 0; %Angular momentum
                        basis{nb}(3) = x0;
                        basis{nb}(4) = y0;
                        basis{nb}(5) = z0;

                        basis{nb}(6) = 1;

                        basis{nb}(7) = S(1,1); %alpha value
                        basis{nb}(8) = S(1,2); %contraction coefficient c
                        basis{nb}(9) = (2*S(1,1)/pi)^.75; %normalization for s-function

                        basis{nb}(10) = S(2,1); %alpha value
                        basis{nb}(11) = S(2,2); %contraction coefficient c
                        basis{nb}(12) = (2*S(2,1)/pi)^.75; %normalization for s-function

                        basis{nb}(13) = S(3,1); %alpha value
                        basis{nb}(14) = S(3,2); %contraction coefficient c
                        basis{nb}(15) = (2*S(3,1)/pi)^.75; %normalization for s-function

                        nb = nb+1;

                        basis{nb}(1,1) = 3; %Number of primitives
                        basis{nb}(2,1) = 0; %Angular momentum
                        basis{nb}(3) = x0;
                        basis{nb}(4) = y0;
                        basis{nb}(5) = z0;

                        basis{nb}(6) = 1;

                        basis{nb}(7) = SP(1,1); %alpha value
                        basis{nb}(8) = SP(1,2); %contraction coefficient c
                        basis{nb}(9) = (2*SP(1,1)/pi)^.75; %normalization for s-function

                        basis{nb}(10) = SP(2,1); %alpha value
                        basis{nb}(11) = SP(2,2); %contraction coefficient c
                        basis{nb}(12) = (2*SP(2,1)/pi)^.75; %normalization for s-function

                        basis{nb}(13) = SP(3,1); %alpha value
                        basis{nb}(14) = SP(3,2); %contraction coefficient c
                        basis{nb}(15) = (2*SP(3,1)/pi)^.75; %normalization for s-function

                        nb = nb+1;

                        basis{nb}(1,1) = 3; %Number of primitives
                        basis{nb}(2,1) = 1; %Angular momentum
                        basis{nb}(3) = x0;
                        basis{nb}(4) = y0;
                        basis{nb}(5) = z0;

                        basis{nb}(6) = 1;
                        basis{nb}(7) = 1;
                        basis{nb}(8) = 1;

                        basis{nb}(9) = SP(1,1); %alpha value
                        basis{nb}(10) = SP(1,3); %contraction coefficient c
                        alpha = SP(1,1);
                        basis{nb}(11) = 2^(1)*(2/pi)^.75*alpha^(0.75+0.5*(1))/(double_factorial(1,0,0))^.5; %normalization for p-function

                        basis{nb}(12) = SP(2,1); %alpha value
                        basis{nb}(13) = SP(2,3); %contraction coefficient c
                        alpha = SP(2,1);
                        basis{nb}(14) = 2^(1)*(2/pi)^.75*alpha^(0.75+0.5*(1))/(double_factorial(1,0,0))^.5; %normalization for p-function

                        basis{nb}(15) = SP(3,1); %alpha value
                        basis{nb}(16) = SP(3,3); %contraction coefficient c
                        alpha = SP(3,1);
                        basis{nb}(17) = 2^(1)*(2/pi)^.75*alpha^(0.75+0.5*(1))/(double_factorial(1,0,0))^.5; %normalization for p-function

                    case 8 %O atom, STO-3G
                        Nel = Nel+8; %number of electrons?
                        S = [130.7093200,0.15432897
                             23.8088610,0.53532814
                             6.4436083,0.44463454];
                        SP = [5.0331513,-0.09996723,0.15591627
                              1.1695961,0.39951283,0.60768372
                              0.3803890,0.70011547,0.39195739]; %3rd column are c coefficients for p

                        nb = nb+1;
                        atomBasisList(index,1) = nb;
                        basis{nb}(1,1) = 3; %Number of primitives
                        basis{nb}(2,1) = 0; %Angular momentum
                        basis{nb}(3) = x0;
                        basis{nb}(4) = y0;
                        basis{nb}(5) = z0;

                        basis{nb}(6) = 1;

                        basis{nb}(7) = S(1,1); %alpha value
                        basis{nb}(8) = S(1,2); %contraction coefficient c
                        basis{nb}(9) = (2*S(1,1)/pi)^.75; %normalization for s-function

                        basis{nb}(10) = S(2,1); %alpha value
                        basis{nb}(11) = S(2,2); %contraction coefficient c
                        basis{nb}(12) = (2*S(2,1)/pi)^.75; %normalization for s-function

                        basis{nb}(13) = S(3,1); %alpha value
                        basis{nb}(14) = S(3,2); %contraction coefficient c
                        basis{nb}(15) = (2*S(3,1)/pi)^.75; %normalization for s-function

                        nb = nb+1;

                        basis{nb}(1,1) = 3; %Number of primitives
                        basis{nb}(2,1) = 0; %Angular momentum
                        basis{nb}(3) = x0;
                        basis{nb}(4) = y0;
                        basis{nb}(5) = z0;

                        basis{nb}(6) = 1;

                        basis{nb}(7) = SP(1,1); %alpha value
                        basis{nb}(8) = SP(1,2); %contraction coefficient c
                        basis{nb}(9) = (2*SP(1,1)/pi)^.75; %normalization for s-function

                        basis{nb}(10) = SP(2,1); %alpha value
                        basis{nb}(11) = SP(2,2); %contraction coefficient c
                        basis{nb}(12) = (2*SP(2,1)/pi)^.75; %normalization for s-function

                        basis{nb}(13) = SP(3,1); %alpha value
                        basis{nb}(14) = SP(3,2); %contraction coefficient c
                        basis{nb}(15) = (2*SP(3,1)/pi)^.75; %normalization for s-function

                        nb = nb+1;

                        basis{nb}(1,1) = 3; %Number of primitives
                        basis{nb}(2,1) = 1; %Angular momentum
                        basis{nb}(3) = x0;
                        basis{nb}(4) = y0;
                        basis{nb}(5) = z0;

                        basis{nb}(6) = 1;
                        basis{nb}(7) = 1;
                        basis{nb}(8) = 1;

                        basis{nb}(9) = SP(1,1); %alpha value
                        basis{nb}(10) = SP(1,3); %contraction coefficient c
                        alpha = SP(1,1);
                        basis{nb}(11) = 2^(1)*(2/pi)^.75*alpha^(0.75+0.5*(1))/(double_factorial(1,0,0))^.5; %normalization for p-function

                        basis{nb}(12) = SP(2,1); %alpha value
                        basis{nb}(13) = SP(2,3); %contraction coefficient c
                        alpha = SP(2,1);
                        basis{nb}(14) = 2^(1)*(2/pi)^.75*alpha^(0.75+0.5*(1))/(double_factorial(1,0,0))^.5; %normalization for p-function

                        basis{nb}(15) = SP(3,1); %alpha value
                        basis{nb}(16) = SP(3,3); %contraction coefficient c
                        alpha = SP(3,1);

                end
                atomBasisList(index,2) = nb;
    end
else
    error('Basis set not recognized');
end

basis = basis(1:nb,1);
end
