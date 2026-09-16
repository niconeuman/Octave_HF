function [Min_Energy,E,ncycle,Da,Db,Dinita,Dinitb,Coulomb,Exchangea,Exchangeb,epsilona,epsilonb,Fa,Fprimea,Fb,Fprimeb] = uSCF(Hcore0,T,Ven,S,Nela,Nelb,basis,pair_data2,Shell_List,Boys_Table,gabcd,Z,atomBasisList)
%This function constructs the unrestricted Fock matrices

maxcycles = 50;
converged = 1;
ncycle = 1;

[V,DV] = eig(S);
X = V*real(DV^(-0.5))*V'; %transforms to an orthonormal basis

%This is the core guess, which is not good
F = T + 0.5*Ven;
Fa = F;
Fb = F;

Fprimea = X'*Fa*X;
Fprimeb = X'*Fb*X;

[Cprimea,epsilona] = eig(Fprimea);
[Cprimeb,epsilonb] = eig(Fprimeb);

Ca = X*Cprimea;
Cb = X*Cprimeb;

%[C,epsilon] = eig(F,S);
[Ca,epsilona] = Sort_Eigs(Ca,epsilona);
[Cb,epsilonb] = Sort_Eigs(Cb,epsilonb);

[Dinita,Dinitb] = uBuildDensity(Ca,Cb,Nela,Nelb);

E = zeros(maxcycles,1);

%D = 0.2*eye(size(H0))+0.8*Dinit;
%Dpreva = Dinita;
%Dprevb = Dinitb;
Dm3a = zeros(size(Hcore0));
Dm3b = zeros(size(Hcore0));
Dm2a = zeros(size(Hcore0));
Dm2b = zeros(size(Hcore0));
Dm1a = 0.0*eye(size(Hcore0))+1.0*Dinita;
Dm1b = 0.0*eye(size(Hcore0))+1.0*Dinitb;
Gprev = zeros(size(Hcore0));
%Damping
Damp = 0.8;

%I consider the last 2 iterations for DIIS
residual = zeros(length(Dinita)*length(Dinita),3);

%[Da,Db] = buildInitDens(basis,Shell_List,Z,atomBasisList,1);

disp('   ncycle   E(current)   E(previous)  Difference DampFactor');

while ((ncycle < maxcycles) && (converged ~= 0))
    ncycle = ncycle + 1;

    %This generates the Coulomb and Exchange matrices, which is the most time consuming step
    [Coulomb,Exchangea,Exchangeb] = uBuildJK(Gprev,Dm1a,Dm1b,Dm2a,Dm2b,basis,pair_data2,Shell_List,Boys_Table,S);

    Fa = Hcore0 + Coulomb + Exchangea;
    Fb = Hcore0 + Coulomb + Exchangeb;

    Fprimea = X'*Fa*X;
    Fprimeb = X'*Fb*X;

    [Cprimea,epsilona] = eig(Fprimea);
    [Cprimeb,epsilonb] = eig(Fprimeb);

    Ca = X*Cprimea;
    Cb = X*Cprimeb;

    [Ca,epsilona] = Sort_Eigs(Ca,epsilona);
    [Cb,epsilonb] = Sort_Eigs(Cb,epsilonb);

    %D0 is the current density, which comes from the last Fock matrix diagonalization.
    [D0a,D0b] = uBuildDensity(Ca,Cb,Nela,Nelb);

    if (ncycle < 3)
        Da = Damp*D0a+(1-Damp)*Dm1a;
        Db = Damp*D0b+(1-Damp)*Dm1b;

    else

        %Calculate the residual
        residual(:,3) = D0a(:)-Dm1a(:);
        residual(:,2) = Dm1a(:)-Dm2a(:);
        residual(:,1) = Dm2a(:)-Dm3a(:);
        %residual(:,1) = Dnewa(:)-Dm2a(:);

        %Equation to solve
        %AugMatrix*cVector = aVector;
        Bmatrix(1,1) = residual(:,1)'*residual(:,1);
        Bmatrix(1,2) = residual(:,1)'*residual(:,2);
        Bmatrix(1,3) = residual(:,1)'*residual(:,3);
        Bmatrix(2,1) = Bmatrix(1,2);
        Bmatrix(2,2) = residual(:,2)'*residual(:,2);
        Bmatrix(2,3) = residual(:,2)'*residual(:,3);
        Bmatrix(3,1) = Bmatrix(1,3);
        Bmatrix(3,2) = Bmatrix(2,3);
        Bmatrix(3,3) = residual(:,3)'*residual(:,3);

        OnesVector = ones(length(Bmatrix),1);
        AugMatrix = [Bmatrix -OnesVector;
                    -OnesVector' 0];

        aVector = zeros(length(OnesVector)+1,1);
        aVector(end) = -1;

        %disp('aVector');
        %disp(aVector);
        cVector = AugMatrix\aVector;

        %disp('cVector');
        %disp(cVector);
        cVector = cVector(1:end-1);

        %Da and Db are the new densities, which will be used for constructing the Fock matrix in the next cycle
        Da = (cVector(3)*D0a+cVector(2)*Dm1a+cVector(1)*Dm2a);
        Db = (cVector(3)*D0b+cVector(2)*Dm1b+cVector(1)*Dm2b);


    end

    Dm3a = Dm2a;
    Dm3b = Dm2b;

    Dm2a = Dm1a;
    Dm2b = Dm1b;

    Dm1a = Da;
    Dm1b = Db;

    E(ncycle) = uFockEnergy(Da,Db,Fa,Fb,Hcore0);

    if ((ncycle > 1) && (abs(E(ncycle)-E(ncycle-1))<1e-2))
        Damp = 1;
    end

    if((ncycle > 1) && (abs(E(ncycle)-E(ncycle-1)) < 1e-5))
        converged = 0;
        disp('cVector');
        disp(cVector);
    end

    format short g
    disp([round(ncycle), E(ncycle),E(ncycle-1),E(ncycle)-E(ncycle-1),Damp]);
end
E = E(1:ncycle);
Min_Energy = E(ncycle);
end
