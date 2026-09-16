function [MinE,E,ncycle,D,Dinit,G,epsilon,F,Fprime] = SCF(H0,T,Nuclear_Attraction,gabcd,S,Nel,Shell_List)
maxcycles = 100;
converged = 1;
ncycle = 0;

[V,DV] = eig(S);
X = V*DV^(-0.5)*V'; %transforms to an orthonormal basis

%This is the core guess, which is not good
F = H0;
%This is an atomic guess.

% Fdiag = [2 2 2/3 2/3 2/3 2 2 2/3 2/3 2/3 1 1];
% F = diag(Fdiag);

%  Fdiag = [2 1 1/3 1/3 1/3 1 1/3 1/3 1/3 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
%  F = diag(Fdiag);

%This is a polarized guess, where nuclear attraction is screened
F = T + 0.3*Nuclear_Attraction;

Fprime = X'*F*X;

[Cprime,epsilon] = eig(Fprime);

C = X*Cprime;

%[C,epsilon] = eig(F,S);
[C,epsilon] = Sort_Eigs(C,epsilon);

Dinit = Build_Density(C,Nel);

E = zeros(maxcycles,1);

%D = 0.2*eye(size(H0))+0.8*Dinit;
Dprev = Dinit;
Dm2 = zeros(size(H0));
Dprev = 0.9*eye(size(H0))+0.1*Dinit;
Gprev = zeros(size(H0));
%Damping
Damp = 0.8;

while ncycle < maxcycles && converged ~= 0
    ncycle = ncycle + 1;
    G = Build_Coulomb_Exchange(Gprev,Dprev,Dm2,gabcd,Shell_List);
    F = H0 +  G;
    Fprime = X'*F*X; %'
    [Cprime,epsilon] = eig(Fprime);
    C = X*Cprime;
    %[C,epsilon] = eig(F,S);
    [C,epsilon] = Sort_Eigs(C,epsilon);


    Dnew = Build_Density(C,Nel);
    if (ncycle > 1) && (abs(E(ncycle)-E(ncycle-1))<1e-2)
        Damp = 1;
    end

    D = Damp*Dnew+(1-Damp)*Dprev;
%     D(3:5,3:5) = D(3:5,3:5).*eye(3);
%     D(7:9,7:9) = D(7:9,7:9).*eye(3);
    Dm2 = Dprev;
    Dprev = D;
    Gprev = G;
    E(ncycle) = Fock_Energy(D,H0,F);
        if(ncycle > 1) && abs(E(ncycle)-E(ncycle-1)) < 1e-6
            converged = 0;
        end
end
E = E(1:ncycle);
MinE = E(ncycle);
end
