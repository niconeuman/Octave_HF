function Ven = Build_OS_EN(basis,pair_data2,Shell_List,AL,Z,Boys_Table)

Ncont = Shell_List(end,2);
Ven = zeros(Ncont,Ncont);

nshells = size(basis,1);

for shella = 1:nshells
    mu_begin = Shell_List(shella,1);
    mu_end   = Shell_List(shella,2);
    basis_a = basis{shella};

    La = basis{shella}(2);
    Dima = (La+1)*(La+2)/2;

    Nprima = basis{shella}(1);
    alphaa = zeros(Nprima,1);
    for na = 1:Nprima
        alphaa(na) = basis{shella}(6+Dima+(na-1)*3+0);
    end
    %disp(alphaa);
    xa = basis{shella}(3);
    ya = basis{shella}(4);
    za = basis{shella}(5);


    for shellb = 1:nshells %This is now calculating the whole matrix, I need to use symmetry
        nu_begin = Shell_List(shellb,1);
        nu_end   = Shell_List(shellb,2);
        basis_b = basis{shellb};

        Lb = basis{shellb}(2);
        Dimb = (Lb+1)*(Lb+2)/2;

        Nprimb = basis{shellb}(1);
        alphab = zeros(Nprimb,1);
        for nb = 1:Nprimb
            alphab(nb) = basis{shellb}(6+Dimb+(nb-1)*3+0);
        end

        %pair_data{a,b} = [pValues PxValues PyValues PzValues KabValues WeightValuesab RPAValues RPBValues pPxValues pPyValues pPzValues];
        KabValues = pair_data2{shella,shellb}(:,5);
        RPAValues = pair_data2{shella,shellb}(:,7:9);
        RPBValues = pair_data2{shella,shellb}(:,10:12);
        pValues = pair_data2{shella,shellb}(:,1);
        WeightValues = pair_data2{shella,shellb}(:,6);

        onesa = ones(Nprima,1);
        onesb = ones(Nprimb,1);

        aValues = kron(alphaa,onesb); %Have to check that these Kronecker products give me the correct ordering
        bValues = kron(onesa,alphab);

        xb = basis{shellb}(3);
        yb = basis{shellb}(4);
        zb = basis{shellb}(5);

        R12sq = (xa-xb)^2+(ya-yb)^2+(za-zb)^2;

        %The function OS1e already produces the contracted shell doublet matrices
        [ENab] = shellEN(basis_a,basis_b,La,Lb,pair_data2,AL,Z,Boys_Table);

        Ven(mu_begin:mu_end,nu_begin:nu_end) = ENab;


    end
end


end
