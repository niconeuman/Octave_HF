function [S,T,Ssph,Tsph] = Build_OS1e(basis,pair_data2,Shell_List,Sph_Shell_List)

Ncont = Shell_List(end,2);
S = zeros(Ncont,Ncont);
T = zeros(Ncont,Ncont);

NcontSph = Sph_Shell_List(end,2);
Ssph = zeros(NcontSph,NcontSph);
Tsph = zeros(NcontSph,NcontSph);

nshells = size(basis,1);

for shella = 1:nshells
    mu_begin = Shell_List(shella,1);
    mu_end   = Shell_List(shella,2);

    mu_sph_begin = Sph_Shell_List(shella,1);
    mu_sph_end   = Sph_Shell_List(shella,2);

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


    for shellb = 1:nshells
        nu_begin = Shell_List(shellb,1);
        nu_end   = Shell_List(shellb,2);

        nu_sph_begin = Sph_Shell_List(shellb,1);
        nu_sph_end   = Sph_Shell_List(shellb,2);

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
        try
                if (size(aValues) ~= size(bValues)) || (size(aValues) ~= size(pValues))
                        disp('aValues,bValues and/or pValues have different sizes');
                end
                if (Lb > La)
                        [SabT,TabT] = OS1e(KabValues,RPBValues,RPAValues,bValues,aValues,pValues,WeightValues,R12sq,Lb,La);
                        Sab = SabT';
                        Tab = TabT';
                else
                        [Sab,Tab] = OS1e(KabValues,RPAValues,RPBValues,aValues,bValues,pValues,WeightValues,R12sq,La,Lb);
                end

                S(mu_begin:mu_end,nu_begin:nu_end) = Sab;
                T(mu_begin:mu_end,nu_begin:nu_end) = Tab;

                Ssph(mu_sph_begin:mu_sph_end,nu_sph_begin:nu_sph_end) = cart2sph(Sab);
                Tsph(mu_sph_begin:mu_sph_end,nu_sph_begin:nu_sph_end) = cart2sph(Tab);
        catch
                disp('Error in shells')
                disp([num2str(shella), ' and ', num2str(shellb)])
                disp('Sab and Tab are:');
                disp([Sab,Tab]);
                % disp('aValues');
                % disp(size(aValues));
                % disp(aValues);
                % disp('bValues');
                % disp(size(bValues));
                % disp(bValues);
                % disp('pValues');
                % disp(size(pValues));
                % disp(pValues);

                shella = nshells+1;
                break;

        end



    end
end


end
