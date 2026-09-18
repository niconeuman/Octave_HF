function [Coulomb,Exchangea,Exchangeb] = uBuildJK_nosym(Gprev,Da,Db,Dm1a,Dm1b,basis,pair_data2,Shell_List,Boys_Table,S)

%This function will do what Build_ERI_OS did, but then it will contract the shell quartets using the density matrix
%to generate the

thrOverlap = 1e-8;
thrD = 1e-10;
nshells = Shell_List(end,3);
Ncont = Shell_List(end,2);
J = zeros(Ncont,Ncont);
Ka = zeros(Ncont,Ncont);
Kb = zeros(Ncont,Ncont);

Jfactor = 1;
Kfactor = 1; %for unrestricted

D = Da + Db;

%Debug function to see if all the matrix is correctly covered
%D = ones(size(D));

for shella = 1:nshells
    mu_begin = Shell_List(shella,1);
    mu_end = Shell_List(shella,2);
    lengthmu = mu_end-mu_begin+1;
    %aa = Shells_List(a,2);
    basis_a = basis{shella};
    La = basis_a(2);
    for shellb = 1:nshells
        nu_begin = Shell_List(shellb,1);
        nu_end = Shell_List(shellb,2);
        lengthnu = nu_end-nu_begin+1;
        %ab = Shells_List(b,2);
        basis_b = basis{shellb};
        Lb = basis_b(2);
        %ab_data = pair_data{a,b};
        if (any(any(abs(S(mu_begin:mu_end,nu_begin:nu_end))>thrOverlap)) && any(any(abs(D(mu_begin:mu_end,nu_begin:nu_end))>thrD)))
        for shellc = 1:nshells
            ka_begin = Shell_List(shellc,1);
            ka_end = Shell_List(shellc,2);
            lengthka = ka_end-ka_begin+1;
            %ac = Shells_List(c,2);
            basis_c = basis{shellc};
            Lc = basis_c(2);

                for shelld = 1:nshells
                    la_begin = Shell_List(shelld,1);
                    la_end = Shell_List(shelld,2);
                    lengthla = la_end-la_begin+1;
                    %ad = Shells_List(d,2);
                    basis_d = basis{shelld};
                    Ld = basis_d(2);
                    %cd_data = pair_data{c,d};

                    if (any(any(abs(S(ka_begin:ka_end,la_begin:la_end))>thrOverlap)) && any(any(abs(D(ka_begin:ka_end,la_begin:la_end))>thrD))) %this later can be changed by a screening option
                        if (La == 0 && Lb == 0 && Lc == 0 && Ld == 0) %SS|SS integrals

                            gSSSSNValues = primitiveFactorsSSSS_3(basis_a,basis_b,basis_c,basis_d,Boys_Table,pair_data2,shella,shellb,shellc,shelld);
                            gmnkl = sum(gSSSSNValues(:,1));
                        else
                            gmnkl = shellOS(basis_a,basis_b,basis_c,basis_d,La,Lb,Lc,Ld,Boys_Table,pair_data2,shella,shellb,shellc,shelld);
                        endif

                        intmu = mu_begin:mu_end;
                        intnu = nu_begin:nu_end;
                        intka = ka_begin:ka_end;
                        intla = la_begin:la_end;

##                        %Coulomb contribution, not making use of symmetry. Some of them will be counted several times
##                        %mn kl
##                        Dtemp = D(intka,intla);
##                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
##                        J(intmu,intnu) = J(intmu,intnu) +  squeeze(sum(sum(Jfactor*(Dres).*gmnkl,4),3));
##
##                        %Exchange contribution, all 2nd and 3rd indices have to be exchanged
##                        #gmknl = permute(gmnkl,[1 3 2 4]);
##                        #According to ChatGPT, only mnkl should be used, but my permute is due to the shape of the matrices.
##
##                        %symFactor = 1;
##
##                        %Construction of the alpha orbital matrix
##                        %mk nl
##                        Dtemp = Da(intnu,intla);
##                        Dres = reshape(Dtemp,[1,size(Dtemp,1),1,size(Dtemp,2)]);
##                        Ka(intmu,intka) = Ka(intmu,intka) -  squeeze(sum(sum(Kfactor*(Dres).*gmnkl,4),2));
##
##                        %Construction of the beta orbital matrix
##                        %mk nl
##                        Dtemp = Db(intnu,intla);
##                        Dres = reshape(Dtemp,[1,size(Dtemp,1),1,size(Dtemp,2)]);
##                        Kb(intmu,intka) = Kb(intmu,intka) -  squeeze(sum(sum(Kfactor*(Dres).*gmnkl,4),2));
                        for imu = 1:length(intmu)
                            mu = intmu(imu);

                            for inu = 1:length(intnu)
                                nu = intnu(inu);

                                for ika = 1:length(intka)
                                    ka = intka(ika);

                                    for ila = 1:length(intla)
                                        la = intla(ila);

                                        eri = gmnkl(imu,inu,ika,ila);

##                                        if rand(1) < 0.0001
##                                          disp("The eri is: ");
##                                          disp(eri);
##                                          disp("the shell indices are: ");
##                                          disp([shella, shellb, shellc, shelld]);
##                                          disp("The D matrices are ");
##                                          disp([Da, Db]);
##                                        end
                                        % Coulomb
                                        J(mu,nu) = J(mu,nu) + ...
                                            (Da(ka,la) + Db(ka,la))*eri;

                                        % Exchange
                                        Ka(mu,ka) = Ka(mu,ka) - Da(nu,la)*eri;
                                        Kb(mu,ka) = Kb(mu,ka) - Db(nu,la)*eri;

                                    end
                                end
                            end
                        end

                    else
                        continue;
                    endif
                endfor

        endfor
        else %line 35
            continue
        endif %line 35
    endfor
endfor

%Because I only have the lowerdiagonal part
%Jtril = tril(J); %to remove upper diagonal elements which occurred due to the use of shells.
%J = tril(Jtril,-1)'+Jtril;

%Ktril = tril(K); %to remove upper diagonal elements which occurred due to the use of shells.
%K = tril(Ktril,-1)'+Ktril;

Coulomb = J;
Exchangea = Ka;
Exchangeb = Kb;


end
