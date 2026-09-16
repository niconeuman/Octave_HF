function [Coulomb,Exchange] = Build_Coulomb_Exchange(Gprev,D,Dm1,basis,pair_data2,Shell_List,Boys_Table)

%This function will do what Build_ERI_OS did, but then it will contract the shell quartets using the density matrix
%to generate the

nshells = Shell_List(end,3);
Ncont = Shell_List(end,2);
J = zeros(Ncont,Ncont);
K = zeros(Ncont,Ncont);

Jfactor = 2;
Kfactor = 0.5;
%Debug function to see if all the matrix is correctly covered
%D = ones(size(D));

for shella = 1:nshells
    mu_begin = Shell_List(shella,1);
    mu_end = Shell_List(shella,2);
    lengthmu = mu_end-mu_begin+1;
    %aa = Shells_List(a,2);
    basis_a = basis{shella};
    La = basis_a(2);
    for shellb = 1:shella
        nu_begin = Shell_List(shellb,1);
        nu_end = Shell_List(shellb,2);
        lengthnu = nu_end-nu_begin+1;
        %ab = Shells_List(b,2);
        basis_b = basis{shellb};
        Lb = basis_b(2);
        %ab_data = pair_data{a,b};

        if shellb == shella
            factormn = 1 - 0.5*ones(lengthmu,lengthnu);
            %factormn = 0.5*ones(lengthmu,lengthnu);
        else
            factormn = 1*ones(lengthmu,lengthnu);
        end


        for shellc = 1:shella
            ka_begin = Shell_List(shellc,1);
            ka_end = Shell_List(shellc,2);
            lengthka = ka_end-ka_begin+1;
            %ac = Shells_List(c,2);
            basis_c = basis{shellc};
            Lc = basis_c(2);
            if (shellc == shella)

            factormk = 1 - 0.5*ones(lengthmu,lengthka);

                for shelld = 1:shellb
                    la_begin = Shell_List(shelld,1);
                    la_end = Shell_List(shelld,2);
                    lengthla = la_end-la_begin+1;
                    %ad = Shells_List(d,2);
                    basis_d = basis{shelld};
                    Ld = basis_d(2);
                    %cd_data = pair_data{c,d};

                    if 1 %this later can be changed by a screening option
                        if (La == 0 && Lb == 0 && Lc == 0 && Ld == 0) %SS|SS integrals

                            gSSSSNValues = primitiveFactorsSSSS_3(basis_a,basis_b,basis_c,basis_d,Boys_Table,pair_data2,shella,shellb,shellc,shelld);
                            gmnkl = sum(gSSSSNValues(:,1));
                        else
                            gmnkl = shellOS(basis_a,basis_b,basis_c,basis_d,La,Lb,Lc,Ld,Boys_Table,pair_data2,shella,shellb,shellc,shelld);
                        endif
                        %Debug function to see if all the matrix is correctly covered
                        %gmnkl = ones(size(gmnkl));

                        intmu = mu_begin:mu_end;
                        intnu = nu_begin:nu_end;
                        intka = ka_begin:ka_end;
                        intla = la_begin:la_end;
                        if shellc == shelld
                            factorkl = 1 - 0.5*ones(lengthka,lengthla);
                            %factorkl = 0.5*ones(lengthka,lengthla);
                        else
                            factorkl = 1*ones(lengthka,lengthla);
                        end

                        if shellb == shelld
                            factornl = 1 - 0.5*ones(lengthnu,lengthla);
                        else
                            factornl = 1*ones(lengthnu,lengthla);
                        end

                        if shelld == shella
                            factorml = 1 - 0.5*ones(lengthmu,lengthla);
                        else
                            factorml = 1*ones(lengthmu,lengthla);
                        end

                        if shellb == shellc
                            factornk = 1 - 0.5*ones(lengthnu,lengthka);
                        else
                            factornk = 1*ones(lengthnu,lengthka);
                        end

                        %To check if gmmmm integrals have to be further divided by some factor
                        if (shella == shellb == shellc == shelld)
                            factormnkl = 0.5*ones(lengthmu,lengthnu,lengthka,lengthla);
                            %indD = [(1:lengthmu)' (1:lengthmu)' (1:lengthmu)' (1:lengthmu)'];
                            %factormnkl(indD) = 0.5*factormnkl(indD);
                            gmnkl = factormnkl.*gmnkl;
                        end
                        %Now it is important to treat matrices differently depending on their symmetry
                        %For example, matrices which lay on diagonals if permuted and added, their diagonal elements would count twice

                        %Contributions to Coulomb matrix
                        try
                            Dtemp = factorkl.*D(intka,intla);
                            Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                            J(intmu,intnu) = J(intmu,intnu) +  factormn.*squeeze(sum(sum(Jfactor*(Dres).*gmnkl,4),3));
                        catch
                            disp(J(intmu,intnu));
                            disp(D(intka,intla));
                            disp(gmnkl);
                        end
                        gklmn = permute(gmnkl,[3 4 1 2]);
                        try
                            Dtemp = factormn.*D(intmu,intnu);
                            Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                            J(intka,intla) = J(intka,intla) +  factorkl.*squeeze(sum(sum(Jfactor*(Dres).*gklmn,4),3));
                        catch
                            disp("J(intka,intla)");
                            disp(J(intka,intla));
                            disp("factorkl");
                            disp(factorkl);
                            disp("D(intmu,intnu)");
                            disp(D(intmu,intnu));
                            disp("gklmn");
                            disp(gklmn);
                            disp("Jkl");
                            disp(squeeze(sum(sum(Jfactor*(Dres).*gklmn,4),3)));
                        end
                        %Contributions to Exchange matrix
                        gmknl = permute(gmnkl,[1 3 2 4]);
                        Dtemp = factornl.*D(intnu,intla);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intmu,intka) = K(intmu,intka) -  factormk.*squeeze(sum(sum(Kfactor*(Dres).*gmknl,4),3)); %It may be necessary to permute this to gmknl, but I don't know yet                        Dtemp = D(intmu,intka);

                        gnmlk = permute(gmnkl,[2 1 4 3]);
                        gnlmk = permute(gnmlk,[1 3 2 4]);
                        Dtemp = factormk.*D(intmu,intka);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intnu,intla) = K(intnu,intla) -  factornl.*squeeze(sum(sum(Kfactor*(Dres).*gnlmk,4),3));

                        gmnlk = permute(gmnkl,[1 2 4 3]);
                        gmlnk = permute(gmnlk,[1 3 2 4]);
                        Dtemp = factornk.*D(intnu,intka);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intmu,intla) = K(intmu,intla) -  factorml.*squeeze(sum(sum(Kfactor*(Dres).*gmlnk,4),3));

                        gnmkl = permute(gmnkl,[2 1 3 4]);
                        gnkml = permute(gnmkl,[1 3 2 4]);
                        Dtemp = factorml.*D(intmu,intla);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intnu,intka) = K(intnu,intka) -  factornk.*squeeze(sum(sum(Kfactor*(Dres).*gnkml,4),3));
                    else
                        continue;
                    endif
                endfor
            else
                factormk = 1*ones(lengthmu,lengthka);
                for shelld = 1:shellc
                    la_begin = Shell_List(shelld,1);
                    la_end = Shell_List(shelld,2);
                    lengthla = la_end-la_begin+1;
                    %ad = Shell_List(d,2);
                    basis_d = basis{shelld};
                    Ld = basis_d(2);
                    %cd_data = pair_data{c,d};

                    if 1 %Later this can be changed for a screening option
                        if (La == 0 && Lb == 0 && Lc == 0 && Ld == 0) %SS|SS integrals

                            gSSSSNValues = primitiveFactorsSSSS_3(basis_a,basis_b,basis_c,basis_d,Boys_Table,pair_data2,shella,shellb,shellc,shelld);
                            gmnkl = sum(gSSSSNValues(:,1));
                        else
                            gmnkl = shellOS(basis_a,basis_b,basis_c,basis_d,La,Lb,Lc,Ld,Boys_Table,pair_data2,shella,shellb,shellc,shelld);
                        endif

                        %Debug function to see if all the matrix is correctly covered
                        %gmnkl = ones(size(gmnkl));

                        intmu = mu_begin:mu_end;
                        intnu = nu_begin:nu_end;
                        intka = ka_begin:ka_end;
                        intla = la_begin:la_end;

                        if shellc == shelld
                            factorkl = 1 - 0.5*ones(lengthka,lengthla);
                            %factorkl = 0.5*ones(lengthka,lengthla);
                        else
                            factorkl = 1*ones(lengthka,lengthla);
                        end

                        if shellb == shelld
                            factornl = 1 - 0.5*ones(lengthnu,lengthla);
                        else
                            factornl = 1*ones(lengthnu,lengthla);
                        end

                        if shelld == shella
                            factorml = 1 - 0.5*ones(lengthmu,lengthla);
                        else
                            factorml = 1*ones(lengthmu,lengthla);
                        end

                        if shellb == shellc
                            factornk = 1 - 0.5*ones(lengthnu,lengthka);
                        else
                            factornk = 1*ones(lengthnu,lengthka);
                        end

                        %To check if gmmmm integrals have to be further divided by some factor
                        if (shella == shellb == shellc == shelld)
                            factormnkl = 0.5*ones(lengthmu,lengthnu,lengthka,lengthla);
                            %indD = [(1:lengthmu)' (1:lengthmu)' (1:lengthmu)' (1:lengthmu)'];
                            %factormnkl(indD) = 0.5*factormnkl(indD);
                            gmnkl = factormnkl.*gmnkl;
                        end

                        %Now it is important to treat matrices differently depending on their symmetry
                        %For example, matrices which lay on diagonals if permuted and added, their diagonal elements would count twice
                        % gklmn = gmnkl;
                        % gmnlk = gmnkl;
                        % glkmn = gmnkl;
                        %Contributions to Coulomb matrix
                        try
                            Dtemp = factorkl.*D(intka,intla);
                            Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                            J(intmu,intnu) = J(intmu,intnu) +  factormn.*squeeze(sum(sum(Jfactor*(Dres).*gmnkl,4),3));
                        catch
                            disp(J(intmu,intnu));
                            disp(D(intka,intla));
                            disp(gmnkl);
                        end
                        gklmn = permute(gmnkl,[3 4 1 2]);
                        try
                            Dtemp = factormn.*D(intmu,intnu);
                            Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                            J(intka,intla) = J(intka,intla) +  factorkl.*squeeze(sum(sum(Jfactor*(Dres).*gklmn,4),3));
                        catch
                            disp("J(intka,intla)");
                            disp(J(intka,intla));
                            disp("factorkl");
                            disp(factorkl);
                            disp("D(intmu,intnu)");
                            disp(D(intmu,intnu));
                            disp("gklmn");
                            disp(gklmn);
                            disp("Jkl");
                            disp(squeeze(sum(sum(Jfactor*(Dres).*gklmn,4),3)));
                        end
                        %Contributions to Exchange matrix
                        gmknl = permute(gmnkl,[1 3 2 4]);
                        Dtemp = factornl.*D(intnu,intla);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intmu,intka) = K(intmu,intka) -  factormk.*squeeze(sum(sum(Kfactor*(Dres).*gmknl,4),3)); %It may be necessary to permute this to gmknl, but I don't know yet                        Dtemp = D(intmu,intka);

                        gnmlk = permute(gmnkl,[2 1 4 3]);
                        gnlmk = permute(gnmlk,[1 3 2 4]);
                        Dtemp = factormk.*D(intmu,intka);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intnu,intla) = K(intnu,intla) -  factornl.*squeeze(sum(sum(Kfactor*(Dres).*gnlmk,4),3));

                        gmnlk = permute(gmnkl,[1 2 4 3]);
                        gmlnk = permute(gmnlk,[1 3 2 4]);
                        Dtemp = factornk.*D(intnu,intka);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intmu,intla) = K(intmu,intla) -  factorml.*squeeze(sum(sum(Kfactor*(Dres).*gmlnk,4),3));

                        gnmkl = permute(gmnkl,[2 1 3 4]);
                        gnkml = permute(gnmkl,[1 3 2 4]);
                        Dtemp = factorml.*D(intmu,intla);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intnu,intka) = K(intnu,intka) -  factornk.*squeeze(sum(sum(Kfactor*(Dres).*gnkml,4),3));

                    else
                        continue;
                    endif
                endfor
            endif
        endfor
    endfor
endfor

%Because I only have the lowerdiagonal part
Jtril = tril(J); %to remove upper diagonal elements which occurred due to the use of shells.
J = tril(Jtril,-1)'+Jtril;

Ktril = tril(K); %to remove upper diagonal elements which occurred due to the use of shells.
K = tril(Ktril,-1)'+Ktril;

Coulomb = J;
Exchange = K;


end
