function [Coulomb,Exchange] = Build_Coulomb_Exchange_2(Gprev,D,Dm1,basis,pair_data2,Shell_List,Boys_Table)

%This function will do what Build_ERI_OS did, but then it will contract the shell quartets using the density matrix
%to generate the

nshells = Shell_List(end,3);
Ncont = Shell_List(end,2);
J = zeros(Ncont,Ncont);
K = zeros(Ncont,Ncont);

Jfactor = 1;
Kfactor = 0.5;
%Debug function to see if all the matrix is correctly covered
%D = ones(size(D));

%If I use the following conditions
% symFactor = 1;
% if (shella == shellb)
%    symFactor = symFactor*0.5;
% end
%
% if (shellc == shelld)
%    symFactor = symFactor*0.5;
% end
%
% if ((shella == shellc) && (shellb == shelld))
%    symFactor = symFactor*0.5;
% end

%all my J matrix elements are the same if D = ones(3D) and gmnkl = ones(4D)

%If I use the following conditions

% symFactorKmknl = 1;
% if (shella == shellc)
%  symFactorKmknl = symFactorKmknl*0.5;
% end
% if (shellb == shelld)
%  symFactorKmknl = symFactorKmknl*0.5;
% end
%
% symFactorKmlnk = 1;
% if (shella == shelld)
%  symFactorKmlnk = symFactorKmlnk*0.5;
% end
% if (shellb == shellc)
%  symFactorKmlnk = symFactorKmlnk*0.5;
% end
%
% if ((shella == shellb) && (shellc == shelld))
%  symFactorKmlnk = symFactorKmlnk*0.5;
%  #symFactorKmknl = symFactorKmknl*0.5;
%  %symFactorK = 0.25;
% end
%
% if ((shella == shellc) && (shellb == shelld))
%  #symFactorKmlnk = symFactorKmlnk*0.5;
%  symFactorKmknl = symFactorKmknl*0.5;
%  %symFactorK = 0.25;
% end

%together with symFactor = 1; after the coulomb equations,
%I get approximately half the K values in the diagonal and different values in ps blocks or pp blocks than in ss blocks
%Commenting the blocks with only one if conditions doesn't fix it (symFactor = 1 still)
%Now commenting symFactor = 1, while keeping only ((shella == shellb) && (shellc == shelld))
%and if ((shella == shellc) && (shellb == shelld)) conditions, gives me close values to half the Coulomb, but still not all of them equal
%If I uncomment the following conditions
% if (shella == shellc)
%    symFactorKmknl = symFactorKmknl*0.5;
% end
% if (shellb == shelld)
%    symFactorKmknl = symFactorKmknl*0.5;
% end
%
%  symFactorKmlnk = 1;
% if (shella == shelld)
%    symFactorKmlnk = symFactorKmlnk*0.5;
% end
% if (shellb == shellc)
%    symFactorKmlnk = symFactorKmlnk*0.5;
% end
%And comment the double conditions, with symFactor = 1 commented,
%I get all equal elements except for non diag elements which are now half (including diagonal pp blocks)
%These elements are lower than half Coulomb
%Setting again with the same conditions symFactor = 1, gives me worse results. So I will keep that commented for now.

%Now I uncomment the following double condition
% if ((shella == shellb) && (shellc == shelld))
%    symFactorKmlnk = symFactorKmlnk*0.5;
%    symFactorKmknl = symFactorKmknl*0.5;
% end

%This only changed the pp blocks
%Uncommenting everything, but with symFactor = 1 still commented, doesn't fix it
%If I uncomment it, I get larger K values than half Coulomb. But perhaps I can now optimize the other conditions


%Exactly half the Coulomb values on all the matrix were obtained when only the same symFactor that was considered for Coulomb integrals was added to
%Exchange integrals. No specific symFactorK was needed (of course I'm explicitly doing the eight-fold contributions)

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

        for shellc = 1:shella
            ka_begin = Shell_List(shellc,1);
            ka_end = Shell_List(shellc,2);
            lengthka = ka_end-ka_begin+1;
            %ac = Shells_List(c,2);
            basis_c = basis{shellc};
            Lc = basis_c(2);
            if (shellc == shella)


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

                        symFactor = 1;
                        if (shella == shellb)
                           symFactor = symFactor*0.5;
                        end

                        if (shellc == shelld)
                           symFactor = symFactor*0.5;
                        end


                        if ((shella == shellc) && (shellb == shelld))
                           symFactor = symFactor*0.5;
                           %symFactorK = 0.25;
                        end

                        symFactorKmknl = 1;
                        symFactorKmlnk = 1;
                       % if (shella == shellc)
                       %    symFactorKmknl = symFactorKmknl*0.5;
                       %    symFactorKmlnk = symFactorKmlnk*0.5;
                       % end
                       % if (shellb == shelld)
                       %    symFactorKmknl = symFactorKmknl*0.5;
                       %    symFactorKmlnk = symFactorKmlnk*0.5;
                       % end

                       %  symFactorKmlnk = 1;
                       % if (shella == shelld)
                       %    symFactorKmlnk = symFactorKmlnk*0.5;
                       % end
                       % if (shellb == shellc)
                       %    symFactorKmlnk = symFactorKmlnk*0.5;
                       % end

                       % if ((shella == shellb) && (shellc == shelld))
                       %    symFactorKmlnk = symFactorKmlnk*0.5;
                       %    symFactorKmknl = symFactorKmknl*0.5;
                       %    %symFactorK = 0.25;
                       % end
                       %
                       % if ((shella == shellc) && (shellb == shelld))
                       %    #symFactorKmlnk = symFactorKmlnk*0.5;
                       %    symFactorKmknl = symFactorKmknl*0.5;
                       %    %symFactorK = 0.25;
                       % end

                        intmu = mu_begin:mu_end;
                        intnu = nu_begin:nu_end;
                        intka = ka_begin:ka_end;
                        intla = la_begin:la_end;

                        %This code may be slow but I hope it will be correct
                        gnmkl = permute(gmnkl,[2 1 3 4]);
                        gmnlk = permute(gmnkl,[1 2 4 3]);
                        gnmlk = permute(gmnkl,[2 1 4 3]);
                        gklmn = permute(gmnkl,[3 4 1 2]);
                        glkmn = permute(gmnkl,[4 3 1 2]);
                        gklnm = permute(gmnkl,[3 4 2 1]);
                        glknm = permute(gmnkl,[4 3 2 1]);

                        %Coulomb contribution, not making use of symmetry. Some of them will be counted several times
                        %mn kl
                        Dtemp = D(intka,intla);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        J(intmu,intnu) = J(intmu,intnu) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*gmnkl,4),3));

                        %nm kl

                           Dtemp = D(intka,intla);
                           Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                           J(intnu,intmu) = J(intnu,intmu) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*gnmkl,4),3));


                        %mn lk

                           Dtemp = D(intla,intka);
                           Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                           J(intmu,intnu) = J(intmu,intnu) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*gmnlk,4),3));


                        %nm lk

                           Dtemp = D(intla,intka);
                           Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                           J(intnu,intmu) = J(intnu,intmu) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*gnmlk,4),3));


                        %kl mn

                           Dtemp = D(intmu,intnu);
                           Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                           J(intka,intla) = J(intka,intla) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*gklmn,4),3));

                        %lk mn

                           Dtemp = D(intmu,intnu);
                           Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                           J(intla,intka) = J(intla,intka) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*glkmn,4),3));

                        %kl nm

                           Dtemp = D(intnu,intmu);
                           Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                           J(intka,intla) = J(intka,intla) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*gklnm,4),3));

                        %lk nm

                           Dtemp = D(intnu,intmu);
                           Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                           J(intla,intka) = J(intla,intka) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*glknm,4),3));


                        %Exchange contribution, all 2nd and 3rd indices have to be exchanged
                        gmknl = permute(gmnkl,[1 3 2 4]);

                        %Now 8-fold symmetry
                        gmlnk = permute(gmknl,[1 4 3 2]);
                        gnkml = permute(gmknl,[3 2 1 4]);
                        gnlmk = permute(gmknl,[3 4 1 2]);

                        % gkmnl = permute(gmknl,[2 1 3 4]);
                        % gmkln = permute(gmknl,[1 2 4 3]);
                         gkmln = permute(gmknl,[2 1 4 3]);
                         glmkn = permute(gmknl,[4 1 2 3]);
                         gknlm = permute(gmknl,[2 3 4 1]);
                         glnkm = permute(gmknl,[4 3 2 1]);
                        % gnlmk = permute(gmknl,[3 4 1 2]);
                        % glnmk = permute(gmknl,[4 3 1 2]);
                        % gnlkm = permute(gmknl,[3 4 2 1]);
                        % glnkm = permute(gmknl,[4 3 2 1]);
                        %symFactor = 1;
                        %mk nl
                        Dtemp = D(intnu,intla);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intmu,intka) = K(intmu,intka) -  symFactorKmknl*symFactor*squeeze(sum(sum(Kfactor*(Dres).*gmknl,4),3));

                        %ml nk
                        Dtemp = D(intnu,intka);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intmu,intla) = K(intmu,intla) -  symFactorKmlnk*symFactor*squeeze(sum(sum(Kfactor*(Dres).*gmlnk,4),3));

                        %nk ml
                        Dtemp = D(intmu,intla);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intnu,intka) = K(intnu,intka) -  symFactorKmlnk*symFactor*squeeze(sum(sum(Kfactor*(Dres).*gnkml,4),3));

                        %nl mk
                        Dtemp = D(intmu,intka);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intnu,intla) = K(intnu,intla) -  symFactorKmknl*symFactor*squeeze(sum(sum(Kfactor*(Dres).*gnlmk,4),3));

                        %km ln (done)
                        Dtemp = D(intla,intnu);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intka,intmu) = K(intka,intmu) -  symFactorKmknl*symFactor*squeeze(sum(sum(Kfactor*(Dres).*gkmln,4),3));

                        %lm kn (done)
                        Dtemp = D(intka,intnu);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intla,intmu) = K(intla,intmu) -  symFactorKmlnk*symFactor*squeeze(sum(sum(Kfactor*(Dres).*glmkn,4),3));

                        %kn lm (done)
                        Dtemp = D(intla,intmu);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intka,intnu) = K(intka,intnu) -  symFactorKmlnk*symFactor*squeeze(sum(sum(Kfactor*(Dres).*gknlm,4),3));

                        %ln km (done)
                        Dtemp = D(intka,intmu);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intla,intnu) = K(intla,intnu) -  symFactorKmknl*symFactor*squeeze(sum(sum(Kfactor*(Dres).*glnkm,4),3));

                        % %mk nl
                        % Dtemp = D(intnu,intla);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intmu,intka) = K(intmu,intka) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*gmknl,4),3));
                        %
                        % %km nl
                        % Dtemp = D(intnu,intla);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intka,intmu) = K(intka,intmu) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*gkmnl,4),3));
                        %
                        % %mk ln
                        % Dtemp = D(intla,intnu);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intmu,intka) = K(intmu,intka) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*gmkln,4),3));
                        %
                        % %km ln
                        % Dtemp = D(intla,intnu);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intka,intmu) = K(intka,intmu) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*gkmln,4),3));
                        %
                        % %nl mk
                        % Dtemp = D(intmu,intka);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intnu,intla) = K(intnu,intla) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*gnlmk,4),3));
                        %
                        % %ln mk
                        % Dtemp = D(intmu,intka);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intla,intnu) = K(intla,intnu) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*glnmk,4),3));
                        %
                        % %nl km
                        % Dtemp = D(intka,intmu);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intnu,intla) = K(intnu,intla) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*gnlkm,4),3));
                        %
                        % %ln km
                        % Dtemp = D(intka,intmu);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intla,intnu) = K(intla,intnu) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*glnkm,4),3));

                    else
                        continue;
                    endif
                endfor
            else

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

                        symFactor = 1;
                        if (shella == shellb)
                           symFactor = symFactor*0.5;
                        end

                        if (shellc == shelld)
                           symFactor = symFactor*0.5;
                        end

                        symFactorK = 1;
                        if ((shella == shellc) && (shellb == shelld))
                           symFactor = symFactor*0.5;
                           %symFactorK = 0.25;
                        end

                         symFactorKmknl = 1;
                         symFactorKmlnk = 1;
                        % if (shella == shellc)
                        %    symFactorKmknl = symFactorKmknl*0.5;
                        %    symFactorKmlnk = symFactorKmlnk*0.5;
                        % end
                        % if (shellb == shelld)
                        %    symFactorKmknl = symFactorKmknl*0.5;
                        %    symFactorKmlnk = symFactorKmlnk*0.5;
                        % end

                        %  symFactorKmlnk = 1;
                        % if (shella == shelld)
                        %    symFactorKmlnk = symFactorKmlnk*0.5;
                        % end
                        % if (shellb == shellc)
                        %    symFactorKmlnk = symFactorKmlnk*0.5;
                        % end

                        % if ((shella == shellb) && (shellc == shelld))
                        %    symFactorKmlnk = symFactorKmlnk*0.5;
                        %    symFactorKmknl = symFactorKmknl*0.5;
                        % end
                        %
                        % if ((shella == shellc) && (shellb == shelld))
                        %    #symFactorKmlnk = symFactorKmlnk*0.5;
                        %    symFactorKmknl = symFactorKmknl*0.5;
                        %    %symFactorK = 0.25;
                        % end

                        intmu = mu_begin:mu_end;
                        intnu = nu_begin:nu_end;
                        intka = ka_begin:ka_end;
                        intla = la_begin:la_end;

                        %This code may be slow but I hope it will be correct
                        gnmkl = permute(gmnkl,[2 1 3 4]);
                        gmnlk = permute(gmnkl,[1 2 4 3]);
                        gnmlk = permute(gmnkl,[2 1 4 3]);
                        gklmn = permute(gmnkl,[3 4 1 2]);
                        glkmn = permute(gmnkl,[4 3 1 2]);
                        gklnm = permute(gmnkl,[3 4 2 1]);
                        glknm = permute(gmnkl,[4 3 2 1]);

                        %Coulomb contribution, not making use of symmetry. Some of them will be counted several times
                        %mn kl
                        Dtemp = D(intka,intla);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        J(intmu,intnu) = J(intmu,intnu) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*gmnkl,4),3));

                        %nm kl

                           Dtemp = D(intka,intla);
                           Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                           J(intnu,intmu) = J(intnu,intmu) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*gnmkl,4),3));


                        %mn lk

                           Dtemp = D(intla,intka);
                           Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                           J(intmu,intnu) = J(intmu,intnu) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*gmnlk,4),3));


                        %nm lk

                           Dtemp = D(intla,intka);
                           Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                           J(intnu,intmu) = J(intnu,intmu) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*gnmlk,4),3));


                        %kl mn

                           Dtemp = D(intmu,intnu);
                           Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                           J(intka,intla) = J(intka,intla) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*gklmn,4),3));

                        %lk mn

                           Dtemp = D(intmu,intnu);
                           Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                           J(intla,intka) = J(intla,intka) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*glkmn,4),3));

                        %kl nm

                           Dtemp = D(intnu,intmu);
                           Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                           J(intka,intla) = J(intka,intla) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*gklnm,4),3));

                        %lk nm

                           Dtemp = D(intnu,intmu);
                           Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                           J(intla,intka) = J(intla,intka) +  symFactor*squeeze(sum(sum(Jfactor*(Dres).*glknm,4),3));


                        %Exchange contribution, all 2nd and 3rd indices have to be exchanged
                        gmknl = permute(gmnkl,[1 3 2 4]);

                        %Now 8-fold symmetry
                        gmlnk = permute(gmknl,[1 4 3 2]);
                        gnkml = permute(gmknl,[3 2 1 4]);
                        gnlmk = permute(gmknl,[3 4 1 2]);

                        gkmln = permute(gmknl,[2 1 4 3]);
                        glmkn = permute(gmknl,[4 1 2 3]);
                        gknlm = permute(gmknl,[2 3 4 1]);
                        glnkm = permute(gmknl,[4 3 2 1]);
                        % gkmnl = permute(gmknl,[2 1 3 4]);
                        % gmkln = permute(gmknl,[1 2 4 3]);
                        % gkmln = permute(gmknl,[2 1 4 3]);
                        % gnlmk = permute(gmknl,[3 4 1 2]);
                        % glnmk = permute(gmknl,[4 3 1 2]);
                        % gnlkm = permute(gmknl,[3 4 2 1]);
                        % glnkm = permute(gmknl,[4 3 2 1]);

                        %symFactor = 1;
                        %mk nl
                        Dtemp = D(intnu,intla);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intmu,intka) = K(intmu,intka) -  symFactorKmknl*symFactor*squeeze(sum(sum(Kfactor*(Dres).*gmknl,4),3));

                        %ml nk
                        Dtemp = D(intnu,intka);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intmu,intla) = K(intmu,intla) -  symFactorKmlnk*symFactor*squeeze(sum(sum(Kfactor*(Dres).*gmlnk,4),3));

                        %nk ml
                        Dtemp = D(intmu,intla);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intnu,intka) = K(intnu,intka) -  symFactorKmlnk*symFactor*squeeze(sum(sum(Kfactor*(Dres).*gnkml,4),3));

                        %nl mk
                        Dtemp = D(intmu,intka);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intnu,intla) = K(intnu,intla) -  symFactorKmknl*symFactor*squeeze(sum(sum(Kfactor*(Dres).*gnlmk,4),3));

                        %km ln (done)
                        Dtemp = D(intla,intnu);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intka,intmu) = K(intka,intmu) -  symFactorKmknl*symFactor*squeeze(sum(sum(Kfactor*(Dres).*gkmln,4),3));

                        %lm kn (done)
                        Dtemp = D(intka,intnu);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intla,intmu) = K(intla,intmu) -  symFactorKmlnk*symFactor*squeeze(sum(sum(Kfactor*(Dres).*glmkn,4),3));

                        %kn lm (done)
                        Dtemp = D(intla,intmu);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intka,intnu) = K(intka,intnu) -  symFactorKmlnk*symFactor*squeeze(sum(sum(Kfactor*(Dres).*gknlm,4),3));

                        %ln km (done)
                        Dtemp = D(intka,intmu);
                        Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        K(intla,intnu) = K(intla,intnu) -  symFactorKmknl*symFactor*squeeze(sum(sum(Kfactor*(Dres).*glnkm,4),3));
                        % %km nl
                        % Dtemp = D(intnu,intla);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intka,intmu) = K(intka,intmu) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*gkmnl,4),3));

                        % %mk ln
                        % Dtemp = D(intla,intnu);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intmu,intka) = K(intmu,intka) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*gmkln,4),3));

                        % %km ln
                        % Dtemp = D(intla,intnu);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intka,intmu) = K(intka,intmu) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*gkmln,4),3));
                        %
                        % %nl mk
                        % Dtemp = D(intmu,intka);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intnu,intla) = K(intnu,intla) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*gnlmk,4),3));
                        %
                        % %ln mk
                        % Dtemp = D(intmu,intka);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intla,intnu) = K(intla,intnu) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*glnmk,4),3));
                        %
                        % %nl km
                        % Dtemp = D(intka,intmu);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intnu,intla) = K(intnu,intla) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*gnlkm,4),3));
                        %
                        % %ln km
                        % Dtemp = D(intka,intmu);
                        % Dres = reshape(Dtemp,[1,1,size(Dtemp)]);
                        % K(intla,intnu) = K(intla,intnu) -  symFactor*squeeze(sum(sum(Kfactor*(Dres).*glnkm,4),3));

                    else
                        continue;
                    endif
                endfor
            endif
        endfor
    endfor
endfor

%Because I only have the lowerdiagonal part
%Jtril = tril(J); %to remove upper diagonal elements which occurred due to the use of shells.
%J = tril(Jtril,-1)'+Jtril;

%Ktril = tril(K); %to remove upper diagonal elements which occurred due to the use of shells.
%K = tril(Ktril,-1)'+Ktril;

Coulomb = J;
Exchange = K;


end
