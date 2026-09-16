function gabcd = Build_ERI_OS_2(basis,Shell_List,Boys_Table,pair_data2)

%This function will do what Build_ERI_OS did, but then it will contract the shell quartets using the density matrix
%to generate the

nshells = Shell_List(end,3);
Ncont = Shell_List(end,2);
gabcd = zeros(Ncont,Ncont,Ncont,Ncont);

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
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = gmnkl;
                        else
                            gmnkl = shellOS(basis_a,basis_b,basis_c,basis_d,La,Lb,Lc,Ld,Boys_Table,pair_data2,shella,shellb,shellc,shelld);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = gmnkl;
                        endif

                        intmu = mu_begin:mu_end;
                        intnu = nu_begin:nu_end;
                        intka = ka_begin:ka_end;
                        intla = la_begin:la_end;
                        gabcd(intnu,intmu,intka,intla) = permute(gabcd(intmu,intnu,intka,intla),[2 1 3 4]);
                        gabcd(intmu,intnu,intla,intka) = permute(gabcd(intmu,intnu,intka,intla),[1 2 4 3]);
                        gabcd(intnu,intmu,intla,intka) = permute(gabcd(intmu,intnu,intka,intla),[2 1 4 3]);
                        gabcd(intka,intla,intmu,intnu) = permute(gabcd(intmu,intnu,intka,intla),[3 4 1 2]);
                        gabcd(intla,intka,intmu,intnu) = permute(gabcd(intmu,intnu,intka,intla),[4 3 1 2]);
                        gabcd(intka,intla,intnu,intmu) = permute(gabcd(intmu,intnu,intka,intla),[3 4 2 1]);
                        gabcd(intla,intka,intnu,intmu) = permute(gabcd(intmu,intnu,intka,intla),[4 3 2 1]);
                        
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
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = gmnkl;
                        else
                            gmnkl = shellOS(basis_a,basis_b,basis_c,basis_d,La,Lb,Lc,Ld,Boys_Table,pair_data2,shella,shellb,shellc,shelld);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = gmnkl;
                        endif

                        intmu = mu_begin:mu_end;
                        intnu = nu_begin:nu_end;
                        intka = ka_begin:ka_end;
                        intla = la_begin:la_end;
                        gabcd(intnu,intmu,intka,intla) = permute(gabcd(intmu,intnu,intka,intla),[2 1 3 4]);
                        gabcd(intmu,intnu,intla,intka) = permute(gabcd(intmu,intnu,intka,intla),[1 2 4 3]);
                        gabcd(intnu,intmu,intla,intka) = permute(gabcd(intmu,intnu,intka,intla),[2 1 4 3]);
                        gabcd(intka,intla,intmu,intnu) = permute(gabcd(intmu,intnu,intka,intla),[3 4 1 2]);
                        gabcd(intla,intka,intmu,intnu) = permute(gabcd(intmu,intnu,intka,intla),[4 3 1 2]);
                        gabcd(intka,intla,intnu,intmu) = permute(gabcd(intmu,intnu,intka,intla),[3 4 2 1]);
                        gabcd(intla,intka,intnu,intmu) = permute(gabcd(intmu,intnu,intka,intla),[4 3 2 1]);


                    else
                        continue;
                    endif
                endfor
            endif
        endfor
    endfor
endfor

end
