function gabcd = Build_ERI_OS(basis,Shell_List,Boys_Table,pair_data2)
%2 jul 2017. This is a function similar to Build_ERI, but which uses 4
%loops because it uses the data of Shell_List (an nb x 3 matrix), instead
%of Shells, (an nb^4 x 12) matrix.

%Size of the contracted basis list (size of gabcd matrix)
Ncont = Shell_List(end,2);
nb = size(basis,1);
gabcd = zeros(Ncont,Ncont,Ncont,Ncont);
shellThr = 1e-7; %Seems to be a reasonable value

for a = 1:nb
    mu_begin = Shell_List(a,1);
    mu_end = Shell_List(a,2);
    %aa = Shells_List(a,2);
    basis_a = basis{a};
    for b = 1:a
        nu_begin = Shell_List(b,1);
        nu_end = Shell_List(b,2);
        %ab = Shells_List(b,2);
        basis_b = basis{b};
        %ab_data = pair_data{a,b};
        %This part is a first attempt at prescreening
        KabMax = max(pair_data2{a,b}(:,5));
        for c = 1:a
            ka_begin = Shell_List(c,1);
            ka_end = Shell_List(c,2);
            %ac = Shells_List(c,2);
            basis_c = basis{c};
            if (c == a)
                for d = 1:b
                    la_begin = Shell_List(d,1);
                    la_end = Shell_List(d,2);
                    %ad = Shells_List(d,2);
                    basis_d = basis{d};
                    %cd_data = pair_data{c,d};
                    KcdMax = max(pair_data2{c,d}(:,5));
                    if abs(KabMax*KcdMax)>shellThr
                        if (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %SS|SS integrals

                            gSSSSNValues = primitiveFactorsSSSS_3(basis_a,basis_b,basis_c,basis_d,Boys_Table,pair_data2,a,b,c,d);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = sum(gSSSSNValues(:,1));

                            % elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %LS|SS integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                            % elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %SL|SS integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                            % elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %SS|LS integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                            % elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %SS|SL integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                            % elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %LaS|LcS integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                            % elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %SLb|LcS integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                            % elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %LaS|SLd integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                            % elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %SLb|SLd integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                            % elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %LaLb|SS integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                            % elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %SS|LcLd integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                            % elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %LaLb|LcS integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                            % elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %LaLb|SLd integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                            % elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %LaS|LcLd integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                            % elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %SLb|LcLd integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                            % elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %LaLb|LcLd integrals
                            %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        else
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
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
                for d = 1:c
                    la_begin = Shell_List(d,1);
                    la_end = Shell_List(d,2);
                    %ad = Shell_List(d,2);
                    basis_d = basis{d};
                    %cd_data = pair_data{c,d};
                    KcdMax = max(pair_data2{c,d}(:,5));
                    if abs(KabMax*KcdMax)>shellThr
                        if (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %SS|SS integrals

                            gSSSSNValues = primitiveFactorsSSSS_3(basis_a,basis_b,basis_c,basis_d,Boys_Table,pair_data2,a,b,c,d);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = sum(gSSSSNValues(:,1));

                        % elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %LS|SS integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        % elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %SL|SS integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        % elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %SS|LS integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        % elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %SS|SL integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        % elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %LaS|LcS integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        % elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %SLb|LcS integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        % elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %LaS|SLd integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        % elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %SLb|SLd integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        % elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %LaLb|SS integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        % elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %SS|LcLd integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        % elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %LaLb|LcS integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        % elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %LaLb|SLd integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        % elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %LaS|LcLd integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        % elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %SLb|LcLd integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        % elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %LaLb|LcLd integrals
                        %     gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
                        else
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,pair_data2,a,b,c,d);
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
