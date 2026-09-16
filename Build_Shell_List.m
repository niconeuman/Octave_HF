function [Shell_List,Sph_Shell_List] = Build_Shell_List(basis)
%This function builds a shell list which has 3 columns,
%[mu_begin mu_end nb]
%It is a much smaller version of Build_Shells, but contains the same
%information in a non-redundant way.

nshells = size(basis,1);
Shell_List = zeros(nshells,3);
Sph_Shell_List = zeros(nshells,3);
mu_begin = 1;

sph_begin = 1;

for t = 1:nshells
  %L = 0, Length_mu = 1
  %L = 1, Length_mu = 3
  %L = 2, Length_mu = 5 Not 6!!!
    Length_sph = (2*basis{t}(2)+1); %basis{t}(2) = L; %spherical armonics
    Length_mu = (basis{t}(2)+1)*(basis{t}(2)+2)/2; %basis{t}(2) = L; %cartesian gto
    sph_end = sph_begin + Length_sph - 1;
    mu_end = mu_begin + Length_mu - 1;
    Sph_Shell_List(t,:) = [sph_begin sph_end t];
    Shell_List(t,:) = [mu_begin mu_end t];
    sph_begin = sph_end + 1;
    mu_begin = mu_end + 1;
end


end
