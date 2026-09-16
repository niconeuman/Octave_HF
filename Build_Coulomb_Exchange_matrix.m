function [G,Coulomb,Exchange] = Build_Coulomb_Exchange_matrix(D,gabcd)


nb = size(D,1);
K2 = zeros(nb,nb);

lengthgabcd = length(gabcd);
lengthgabcdSq = lengthgabcd*lengthgabcd;
gabcd2D = reshape(gabcd,[lengthgabcdSq,lengthgabcdSq]);

%This seems to always hold
%gabcd4D(i,j,k,l) == gabcd2D((j-1)*4+i,(l-1)*4+k)

ind = (1:lengthgabcdSq)'; %'
indD = [ind ind];

%This code works!
Dtemp = D;
Dtemp(indD) = 0.5*Dtemp(indD);
%gabcdPerm = permute(gabcd,[3,4,1,2]);
J2 = squeeze(2*sum(sum(Dtemp.*gabcd,2),1));

%This also works
gabcdPerm = permute(gabcd,[1,4,3,2]);
K2 = squeeze(-0.5*sum(sum(D.*gabcdPerm,2),1));
%disp('K2 by matrix broadcast method is');
%disp(K2(1:4,1:4));

%disp('J2 by matrix broadcast method is');
%disp(J2(1:4,1:4));

%J2 = zeros(Ncont,Ncont);

%for i = 1:nb
%    for j = 1:i
%        gabcd_ij = squeeze(gabcd(i,j,:,:));
%        gabcd_ij(indD) = 0.5*gabcd_ij(indD); %because all elements are multiplied by two, but diagonals should be multiplied by one

%        J2(i,j) = J2(i,j)+2*sum(sum(D.*gabcd_ij));

%        for k = 1:nb
%            for l = 1:k-1
%                %J2(i,j) = J2(i,j) + 2*D(k,l)*gabcd(i,j,k,l);
%                J2(i,j) = J2(i,j) + 2*D(k,l)*gabcd2D((j-1)*lengthgabcd+i,(l-1)*lengthgabcd+k);
%            end
%            l = k;
%                %J2(i,j) = J2(i,j) + D(k,l)*gabcd(i,j,k,l);
%                J2(i,j) = J2(i,j) + D(k,l)*gabcd2D((j-1)*lengthgabcd+i,(l-1)*lengthgabcd+k);
%        end
%    end
%end

%for i = 1:nb
%    for j = 1:i
%        J2(j,i) = J2(i,j);
%    end
%end

%disp('J2 by loop over i,j and matrix operation within is');
%disp(J2(1:4,1:4));

%K2 = zeros(Ncont,Ncont);
%for i = 1:nb
%    for l = 1:i
%        gabcd_il = squeeze(gabcd(i,:,:,l));
%        K2(i,l) = K2(i,l)-0.5*sum(sum(D.*gabcd_il));
%%        for k = 1:nb
%%            for j = 1:nb
%%                %K2(i,l) = K2(i,l) - 0.5*D(k,j)*gabcd(i,j,k,l);
% %               K2(i,l) = K2(i,l) - 0.5*D(k,j)*gabcd2D((j-1)*lengthgabcd+i,(l-1)*lengthgabcd+k);
%%            end
%%       end
%    end
%end

%disp('K2 by loop over i,j and matrix operation within is');
%disp(K2(1:4,1:4));

%{
for i = 1:nb
%    mu_begin = Shell_List(i,1);
%    mu_end = Shell_List(i,2);
    for l = 1:i
%        la_begin = Shell_List(l,1);
%        la_end = Shell_List(l,2);


        for k = 1:nb
%            ka_begin = Shell_List(k,1);
%            ka_end = Shell_List(k,2);
            for j = 1:nb
%                nu_begin = Shell_List(j,1);
%                nu_end = Shell_List(j,2);
                %K2(i,l) = K2(i,l) - 0.5*D(k,j)*gabcd(i,j,k,l);
                K2(i,l) = K2(i,l) - 0.5*D(k,j)*gabcd2D((j-1)*lengthgabcd+i,(l-1)*lengthgabcd+k);
            end
        end
    end
end
%}


%G = J2 + K2;

%for i = 1:nb
%    for j = 1:i
%        K2(j,i) = K2(i,j);
%    end
%end
Coulomb = J2;
Exchange = K2;
G = J2 + K2;

J = zeros(nb,nb);
K = zeros(nb,nb);

for i = 1:nb
    for j = 1:nb
        for k = 1:nb
            for l = 1:nb
                J(i,j) = J(i,j) + D(k,l)*gabcd(i,j,k,l);
                %K(i,l) = K(i,l) - 0.5*D(k,j)*gabcd(i,j,k,l);
                K(i,j) = K(i,j) - 0.5*D(k,l)*gabcd(i,k,j,l);
            end
        end
    end
end

Coulomb = J;
Exchange = K;
G = J + K;


endfunction
