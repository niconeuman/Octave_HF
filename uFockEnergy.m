function out = uFockEnergy(Da,Db,Fa,Fb,H0)
nb = size(Da,1);
DT = Da+Db;
E = 0;
for n = 1:nb
    for m = 1:nb
        E = E+0.5*(DT(n,m)*H0(n,m)+Da(n,m)*Fa(n,m)+Db(n,m)*Fb(n,m)); %The factor multiplying D is 1, not 1/2. Checked from D. Crawford's page and James Johns program
    end
end
out = E;
