function [Da,Db] = uBuildDensity(Ca,Cb,Nela,Nelb)

nb = size(Ca,1);
Da = zeros(nb,nb);
Db = zeros(nb,nb);

for i = 1:Nela
    Da = Da + 1*Ca(:,i)*Ca(:,i)';
end

for i = 1:Nelb
    Db = Db + 1*Cb(:,i)*Cb(:,i)';
end
