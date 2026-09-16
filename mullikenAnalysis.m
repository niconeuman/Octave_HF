function [MullikenCharge,MullikenSpin] = mullikenAnalysis(S,Da,Db,Z,atomBasisList,Shell_List)

MullikenCharge = zeros(size(Z'));
MullikenSpin = zeros(size(Z'));

for iatom = 1:length(Z)
    basisStart = atomBasisList(iatom,1);
    basisEnd = atomBasisList(iatom,2);
    muStart  = Shell_List(basisStart,1);
    muEnd    = Shell_List(basisEnd,2);

    DaSBlock = (Da*S)(muStart:muEnd,muStart:muEnd);
    DbSBlock = (Db*S)(muStart:muEnd,muStart:muEnd);

    %disp([basisStart,basisEnd,muStart,muEnd]);

    DaBlock = Da(muStart:muEnd,muStart:muEnd);
    DbBlock = Db(muStart:muEnd,muStart:muEnd);
    SBlock  =  S(muStart:muEnd,muStart:muEnd);

    MullikenCharge(iatom) = Z(iatom)-sum(diag(DaSBlock))-sum(diag(DbSBlock));
    MullikenSpin(iatom) = sum(diag(DaSBlock))-sum(diag(DbSBlock));

end




end
