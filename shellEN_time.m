function [VEN,tInitVectors,tLoopVectors] = shellEN_time(basis_a,basis_b,La,Lb,pair_data2,AL,Z,Boys_Table)


[RPAValues,RPBValues,RPCValues,RAB,pValues,VssNValues,tInitVectors,tLoopVectors] = primitiveFactorsNuc_time(basis_a,basis_b,La,Lb,Boys_Table,AL,Z);

if (Lb > La)
    Vswap = shellEN(basis_b,basis_a,Lb,La,pair_data2,AL,Z,Boys_Table);
    VEN = permute(Vswap,[2 1]);

elseif (La == 0 && Lb == 0) %[s|s]

%12may2019. I have to check that primitiveFactorsNuc works fine, comparing it to primitiveFactors3, etc.
    %[RPAValues,RPBValues,RPCValues,RAB,pValues,VssNValues] = primitiveFactorsNuc(basis_a,basis_b,La,Lb,Boys_Table,AL,Z);

    Vout = ENss(RPAValues,RPBValues,RPCValues,pValues,VssNValues);
    VEN = Vout;

elseif (La == 1 && Lb == 0) %[p|s]

    %[RPAValues,RPBValues,RPCValues,RAB,pValues,VssNValues] = primitiveFactorsNuc(basis_a,basis_b,La,Lb,Boys_Table,AL,Z);
    Vout = ENps(RPAValues,RPBValues,RPCValues,pValues,VssNValues);
    VEN = Vout;

elseif (La == 2 && Lb == 0) %[d|s]

    %[RPAValues,RPBValues,RPCValues,RAB,pValues,VssNValues] = primitiveFactorsNuc(basis_a,basis_b,La,Lb,Boys_Table,AL,Z);
    Vout = ENds(RPAValues,RPBValues,RPCValues,pValues,VssNValues);
    VEN = Vout;

elseif (La == 0 && Lb == 1) %[s|p]

    %[RPAValues,RPBValues,RPCValues,RAB,pValues,VssNValues] = primitiveFactorsNuc(basis_a,basis_b,La,Lb,Boys_Table,AL,Z);
    Vout = ENsp(RPAValues,RPBValues,RPCValues,pValues,VssNValues);
    VEN = Vout;

elseif (La == 0 && Lb == 2) %[s|d]

    %[RPAValues,RPBValues,RPCValues,RAB,pValues,VssNValues] = primitiveFactorsNuc(basis_a,basis_b,La,Lb,Boys_Table,AL,Z);
    Vout = ENsd(RPAValues,RPBValues,RPCValues,pValues,VssNValues);
    VEN = Vout;

elseif (La == 1 && Lb == 1) %[p|p]

    %[RPAValues,RPBValues,RPCValues,RAB,pValues,VssNValues] = primitiveFactorsNuc(basis_a,basis_b,La,Lb,Boys_Table,AL,Z);
    Vout = ENpp(RPAValues,RPBValues,RPCValues,pValues,VssNValues);
    VEN = Vout;

elseif (La == 2 && Lb == 1) %[d|p]

    %[RPAValues,RPBValues,RPCValues,RAB,pValues,VssNValues] = primitiveFactorsNuc(basis_a,basis_b,La,Lb,Boys_Table,AL,Z);
    Vout = ENdp(RPAValues,RPBValues,RPCValues,pValues,VssNValues);
    VEN = Vout;

elseif (La == 2 && Lb == 2) %[d|d]

    %[RPAValues,RPBValues,RPCValues,RAB,pValues,VssNValues] = primitiveFactorsNuc(basis_a,basis_b,La,Lb,Boys_Table,AL,Z);
    Vout = ENdd(RPAValues,RPBValues,RPCValues,pValues,VssNValues);
    VEN = Vout;

end

end
