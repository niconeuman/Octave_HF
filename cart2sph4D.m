function Asph = cart2sph4D(A)


Aperm = permute(A,[3,4,1,2]);

%I transform the 4D array into a cell array of dimensions dim3 by dim4,
%where each element is the matrix corresponding to the first two dimensions
Acell = mat2cell(Aperm,size(Aperm,1),size(Aperm,2),ones(size(Aperm,3),1),ones(size(Aperm,4),1));

AsphPartial = cellfun(@cart2sph,Acell,"UniformOutput",0);

AsphPerm = cell2mat(AsphPartial);

Asph = permute(AsphPerm,[3,4,1,2]);

Acell = mat2cell(Aperm,size(Aperm,1),size(Aperm,2),ones(size(Aperm,3),1),ones(size(Aperm,4),1));

end
