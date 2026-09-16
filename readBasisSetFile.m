function [coeffMatrices,matrixTypes] = readBasisSetFile(fileName,elementName)

%This function reads a text file (fileName = 'example.ext') which contains basis set information in
%nw-chem format, as the following example, and returns the matrices in a
%cell array, for elementName ('H', 'He', 'Cl', etc)
%matrix types are 'S' 'SP' 'P' 'D' etc
% #----------------------------------------------------------------------
% # Basis Set Exchange
% # Version v0.8.2
% # https://www.basissetexchange.org
% #----------------------------------------------------------------------
% #   Basis set: STO-3G
% # Description: STO-3G Minimal Basis (3 functions/AO)
% #        Role: orbital
% #     Version: 1  (Data from Gaussian09)
% #----------------------------------------------------------------------
%
%
% BASIS "ao basis" PRINT
% #BASIS SET: (9s,6p) -> [3s,2p]
% Na    S
%       0.2507724300E+03       0.1543289673E+00
%       0.4567851117E+02       0.5353281423E+00
%       0.1236238776E+02       0.4446345422E+00
% Na    SP
%       0.1204019274E+02      -0.9996722919E-01       0.1559162750E+00
%       0.2797881859E+01       0.3995128261E+00       0.6076837186E+00
%       0.9099580170E+00       0.7001154689E+00       0.3919573931E+00
% Na    SP
%       0.1478740622E+01      -0.2196203690E+00       0.1058760429E-01
%       0.4125648801E+00       0.2255954336E+00       0.5951670053E+00
%       0.1614750979E+00       0.9003984260E+00       0.4620010120E+00
% #BASIS SET: (9s,6p) -> [3s,2p]
% P    S
%       0.4683656378E+03       0.1543289673E+00

fileContentsStr = fileread(fileName);

%if elementName is one letter long, like S or P, there can be confusions
%with the name of the matrices S and SP. So if that is the case we padd it
%with a following whitespace.
%I actually need 2 whitespaces, because some words end up with S for
%example

if length(elementName) == 1
    elementName = [elementName, '  '];
end


firstCharacterPositions = strfind(fileContentsStr,elementName);
%disp('firstCharacterPositions');
%disp(firstCharacterPositions);
endBlock = strfind(fileContentsStr,'#BASIS');
%disp(endBlock);
nMatrices = length(firstCharacterPositions);

coeffMatrices = cell(nMatrices,1);
matrixTypes = cell(nMatrices,1);

for k = 1:nMatrices

    if k < nMatrices
        partialStr = fileContentsStr(firstCharacterPositions(k):firstCharacterPositions(k+1)-1); %the -1 is because I want to end on the character before the next elementName
        %disp('partialStr');
        %disp(partialStr);
    else
        %I need to find the closest endBlock element that is higher than
        %firsCharacterPositions(k);
        nextBlock = endBlock > firstCharacterPositions(k);
        nextBlockPos = endBlock(nextBlock);
        %disp(nextBlockPos(1));

        partialStr = fileContentsStr(firstCharacterPositions(k):nextBlockPos(1)-1);
    end

    %disp(partialStr);
    %partialLines = splitlines(partialStr); %MATLAB version
    partialLines = strsplit(partialStr,'\n');
    %disp('partialLines');
    %disp(partialLines);
    %nLines = size(partialLines,1); %This worked on Matlab
    nLines = size(partialLines,2); %This works in Octave
    %disp('nLines');
    %disp(nLines);
    tempMatrix = [];

    titleLine = partialLines{1}; %This lines contain the matrixType information

    try %Matlab syntax
        if contains(titleLine,' S')
            matrixTypes{k} = 'S';
        elseif contains(titleLine,' P')
            matrixTypes{k} = 'P';
        elseif contains(titleLine,' SP')
            matrixTypes{k} = 'SP';
        elseif contains(titleLine,' D')
            matrixTypes{k} = 'D';
        end
    catch
        if ~isempty(strfind(titleLine,' S'))
            matrixTypes{k} = 'S';
        elseif ~isempty(strfind(titleLine,' P'))
            matrixTypes{k} = 'P';
        elseif ~isempty(strfind(titleLine,' SP'))
            matrixTypes{k} = 'SP';
        elseif ~isempty(strfind(titleLine,' D'))
            matrixTypes{k} = 'D';
        end
    end
    %disp(partialLines{2});
    for p = 2:nLines
        %disp('p is');
        %disp(p);
        %disp('partialLines{p}');
        %disp(partialLines{p});
        temp = partialLines{p};  %Matlab syntax
        %temp = partialLines(p);  %Matlab syntax
        %disp(temp);
        tempMatrix = vertcat(tempMatrix,str2num(temp));
    end
    %str2num(temp{1});

    coeffMatrices{k} = tempMatrix;

end













end
