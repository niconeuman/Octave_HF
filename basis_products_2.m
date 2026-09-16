function [pair_data,RABdata,pOverlap] = basis_products_2(basis,primitiveThr,overlapThr)

%This function calculates shell pair information for further use in the two
%electron integral loops.

%The function takes a very small amount of time compared with the electron
%repulsion and even the one electron integrals, so there is not much point
%in optimizing it.
nshells = size(basis,1);

%The order of the elements in pair_data(i,j) differs from the order in pair_data(j,i), but the elements are the same.
%Need to check later if I can just save the lower triangular form.


pair_data = cell(nshells,nshells); %I need a cell structure because I have many fields
RABdata = cell(nshells,nshells);
pOverlap = zeros(nshells,nshells);
%basis{a}(1) = number of primitives
%basis{a}(2) = angular momentum L
%basis{a}(3:5) = x0,y0,z0
%basis{a}(6:6+DimL-1) = relative normalization factors (s -> DimL = 1, p -> DimL = 3, d -> DimL = 6, etc)
%basis{a}(6+DimL+3*(na-1)+0) = alpha (na-1 = 0,1,2, number of primitives - 1)
%basis{a}(6+DimL+3*(na-1)+1) = coefficient c (na-1 = 0,1,2, number of primitives - 1)
%basis{a}(6+DimL+3*(na-1)2) = normalization N (na-1 = 0,1,2, number of primitives - 1)

for a = 1:nshells
    for b = 1:nshells
        Nintsab = basis{a}(1)*basis{b}(1);
        KabValues = zeros(Nintsab,1);
        pValues = zeros(Nintsab,1);
        WeightValuesab = zeros(Nintsab,1);

        PxValues = zeros(Nintsab,1);
        PyValues = zeros(Nintsab,1);
        PzValues = zeros(Nintsab,1);

        RPAValues = zeros(Nintsab,3); %This will be row vectors
        RPBValues = zeros(Nintsab,3);

        xa = basis{a}(3);
        ya = basis{a}(4);
        za = basis{a}(5);

        xb = basis{b}(3);
        yb = basis{b}(4);
        zb = basis{b}(5);

        RABdata{a,b} = [xa-xb;ya-yb;za-zb];

        t = 1;
        for na = 1:basis{a}(1)
                La = basis{a}(2);
                Dima = (La+1)*(La+2)/2;
                alphaa = basis{a}(6+Dima+(na-1)*3+0);
                ca = basis{a}(6+Dima+(na-1)*3+1);
                Na = basis{a}(6+Dima+(na-1)*3+2);
                for nb = 1:basis{b}(1)
                        Lb = basis{b}(2);
                        Dimb = (Lb+1)*(Lb+2)/2;
                        alphab = basis{b}(6+Dimb+(nb-1)*3+0);
                        cb = basis{b}(6+Dimb+(nb-1)*3+1);
                        Nb = basis{b}(6+Dimb+(nb-1)*3+2);

                        p = alphaa+alphab;

                        Px = (alphaa*xa + alphab*xb)/p;
                        Py = (alphaa*ya + alphab*yb)/p;
                        Pz = (alphaa*za + alphab*zb)/p;

                        pValues(t) = p;
                        PxValues(t) = Px;
                        PyValues(t) = Py;
                        PzValues(t) = Pz;

                        RPAValues(t,:) = [Px-xa,Py-ya,Pz-za];
                        RPBValues(t,:) = [Px-xb,Py-yb,Pz-zb];

                        rhoAB = alphaa*alphab/p;
                        Kab = exp(-rhoAB*((RABdata{a,b}(1))^2+(RABdata{a,b}(2))^2+(RABdata{a,b}(3))^2));

                        KabValues(t) = Kab;
                        WeightValuesab(t) = ca*Na*cb*Nb;

                        t = t + 1;
                end %nb
        end %na

        pPxValues = pValues.*PxValues;
        pPyValues = pValues.*PyValues;
        pPzValues = pValues.*PzValues;

        pair_data{a,b} = [pValues PxValues PyValues PzValues KabValues WeightValuesab RPAValues RPBValues pPxValues pPyValues pPzValues];
        % if isempty(pair_data{a,b})
        %         pair_data{a,b} = zeros(1,15);
        %         pOverlap(a,b) = 0;
        % else
        %         pOverlap(a,b) = sum(pi^1.5*(pValues.^(-1.5)).*KabValues.*WeightValuesab);
        % end

        % if abs(pOverlap(a,b)) < overlapThr
        %         %Here I'm dividing the matrix pair_data into halves, according to large and low pValues.
        %         pMiddle = mean(pValues);
        %         indHigh = pValues > pMiddle;
        %         indLow  = pValues <= pMiddle;
        %
        %         pOverlapLow = (pi^1.5*(pValues(indLow).^(-1.5)).*KabValues(indLow).*WeightValuesab(indLow));
        %         pValueLow = sum(pValues(indLow).*pOverlapLow)/sum(pOverlapLow);
        %         PxValueLow = sum(PxValues(indLow).*pOverlapLow)/sum(pOverlapLow);
        %         PyValueLow = sum(PyValues(indLow).*pOverlapLow)/sum(pOverlapLow);
        %         PzValueLow = sum(PzValues(indLow).*pOverlapLow)/sum(pOverlapLow);
        %         RPAValueLow = sum(RPAValues(indLow,:).*pOverlapLow,1)/sum(pOverlapLow); %This works due to broadcasting, remember RPAValues is a (Nintsab,3) matrix
        %         RPBValueLow = sum(RPBValues(indLow,:).*pOverlapLow,1)/sum(pOverlapLow);
        %         pPxValueLow = sum(pPxValues(indLow).*pOverlapLow)/sum(pOverlapLow);
        %         pPyValueLow = sum(pPyValues(indLow).*pOverlapLow)/sum(pOverlapLow);
        %         pPzValueLow = sum(pPzValues(indLow).*pOverlapLow)/sum(pOverlapLow);
        %         %I think I need to sum WeightValues, and also Kab
        %
        %         KabValueLow = sum(KabValues(indLow));
        %         WeightValueabLow = sum(WeightValuesab(indLow));
        %
        %         pOverlapHigh = (pi^1.5*(pValues(indHigh).^(-1.5)).*KabValues(indHigh).*WeightValuesab(indHigh));
        %         pValueHigh = sum(pValues(indHigh).*pOverlapHigh)/sum(pOverlapHigh);
        %         PxValueHigh = sum(PxValues(indHigh).*pOverlapHigh)/sum(pOverlapHigh);
        %         PyValueHigh = sum(PyValues(indHigh).*pOverlapHigh)/sum(pOverlapHigh);
        %         PzValueHigh = sum(PzValues(indHigh).*pOverlapHigh)/sum(pOverlapHigh);
        %         RPAValueHigh = sum(RPAValues(indHigh,:).*pOverlapHigh,1)/sum(pOverlapHigh);
        %         RPBValueHigh = sum(RPBValues(indHigh,:).*pOverlapHigh,1)/sum(pOverlapHigh);
        %         pPxValueHigh = sum(pPxValues(indHigh).*pOverlapHigh)/sum(pOverlapHigh);
        %         pPyValueHigh = sum(pPyValues(indHigh).*pOverlapHigh)/sum(pOverlapHigh);
        %         pPzValueHigh = sum(pPzValues(indHigh).*pOverlapHigh)/sum(pOverlapHigh);
        %         %I think I need to sum WeightValues, and also Kab
        %
        %         KabValueHigh = sum(KabValues(indHigh));
        %         WeightValueabHigh = sum(WeightValuesab(indHigh));
        %         try
        %         pair_data{a,b} = [pValueLow PxValueLow PyValueLow PzValueLow KabValueLow WeightValueabLow RPAValueLow RPBValueLow pPxValueLow pPyValueLow pPzValueLow;...
        %                           pValueHigh PxValueHigh PyValueHigh PzValueHigh KabValueHigh WeightValueabHigh RPAValueHigh RPBValueHigh pPxValueHigh pPyValueHigh pPzValueHigh];
        %         catch
        %                 disp(indHigh);
        %                 disp(' ');
        %                 disp(indLow);
        %                 disp(' ');
        %                 disp(RPAValues);
        %                 disp(' ');
        %                 disp(RPBValues);
        %                 disp(' ');
        %                 disp(RPAValueLow);
        %                 disp(' ');
        %                 disp(RPBValueLow);
        %                 disp(' ');
        %                 disp(RPAValueHigh);
        %                 disp(' ');
        %                 disp(RPBValueHigh);
        %         end_try_catch
        % endif
        %uncroppedSize = length(KabValues);
        % indKab_aboveThr = find(pair_data{a,b}(:,5) > primitiveThr);
        % pair_data{a,b} = pair_data{a,b}(indKab_aboveThr,:);
        % croppedSize = length(indKab_aboveThr);
        %
        % if isempty(pair_data{a,b})
        %         pair_data{a,b} = zeros(1,15);
        %         pOverlap(a,b) = 0;
        % else
        %         pOverlap(a,b) = sum(pi^1.5*(pValues.^(-1.5)).*KabValues.*WeightValuesab);
        % end
        %if uncroppedSize-croppedSize > 0
        %        disp(['pair_data(' num2str(a) ',' num2str(b) ') has been cropped by ' num2str(uncroppedSize-croppedSize) ' primitives']);
        %end
    end %b
end %a




end %function
