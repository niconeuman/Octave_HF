function [xValues,KabValues,RPAValues,RPBValues,RPCValues,pValues,WeightValues,ZValues] = loopPrimitivesNuc_Octave(basis_a,basis_b,Z,AL)

nAtom = size(Z,2);

Nints = basis_a(1)*basis_b(1)*nAtom;

xValues = zeros(Nints,1);

KabValues = zeros(Nints,1);
RPAValues = zeros(Nints,3);
RPBValues = zeros(Nints,3);
RPCValues = zeros(Nints,3);
pValues = zeros(Nints,1);

WeightValues = zeros(Nints,1);
ZValues = zeros(Nints,1);

  t = 1;

  xa = basis_a(3);
  ya = basis_a(4);
  za = basis_a(5);

  xb = basis_b(3);
  yb = basis_b(4);
  zb = basis_b(5);

  RAB = [xa-xb;ya-yb;za-zb]; %column vector

  for na=1:basis_a(1) %loops over the number of primitives in the 1st contracted basis function
      La = basis_a(2);
      Dima = (La+1)*(La+2)/2;
      alphaa = basis_a(6+Dima+(na-1)*3+0);
      ca = basis_a(6+Dima+(na-1)*3+1);
      Na = basis_a(6+Dima+(na-1)*3+2);



      for nb=1:basis_b(1) %loops over number of primitives in the 2nd contracted basis function
          Lb = basis_b(2);
          Dimb = (Lb+1)*(Lb+2)/2;
          alphab = basis_b(6+Dimb+(nb-1)*3+0);
          cb = basis_b(6+Dimb+(nb-1)*3+1);
          Nb = basis_b(6+Dimb+(nb-1)*3+2);
          % disp(alphab);
          % disp(cb);
          % disp(Nb);



          p = alphaa+alphab;

          Px = (alphaa*xa + alphab*xb)/p;
          Py = (alphaa*ya + alphab*yb)/p;
          Pz = (alphaa*za + alphab*zb)/p;

          RPA = [Px-xa,Py-ya,Pz-za];
          RPB = [Px-xb,Py-yb,Pz-zb];



          rhoAB = alphaa*alphab/p;
          Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));

          for N = 1:nAtom

              RPC = [Px-AL(N,1),Py-AL(N,2),Pz-AL(N,3)]; %column vector
              RPC2 = RPC(1)^2+RPC(2)^2+RPC(3)^2;

              x = p*RPC2;

              ZValues(t) = Z(N);
              xValues(t) = x;
              KabValues(t) = Kab;
              RPAValues(t,:) = [RPA(1),RPA(2),RPA(3)];
              RPBValues(t,:) = [RPB(1),RPB(2),RPB(3)];
              RPCValues(t,:) = [RPC(1),RPC(2),RPC(3)];
              pValues(t) = p;
              WeightValues(t) = ca*Na*cb*Nb;

              t = t + 1;


          end
      end
  end




end
