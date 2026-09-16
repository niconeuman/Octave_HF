function Asph = cart2sph(A)

%permA = permute(A,[]);

%        dxx dxy dxz dyy dyz dzz
%dx2-y2
%dz2
%dxy
%dxz
%dyz

cart2sphP = [0 0 1
            1 0 0
            0 1 0];

cart2sphD = [-1/2 0 0 -1/2 0 1
              0 0 sqrt(3)          0 0 0
              0 0 0 0 sqrt(3)          0
              sqrt(3)/2 0 0 -sqrt(3)/2 0 0
              0 sqrt(3) 0 0          0  0];



if size(A,1) == 6
    Atemp = cart2sphD*A;
elseif size(A,1) == 3
    Atemp = cart2sphP*A;
else
    Atemp = A;
end

if size(Atemp,2) == 6
    Asph = Atemp*cart2sphD';
elseif size(Atemp,2) == 3
    Asph = Atemp*cart2sphP';
else
    Asph = Atemp;
end

end
