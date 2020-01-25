function SOL=tug_DetUCoord(hd,n)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_DetUCoord(avc,n)
% Define variables:
%  output:
%  SOL        -- Returns the unanimity coordinates.
%                Field variable gives result in Matlab and Mathematica format.
%  core_vert  -- Matrix of core vertices. Output is numeric or a string.
%  crst       -- The core constraints.
%  input:
%  hd         -- A list of positive unanimity coordinates of length binom(n,2) or
%                a vector of unanmity coordinates of length 2^n-1. The coordinates
%                of coalition size 2 will be extracted.
%  n          -- The number of persons involved.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/06/2011        0.1 beta        hme
%

% Here we assume that the user has represented the game correctly.

if nargin<1
    error('At least the unanimity coordinates must be given!');
elseif nargin <2
    N=length(hd);
    gr=dec2bin(N);
    n=length(gr);
    if (2^n-1)~=N
      error('Unanimity coordninates have not the correct size!');
    end
    whd=gameToMama(hd);
    bn=binomial(n,2);
    sc2=n+bn;
    fc2=n+1;
    wuc=whd(:,fc2:sc2);
    lc=wuc<0;
    if any(lc)
    warning('Some unanimity coordinates of coalition size 2 have negative values!');
    warning('No average-convex game can be constructed from the set of unanmity coordinates!')
    else
    end
    SOL=avcoord(wuc,n);
else
    N=length(hd);
    gr=dec2bin(N);
    bn=binomial(n,2);
    switch N
      case bn
         lc=hd<0;
         if any(lc)
            warning('Some unanimity coordinates of coalition size 2 have negative values!');
            warning('No average-convex game can be constructed from the set of unanmity coordinates!')
         else
         end
         SOL=avcoord(hd,n);
      case (2^n-1)
         whd=gameToMama(hd);
         bn=binomial(n,2);
         sc2=n+bn;
         fc2=n+1;
         wuc=whd(:,fc2:sc2);
         lc=wuc<0;
         if any(lc)
            warning('Some unanimity coordinates of coalition size 2 have negative values!');
            warning('No average-convex game can be constructed from the set of unanmity coordinates!')
         else
         end
         SOL=avcoord(wuc,n);
      otherwise
       error('Unanimity coordninates have not the correct size!');
    end
end

function SOL=avcoord(wuc,n)

math('quit')
pause(1)
math('$Version')
math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }');
disp('Passing unanimity coordninates to Mathematica ...')
math('matlab2math','n1',n);
math('matlab2math','mhd',wuc);
math('coord=Flatten[mhd,1]');
math('bds=Flatten[n1][[1]]');
math('T=Flatten[Range[n1]]');
msvg=math('mbv=DetUCoord[coord,T]');
svg=math('math2matlab','mbv');
svg(:,1)=[];
sv_g=gameToMatlab(svg);
SOL=struct('UnanimityCoordninates',sv_g,'MUnanimityCoordninates',msvg);
math('quit')
