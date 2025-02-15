function [sS,imat,lxS]=CanonicalOrder(clm,n);
% CANONICALORDERORDER orders a set of coalitions represented by a cell, matrix, or array of integers
% with shortest size first, second shortest size second and so on. Entities having
% the same size are ordered lexicographically.    
%
% Usage: sS=CanonicalOrder(clm,n);
%
% Define variables:
%  output:
%  sS       -- The list of coalitions in unique integer representation sorted 
%              by a canonical order.
%  imat     -- The corresponding incidence matrix of sS.
%  lxS      -- The canonical numbers associated with coalitions.
%
%  input:
%  clm      -- The cell or matrix, which contains the coalition information
%              in generic power set representation (e.g. Mathematica format).   
%  n        -- Number of player involved, must be an integer.
%
%    
% Examples:
% E1.)    
% Coalitions of unequal size must be represented as a cell input:    
% clm={[1 2 5 6] [1 3 4], [1 6 7], [2 3 6 7], [2 4 7], [2 5 7], [4 5 6 7]};
%
%   [sS,imat]=CanonicalOrder(clm,7)
% which returns sS =
%
%    13    97    74    82    51   102   120
%
% the set of coalitions into a canonical order given by its unique integer representation.    
%
% The canonical order of the set of coalitions is grasped by its incidence matrix given by
%
% imat =
%
%   1   0   3   4   0   0   0
%   1   0   0   0   0   6   7
%   0   2   0   4   0   0   7
%   0   2   0   0   5   0   7
%   1   2   0   0   5   6   0
%   0   2   3   0   0   6   7
%   0   0   0   4   5   6   7
%    
% E2.)   
% Coalitions of equal size can be represented as a matrix
% clmat=[1 2 5; 1 3 4; 1 6 7; 2 3 6; 2 4 7; 2 5 7; 4 5 6;]
%  
%    [sS2,imat2]=CanonicalOrder(clmat,7);
% which returns sS2 =
%
%    19   13   97   38   74   82   56
%
% Again we get the set of coalitions into canonical order given by its unique integer representation.
% And as a second output the incidence matrix imat2 in canonical order, which returns
%
% imat2 =
%
%   1   2   0   0   5   0   0
%   1   0   3   4   0   0   0
%   1   0   0   0   0   6   7
%   0   2   3   0   0   6   0
%   0   2   0   4   0   0   7
%   0   2   0   0   5   0   7
%   0   0   0   4   5   6   0
% 
% E3.)
% Consider the following array of coalitions represented by its unique integer representation:
% sS5 =
%
%    1    2    3   13   21   25   14   22    
%
% then invoke 
%
%   [lS5,imat5,lxs5]=LexOrder2(sS5,5)    
%
% to get the canonical sorted data array of coalitions given by
%    
% lS5 =
%
%    3   13   21   25    1   14   22    2
%
% The incidence matrix is given by
%    
% imat5 =
%
%   1   2   0   0   0
%   1   0   3   4   0
%   1   0   3   0   5
%   1   0   0   4   5
%   1   0   0   0   0
%   0   2   3   4   0
%   0   2   3   0   5
%   0   2   0   0   0
%
% and the array of sorted canonical numbers (descend) is listed below.
%
% lxs5 =
%
%   24   22   21   19   16   14   13    8    
%
%    
%    
    

%    
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/02/2022        1.9.1           hme
%                    
pl=1:n;   
if iscell(clm)      
   sS=clToMatlab2(clm);
   csz=numel(clm);
   clm2=cell(1,csz);
   for k=1:csz
       vc=zeros(1,n);
       vc(clm{k})=clm{k};
       cl(k)=sum(2.^(n-vc));
   end
   [lxS,sidx]=sort(cl,'descend');
elseif ismatrix(clm)
   [s1,s2]=size(clm);
   nz=nnz(clm);
   if s1*s2>nz
      for kk=1:s1
          cl{kk}=pl(clm(kk,:)>0);
      end
      [sS,imat,lxS]=CanonicalOrder(cl,n);
      return;      
   elseif s1>1 && s1*s2==nz
      sS=clToMatlab2(clm);
      cl=zeros(1,s1);
      for kk=1:s1
          vc=zeros(1,n); 
          vc(clm(kk,:))=clm(kk,:);
          rv=n-vc;
          rv=2.^rv;
          cl(kk)=sum(rv);
      end    
      [lxS,sidx]=sort(cl,'descend');
   else
      sS=clm; 
      lxS=[];
      clm=[];
   end    
else
   error('The input format must either be a cell or matrix!') 
end    

if isempty(lxS)
   it=0:-1:1-n;
   indM=rem(floor(sS(:)*pow2(it)),2);
   lS=length(sS);
   expl=pl.*ones(lS,n);
   imat=indM.*expl;
   clm=n-imat;
   clm=2.^clm;
   cl=sum(clm,2)'; %% canonical numbers of coalitions S.
   %% canonical order R > T iff d(R) > d(T).
   [st,sidx]=sort(cl,'descend');
   %% Determining canonical order of coalitions S.  
   sS=sS(sidx);
   imat=imat(sidx,:);
   lxS=cl;
else    
  sS=sS(sidx);
  lS=length(sS);
  pl=1:n;    
  it=0:-1:1-n;
  indM=rem(floor(sS(:)*pow2(it)),2);
  expl=pl.*ones(lS,n);
  imat=indM.*expl;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sS=clToMatlab2(clm)
% CLTOMATLAB2 same as clToMatlab but unsorted.  
%
%    
if iscell(clm)
   csz=numel(clm);
   clm2=cell(1,csz);
   for k=1:csz
       sS(k)=sum(2.^(clm{k}-1));
   end
 elseif ismatrix(clm) 
  clm=clm-1;
  clm=2.^clm;
  sS=sum(clm,2)';
else
  error('The input format must either be a cell or matrix!')
end
