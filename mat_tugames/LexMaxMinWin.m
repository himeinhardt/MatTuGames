function [sS,imat,lxS]=LexMaxMinWin(clm,n)
% LexMaxMinWin determines the lexicographical maximal coalition from the clm represented as a cell, matrix, or array of integers.
% The input clm is first sorted in accordance with the shortest size, and within the block of shortest size the lexicographical maximum is selected.    
%
% Usage: sS=LexMaxMinWin(clm,n);
%
% Define variables:
%  output:
%  sS       -- The lex-max coalition of clm.
%  imat     -- The corresponding incidence vector of sS.
%  lxS      -- The dual number associated with that coalition.
%
%  input:
%  clm      -- A cell, matrix, or array (vector) which contains the information about a set of coalitions,
%              for instance, the set of minimal winning coalitions.    
%  n        -- Number of player involved, must be an integer.
%
%    
% Examples:
% E1.)
% Coalitions of unequal size must be represented as a cell input:    
% clm={[1 2 5 6] [1 3 4], [1 6 7], [2 3 6 7], [2 4 7], [2 5 7], [4 5 6 7]};
%
%   [sS,imat,lxS]=LexMaxMinWin(clm,7)
% which returns sS =
%
%    13
%
% the coalition of lexicographical maximum.
%
% The type of the set of coalitions is grasped by its incidence vector given by
%
% imat =
%
%   1   0   3   4   0   0   0
% 
% and the dual number is 
% 
% lxS = 
%
%   102    
%    
% E2.)
% Coalitions of equal size can be represented as a matrix
% clmat=[1 2 5; 1 3 4; 1 6 7; 2 3 6; 2 4 7; 2 5 7; 4 5 6;]
%  
%    [sS2,imat2,lxS2]=LexMaxMinWin(clmat,7);
% which returns sS2 =
%
%    19
%
% which is the lexicographical maximal coalition of clmat. 
% And as a second output the incidence vector imat2 is returned
%
% imat2 =
%
%   1   2   0   0   5   0   0
%
% and the dual number is 
%
%  lxS2 = 
%
%        100    
%
%    
% E3.)
% Consider the following array of coalitions represented by its unique integer representation:
% sS5 =
%
%    1    2    3   13   21   25   14   22    
%
% then invoke 
%
%   [lS5,imat5,lxs5]=LexMaxMinWin(sS5,5)    
%
% to get again the lexicographical maximal coalition
%    
% lS5 =
%
%    1
%
% And the incidence vector is given by
%    
% imat5 =
%
%   1   0   0   0   0
%
% and finally the dual number is listed by
%
% lxs5 =
%
%   16
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
%   06/07/2022        1.9.1           hme
%                    
    
pl=1:n;
if iscell(clm)      
   sS=clToMatlab2(clm);
   csz=numel(clm);
   clm2=cell(1,csz);
   for k=1:csz
       cl(k)=sum(2.^(n-clm{k}));
   end
   [lxS,sidx]=max(cl);
elseif ismatrix(clm)
   [s1,s2]=size(clm);
   nz=nnz(clm);
   if s1*s2>nz
      for kk=1:s1
          cl{kk}=pl(clm(kk,:)>0);
      end
      [sS,imat]=LexMaxMinWin(cl,n);
      return;
   elseif s1>1 && s1*s2==nz 
      sS=clToMatlab2(clm);
      clm=n-clm;
      clm=2.^clm;
      [lxS,sidx]=max(sum(clm,2)');
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
   bd=length(sS);
   ov=ones(n,1);
   clsize=indM*ov;
   mcl=min(clsize);
   eqm=find(clsize==mcl);
   lc=length(eqm);
   if lc~=bd
      sS=sS(eqm);
      indM=indM(eqm,:);
   end
   expl=pl.*ones(lc,n);
   imat=indM.*expl;
   cl=zeros(1,bd);
   for kk=1:lc
       clm=pl(imat(kk,:)>0);
       clm=n-clm;
       clm=2.^clm;
       cl(kk)=sum(clm); %% dual numbers of coalitions sS.
   end
   [lxS,sidx]=max(cl);
   %% lexicographic order R lex T iff d(R) > d(T).
   %% Determining lexicographic maximal of coalitions sS.  
   sS=sS(sidx);
   imat=imat(sidx,:);
else    
  sS=sS(sidx);
  bd=length(sS);
  pl=1:n;    
  it=0:-1:1-n;
  indM=rem(floor(sS(:)*pow2(it)),2);
  ov=ones(n,1);
  clsize=indM*ov;
  mcl=min(clsize);
  eqm=find(clsize==mcl);
  lc=length(eqm);
  if lc~=bd
     sS=sS(eqm);
    indM=indM(eqm,:);
  end
  expl=pl.*ones(1,n);
  imat=indM.*expl;
  imat=imat(1,:);
  sS=sS(1);
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
    
    
    
