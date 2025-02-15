function PRK=GetPrk(v,x,tol,str)
% GETPRK tries to find an additional pre-kernel elements for game v
% from the DM-reduced game v_x of a set of partitions of the player set. 
% 
%  Usage: PRKS=GetPrk(v,x)
%
%  Inspired by J. Arin and I. Katsev (2013),The coincidence of the kernel and nucleolus of a convex game: An alternative proof
%
%
% Define variables:
%  output: Fields
%  prks     -- Returns a pre-kernel matrix of size (tmax,n) for
%              game v.
%  pkQ      -- Returns a list of true/false values of size (1,tmax) to indicate whether the row belongs to 
%              the pre-kernel (1) or not (0).
%  x        -- Returns the original pre-kernel element.
%  dmrd     -- Returns the list of Davis-Maschler reduced games w.r.t. the partition set and x.
%  ptn      -- Retunrs the set of partitions which are included in D(x), the set of maximal excesses.
%  tol      -- Tolerance value. By default, it is set to 10^7*eps.
%             (optional)
%  str      -- String value to replicate the original solution. Admissible values are: 
%               'repl' to replicate the orignal pre-kernel element,
%               'mult' to try to find another pre-kernel,
%              default is "repl".
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- Payoff vector of size(1,n), pre-kernel vector.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/21/2022        1.9.1           hme
%


if nargin<2
   x=PreKernel(v);
   tol=10^6*eps;
   str='repl';
elseif nargin<3
   tol=10^6*eps;
   str='repl';
elseif nargin<4
   if isempty(tol)
      tol=10^6*eps;
   end	   
   str='repl';  	
end	

N=length(v);
[~, n]=log2(N);
pl=1:n;

SATIS=satisfaction(v,x);
ptn=SATIS.ptn';
if isempty(ptn)
   disp('No Partition detected; no additional pre-kernel element can be computed!');	
   return;
end	
sz=size(ptn);
v_t=cell(sz(1),sz(2));
pkQ=false(1,sz(1));
z=zeros(sz(1),n);

for ll=1:sz(1)
  for k=1:sz(2)
      sS=ptn(ll,k);	
      inS=bitget(sS,pl)==1;
      if strcmp('repl',str); % must replicate the original pre-kernel x.      
         v_t{ll,k}=RedGame(v,x,sS);
      elseif strcmp('mult',str)  % trying to get another pre-kernel element.
         v_t{ll,k}=CvRedGame(v,x,sS);
      else  % default is replicating.
         v_t{ll,k}=RedGame(v,x,sS);
      end	      
      pk_vx=PreKernel(v_t{ll,k});
      z(ll,inS)=pk_vx;
  end 
  pkQ(ll)=PrekernelQ(v,z(ll,:));  
end

PRK.prks=z;
PRK.pkQ=pkQ;
PRK.x=x;
PRK.dmrd=v_t;
PRK.ptn=ptn;
