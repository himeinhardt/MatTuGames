function lexd=LED(v,x)
% LED computes the large excess difference w.r.t. the payoff x.
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
%          Sudhoelter (1997), The modified nucleolus: Properties and axiomatizations. International Journal of Game Theory, 26
%          (2):147–182, Jun 1997. ISSN 1432-1270. doi: 10.1007/BF01295846. URL https://doi.org/10.1007/BF01295846.
%
% Usage: lexd=LED(v,x)
%
% Define variables:
%
%  output:
%  lexd     -- Returns large excess difference,
%              
%  input:
%  v        -- A Tu-Game v of length 2^n-1.
%  x        -- payoff vector of size(1,n). Must be efficient.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%              (optional)
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/01/2018        1.0             hme
%

N=length(v);
[~, n]=log2(N);
dv=dual_game(v);
ex_v=excess(v,x);
ex_dv=excess(dv,x);

if N==1
   lexd=inf; %% inf - (-inf)
  return
end


N1=N+1;
n1=2*n;
N2=2^n1-1;
ii=1;
vs1=-inf(N2,1);
for k=1:N1-1
    for jj=1:N1-1
          if k>1 && jj >1
             ii=(k-1)+(jj-1)*N1;
             vs1(ii)=min(ex_v(k-1),ex_v(jj-1))-ex_v(k-1)-ex_dv(jj-1);
           elseif k==1 && jj >1
             ii=N1*(jj-1);
             vs1(ii)=min(ex_v(1),ex_v(jj-1))-ex_v(1)-ex_dv(jj-1) ;
           elseif k>1 && jj==1
             ii=k-1;
             vs1(ii)=min(ex_v(k-1),ex_v(1))-ex_v(k-1)-ex_dv(1);
           end
    end
end
sS=find(vs1==-inf);
vs1(sS)=[];
lexd=min(vs1);

