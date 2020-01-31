function sexd=SED(v,x)
% SED computes the small excess difference w.r.t. the payoff x.
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
% Usage: sexd=SED(v,x)
%
% Define variables:
%
%  output:
%  sexd     -- Returns small excess difference,
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
   sexd=-inf; %% -inf - (inf)
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
             vs1(ii)=max(ex_v(k-1),ex_v(jj-1))-ex_v(k-1)-ex_dv(jj-1);
           elseif k==1 && jj >1
             ii=N1*(jj-1);
             vs1(ii)=max(ex_v(1),ex_v(jj-1))-ex_v(1)-ex_dv(jj-1) ;
           elseif k>1 && jj==1
             ii=k-1;
             vs1(ii)=max(ex_v(k-1),ex_v(1))-ex_v(k-1)-ex_dv(1);
           end
    end
end
sS=find(vs1==-inf);
vs1(sS)=[];
sexd=max(vs1);

