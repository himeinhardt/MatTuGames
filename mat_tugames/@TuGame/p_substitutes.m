function sbs_v=p_substitutes(clv,tol)
% P_SUBSTITUTES establishes which pair of players are substitutes using
% Matlab's PCT.
%
% Usage: sbs_v=clv.p_substitutes(tol)
%
% Define variables:
% output:
% sbs        -- A matrix of maximal size(binom(n,2),2). Shows
%               in each row the pair that are substitutes.
%
%  input:
%  clv      -- TuGame class object.
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/05/2015        0.6             hme
%

if nargin <2
   tol=10^6*eps;
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
S=1:N;
sbs_v=cell(n,1);

parfor ii=1:n
    myTemp = zeros(n,2);
    for jj=1:n
       if ii < jj
         Sni=bitget(S,ii)==0;
         Snj=bitget(S,jj)==0;
         Snij=S(Sni & Snj);
         ci=2^(ii-1);
         cj=2^(jj-1);
         Si=bitor(Snij,ci);
         Sj=bitor(Snij,cj);
         eqQ=all(abs(v(Si)-v(Sj))<tol);
         if eqQ==1
            myTemp(jj,:)=[ii,jj];    
         end
       end  
    end
    sbs_v{ii} = myTemp; 
end   
%
% Reformatting result
% 
sbs_v = unique(cell2mat(sbs_v),'rows');
sbs_v(1,:)=[];
