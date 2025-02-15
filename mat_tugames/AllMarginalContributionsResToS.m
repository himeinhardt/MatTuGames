function [Mgc,shv]=AllMarginalContributionsResToS(v,S)
% ALLMARGINALCONTRIBUTIONSRESTOS computes all marginal worth vectors 
% of a TU-game v restricted to coalition S.
%
% Usage: [Mgc shv]=AllMarginalContributionsResToS(v,S)
% Define variables:
%  output:
%  Mgc      -- The matrix of marginal contributions restricted to S.
%  shv      -- Shapley value of the sub-game vS.
%    
%  input:
%  v        -- A TU-Game of length 2^n-1.
%  S        -- Coalition S given by its unique integer representation.
%              Thus, an integer must be given.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/11/2021        1.9.1           hme
%

N=length(v);
[~, n]=log2(N);
pl=1:n;
spm=perms(pl);
[s1,s2]=size(spm);
a=pl(logical(bitget(S,pl)));
la=length(a);
Mgc=zeros(s1,la);

for k=1:s1
    Si=[];
    rw=spm(k,:);
    for ii=1:la
	Soi=[];    
	for jj=1:la    
           if rw(a(jj))<rw(a(ii))
              Soi=[Soi,a(jj)];
           end 	   
        end
     Si=[Soi,a(ii)];
     soi=sum(2.^(Soi-1));
     si=sum(2.^(Si-1));
     if soi>0
       Mgc(k,ii)=v(si)-v(soi);
      else
       Mgc(k,ii)=v(si);	    
     end   
    end
end
shv=sum(Mgc)/s1;
