function ecg=p_ECGValue(v)
% P_ECGVALUE computes the Equal Collective Gains value of a TU-game v using Matlab's PCT.
%
% Source: Emilio Calvo Ramón and Esther Gutiérrez-López (2020), 
%         The Equal Collective Gains value in Cooperative Games
%
% Usage: ecg_vl=p_ECGValue(v)
% Define variables:
%  output:
%  ecg      -- The Equal Collective Gains value of a TU-game v.
%
%  input:
%  v        -- A TU-Game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/23/2020        1.9             hme
%


N=length(v);
[~, n]=log2(N);
ecg=zeros(1,n);
S=1:N;
pl=1:n;
dvS=zeros(1,N);
fc=zeros(1,N);
parfor ss=1:N
    kS=pl(bitget(ss,pl)==1);
    lkS=length(kS);
    fc(ss)=(factorial(n-lkS)*factorial(lkS-1))/(factorial(n-1)*lkS);
    tS=2.^(kS-1);
    if ss==tS %empty set Swk
       dvS(ss)=v(ss);
    else   
       dvS(ss)=v(ss) - sum(v(ss-tS))/(lkS-1);
    end   
end

parfor kk=1:n
    cl=S(bitget(S,kk)==1);
    afc=fc(cl);
    delv=dvS(cl);
    ecg(kk)=afc*delv';
end	
