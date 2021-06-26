function ecg=p_ECGValue(clv)
%P_ECGVALUE computes the Equal Collective Gains value of a TU-game v using Matlab's PCT.
%
% Source: Emilio Calvo Ramón and Esther Gutiérrez-López (2020), 
%         The Equal Collective Gains value in Cooperative Games
%
% Usage: ecg_vl=clv.p_ECGValue()
% Define variables:
%  output:
%  ecg      -- The Equal Collective Gains value of a TU-game v.
%
%  input:
%  clv      -- TuGame class object.
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

N=clv.tusize;
n=clv.tuplayers;
v=clv.tuvalues;
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
