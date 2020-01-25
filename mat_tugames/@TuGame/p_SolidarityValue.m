function sd_vl=p_SolidarityValue(clv)
% P_SOLIDARITYVALUE computes the Solidarity value of a TU-game v
% using Matlab's PCT.
%
% Usage: sd_vl=clv.p_SolidarityValue()
% Define variables:
%  output:
%  sd_vl    -- The Solidarity value of a TU-game v.
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
%   07/26/2013        0.4             hme
%                
n=clv.tuplayers;
wv=ones(1,n);
sd_vl=clv.p_weightedSolidarity(wv);    
    
