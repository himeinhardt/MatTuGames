function sd_vl=SolidarityValue(clv)
% SOLIDARITYVALUE computes the Solidarity value of a TU-game v.
%
% Usage: sd_vl=SolidarityValue(clv)
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
sd_vl=weightedSolidarity(clv,wv);    
    
