function sd_vl=SolidarityValue(v)
% SOLIDARITYVALUE computes the Solidarity value of a TU-game v.
%
% Usage: sd_vl=SolidarityValue(v)
% Define variables:
%  output:
%  sd_vl    -- The Solidarity value of a TU-game v.
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
%   07/19/2013        0.4             hme
%                
N=length(v);
[~, n]=log2(N);
wv=ones(1,n);
sd_vl=weightedSolidarity(v,wv);    
    
