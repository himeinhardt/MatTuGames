function pc=parity_coeff(clv)
% PARITY_COEFF computes the parity transform of the TU-game v.
% Source: U. Faigel and M. Grabisch (2014). 
%
% Usage: pc=parity_coeff(clv)
%

% Define variables:
%  output:
%  pc       -- The parity coefficients of the TU-game v.
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
%   02/25/2014        0.5             hme
%                
    
v=clv.tuvalues;
n=clv.tuplayers;
pb=parity_basis(n);
pc=v*pb;

