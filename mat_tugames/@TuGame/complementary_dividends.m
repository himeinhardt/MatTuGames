function cd=complementary_dividends(clv);
% COMPLEMENTARY_DIVIDENDS computes the complementary dividends.
% For n>14 this function needs some time to complete.
%
% Usage: cd=clv.complementary_dividends()
%
% Define variables:
%  output:
%  cd       -- Complementary dividends.
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
%   02/16/2022        1.9.1           hme
%                

v.tuvalues;
n=clv.tuplayers;

cb=complementary_basis(n,'sparse');
cd=(inv(cb)*v')';
