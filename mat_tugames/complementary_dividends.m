function cd=complementary_dividends(v);
% COMPLEMENTARY_DIVIDENDS computes the complementary dividends.
% For n>14 this function needs some time to complete.
%
% Usage: cd=complementary_dividends(v)
%
% Define variables:
%  output:
%  cd       -- Complementary dividends.
%
%  input:
%  v        -- A TU-game of length 2^n-1.
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


N=length(v);
[~, n]=log2(N);


cb=complementary_basis(n,'sparse');
cd=(inv(cb)*v')';
