function vmat=value_matrix(asm)
% VALUE_MATRIX computes from an assignment matrix the corresponding value matrix for a permutation game. 
%
% Usage: vmat=value_matrix(asm)
%
% Define variables:
%  output:
%  vmat     -- A value matrix for a permutation game.
%            
%  input:
%
%  asm      -- A square assignment matrix.
%              To construct one invokes for instance,
%              asm=magic(5);
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/25/2018        0.9             hme
%


[a1,a2]=size(asm);
m=a1+a2;
vmat=[zeros(a1,a1),asm;zeros(a2,m)];

