function iacv=IrredAntiCore(cm,method)
% IRREDANTICORE computes from a cost matrix the corresponding extreme points of 
% the irreducible anti-core of the associated m.c.s.t. game. 
% Using Prim's or Kruskal’s algorithm.
%
% Usage: iacv=IrredAntiCore(cm,method)
%
% Define variables:
%  output:
%  iacv     -- Extreme points of the irreducible anti-core of the associated m.c.s.t. game
%            
%  input:
%
%  cm       -- A square cost matrix (n+1xn+1) derived from a cost
%              spanning graph problem. For instance, for a four
%              person game the size of the matrix must be (5x5). The source
%              is player 1.
%  method   -- A string to define the method. Admissible methods are:
%              'prim' 
%              'kruskal'
%              The default method is 'prim'.
%

%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/22/2021        1.9             hme
%

if nargin < 2
   method='prim';
end


vc=mcst_game(cm,method);
imvc=IrredCostMatrix(cm,vc.trm{end});
irdv=mcst_game(imvc.icm,method);
iacv=CddAntiCoreVertices(irdv.c_v);
