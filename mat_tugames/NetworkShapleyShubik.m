function nSSI=NetworkShapleyShubik(E,th,w)
% NETWORKSHAPLEYSHUBIK computes the network Shapley-Shubik power index from the set of winning coalitions of a network E
% while imposing a threshold of th.
%
% Inspired by the Source: M. Holler and F. Rupp, Power in Networks: The Medici,2020.
%
% Usage: nJI=NetworkShapleyShubik(E,th)
% Define variables:
%  output:
%  nSSI     -- Network Shapley-Shubik Powr Index.
%
%  input:
%  E        -- An edge matrix of size (lx2) or a cell of numel l.
%              The source must be given by 1, and the sink by the
%              number of the player set. However, the edge matrix
%              can also be of size (lx3) then c can be empty.
%  w        -- A vector of power weights.
%  th       -- Threshold to pass a bill (positive number) greater than 2 and 
%              not greater than or equal to n.
%
% Example:
% Define a matrix of edges given by
% E =
%   1   1   1   1   2   2   3   4
%   2   3   4   5   3   4   4   5
%
% or equivalently
%
% E={[1 2];[1 3];[1 4];[1 5];[2 3];[2,4];[3 4];[4 5]}
% Do not forget to set the separator by semicolon (;) not by comma.
%
% We have 8 edges here. Furthermore, in total we have 5 vertices, and one
% source given by number 1, and a sink by number 5. 
% Set th=3 for the threshold to pass a bill.
% Optinally, one can specify the voting weight of each player by a weighted vector, for instance,
% by w=[2 5 6 2 4]. Otherwise, the default of equal weight is assumed.
%
% Here, we have a player set of {1,2,3,4,5}, hence n=5.
%
% Now, invoke
%  nSSI=NetworkShapleyShubik(E,th)
% to get
%
%  0.2500    0.1667    0.1667    0.2500    0.1667
%
% For the example with weights vector w and th=10, we get
% nSSI=NetworkShapleyShubik(E,10,w)
%
%  0.1000    0.2667    0.3500    0.1000    0.1833
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/30/2020        1.9             hme
%

if nargin <2
   error('The threshold must be given');
elseif nargin < 3
   v=NetworkMajorityGame(E,th);
else
   v=NetworkMajorityGame(E,th,w);
end

nSSI=ShapleyValue(v);
