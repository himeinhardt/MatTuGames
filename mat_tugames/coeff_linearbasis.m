function [lcf, lb]=coeff_linearbasis(v,method)
% COEFF_LINEARBASIS computes the coeffients of the linear basis of
% game v.
% For n>12 this function needs some time to complete.
%
% Usage: [lcf lb]=coeff_linearbasis(v,method)
%
% Define variables:
%  output:
%  lcf      -- Coeffients of the linear basis of a game v.
%  lb       -- linear basis of game v.
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%  method   -- A string to format the matrix. Permissible methods
%              to format the matrices are 'full','sparse' or the
%              empty string '', to invoke the default, which is sparse.

%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/21/2013        0.4             hme
%   02/15/2014        0.5             hme
%                

N=length(v);
[~, n]=log2(N);
if nargin < 2
   method='sparse';
else
   if isempty(method)
      method='sparse';
   end
end

lb=linear_basis(n,method);
lcf=(lb\v')';
if nargout==2
   if strcmp(method,'full')
      lb=full(lb);
   end
end
