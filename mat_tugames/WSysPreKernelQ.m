function [cQ,wpk_v,wpkQ,pkQ_wex]=WSysPreKernelQ(v,pS,x,tol)
% WSYSPRENUCLQ checks the correctness of the pre-kernel element of a weight system pS.
% Requires Matlab's Optimization toolbox (default), otherwise CPLEX.
%
% Source: Remark 2.4 (Kleppe, J.; Reijnierse, J.H.; Sudhölter, P., 2013)
%
% Usage: [cQ,wpn_v,nwex,wexQ]=WSysPreKernelQ(v,pS,tol) 
%
% Define variables:
%  output:
%  cQ       -- Returns true (1) whenever the weighted pre-kernel element w.r.t.
%              the Tu-game v and a weight system pS was correctly evaluated, 
%              otherwise false (0).
%  wpk_v    -- Returns the weighted pre-kernel point w.r.t. v and pS
%  wpkQ     -- Returns true (1) if the pre-imputation is a pre-kernel element 
%              of the wight system pS, otherwise false (0). 
%  pkQ_wex  -- Returns true (1) if the pre-kernel imputation of the excess game
%              is a null vector, otherwise false (0).
%  
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  pS       -- A vector of weights of length 2^n-1.
%  x        -- An imputation of length n.
%  tol      -- Specifies a tolerance value, default is 10^7*eps.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/14/2014        0.5               hme
%   01/27/2024        1.9.1             hme
%
%
narginchk(2,3); % check for legal number of input arguments.

if nargin<3
  N=length(v);
  [~, n]=log2(N);
  if (2^n-1)~=N
    error('Game has not the correct size!');
  end
  np=length(pS);
  if np~=N
    error('Weight vector has not the correct size!');
  end
  tol=10^7*eps;
  x='';
elseif nargin<4
  N=length(v);
  [~, n]=log2(N);
  if (2^n-1)~=N
    error('Game has not the correct size!');
  end
  np=length(pS);
  if np~=N
    error('Weight vector has not the correct size!');
  end
  tol=10^7*eps;
else
  N=length(v);
  [~, n]=log2(N);
  if (2^n-1)~=N
    error('Game has not the correct size!');
  end
  np=length(pS);
  if np~=N
    error('Weight vector has not the correct size!');
  end
end

if isempty(x)
  wpk_v=weightedPreKernel(v,pS);
else
  wpk_v=x;
end

wpkQ=weightedPrekernelQ(v,wpk_v,pS);
wex=weightedExcess(v,wpk_v,pS);
zv=zeros(1,n);
pkQ_wex=PrekernelQ(wex,zv);
cQ=pkQ_wex & wpkQ;
