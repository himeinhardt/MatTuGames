function [cQ,wpk_v,wpkQ,pkQ_wex]=WSysKernelQ(clv,pS,x,tol)
% WSYSPRENUCLQ checks the correctness of the kernel of a weight system pS.
% Requires Matlab's Optimization toolbox (default), otherwise CPLEX.
%
% Source: Remark 2.4 (Kleppe, J.; Reijnierse, J.H.; Sudhölter, P., 2013)
%
% Usage: [cQ,wpn_v,nwex,wexQ]=WSysKernelQ(clv,pS,tol) 
%
% Define variables:
%  output:
%  cQ       -- Returns true (1) whenever the weighted kernel element w.r.t.
%              the Tu-game v and a weight system pS was correctly evaluated, 
%              otherwise false (0).
%  wpk_v    -- Returns the weighted kernel element w.r.t. v and pS
%  wpkQ     -- Returns true (1) if the pre-imputation is a kernel point
%              of the wight system pS, otherwise false (0). 
%  pkQ_wex  -- Returns true (1) if the kernel imputation of the excess game
%              is a null vector, otherwise false (0).
%  
%  input:
%  clv      -- TuGame class object.
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
%   07/24/2014        0.5             hme
%
%
if nargin<3
  N=clv.tusize;
  n=clv.tuplayers;
  np=length(pS);
  if np~=N
    error('Weight vector has not the correct size!');
  end
  tol=10^7*eps;
  x='';
elseif nargin<4
  N=clv.tusize;
  n=clv.tuplayers;
  np=length(pS);
  if np~=N
    error('Weight vector has not the correct size!');
  end
  tol=10^7*eps;
else
  N=clv.tusize;
  n=clv.tuplayers;
  np=length(pS);
  if np~=N
    error('Weight vector has not the correct size!');
  end
end

if isempty(x)
  wpk_v=clv.weightedKernel(pS);
else
  wpk_v=x;
end

wpkQ=clv.weightedKernelQ(wpk_v,pS);
wex=clv.weightedExcess(wpk_v,pS);
zv=zeros(1,n);
pkQ_wex=kernelQ(wex,zv);
cQ=pkQ_wex & wpkQ;
