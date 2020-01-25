function x=Talmudic_Rule(E,cl_vec,tol)
% TALMUDIC_RULE computes for a bankruptcy situation (E,cl_vec)
% the corresponding generalized contested garment or n-creditor solution.
% This solution coincides with the nucleolus of the corresponding
% modest bankruptcy game.
%
% Usage: x=Talmudic_Rule(E,cl_vec,tol)
% Define variables:
%  output:
%  x        -- The contested garment solution of the bankruptcy 
%              problem (E,cl_vec).
%  input:
%  E        -- Estate E (positive number).
%  cl_vec   -- Vector of claims of the claimants.
%  tol      -- Tolerance value. By default, it is set to 10^3*eps.
%              (optional) 



%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/17/2010        0.1 beta        hme
%                

if nargin<3
 tol=10^6*eps;
 else
end

n=length(cl_vec);
min_cl=min(cl_vec);
th_cl=n/2*min_cl;

if th_cl>E
   x=ones(1,n)*E/n;
else
sum_cl=cl_vec*ones(n,1);
hlf_cl=(1/2)*sum_cl;
lam=E/n;
x=ones(1,n);
x=find_sol(E,cl_vec,x,lam,hlf_cl,n,tol);
end

%-----------------------
function x=find_sol(E,cl_vec,x,lam,hlf_cl,n,tol);

if hlf_cl>E
   x=min(lam,(1/2)*cl_vec);
   lam=(E+max(0,lam-(1/2)*cl_vec)*ones(n,1))/n;
 else
   x=max(cl_vec-lam,(1/2)*cl_vec);
   lam=(max(cl_vec,lam+(1/2)*cl_vec)*ones(n,1)-E)/n;
end

sum_x=x*ones(n,1);

 if abs(sum_x-E)<tol
   x;
 else
%   warning('Solution not consistent!');
   x=find_sol(E,cl_vec,x,lam,hlf_cl,n,tol);
 end

