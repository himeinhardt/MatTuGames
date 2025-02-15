function ADC=Anti_CPCore(clv,x,tol)
% ANTI_CPCORE computes the closest point of the anti-core to x. Default is cplexmex
% otherwise Matlab's Optimization Toolbox. 
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.5.1 and higher) 
%
% Usage: ADC=clv.Anti_CPCore(x,tol)
% Define variables:
%  output:
%  ACp       -- Closest point of the anti-core to x.
%  D         -- The distance between the points.
%  acr_vaild -- Indicates if ACp is in the anti-core (1 otherwise 0). 
%  resid     -- The residual.
%  ef        -- Exitflag.
%  lambda    -- Containing the Lagrange multipliers.
%  x         -- The reference point from which the distance to core should be drawn.
%
%  input:
%  clv      -- TuGame class object.
%  x        -- The reference point from which the distance to core should be drawn.
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/21/2021        1.9.1           hme
%                


if nargin<3
 tol=10^8*eps;
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

S=1:N;
A=zeros(N,n);
for k=1:n, A(:,k) = bitget(S,k)==1;end
v1(S)=v(S)+tol;
v1(N)=v(N);
%A=-A;
B=v1';
Aeq=ones(1,n);
beq=v(N);
x=x';
try
    warning('off','all');
    options = cplexoptimset('cplex');
    options.barrier.convergetol=1e-12;
    options.simplex.tolerances.feasibility=1e-9;
    options.simplex.tolerances.optimality=1e-9;
    options.emphasis.numerical=0;
    options.threads=8;
    options.barrier.display=0;
    options.feasopt.tolerance=1e-12;
    warning('on','all');	
   [Cp,D,residual,exitflag,output,lambda] = cplexlsqlin(eye(n),x,A,B,Aeq,beq,[],[],[],options);
catch
    opts = optimoptions('lsqlin','Algorithm','active-set','Display','off','TolFun',1e-12);
    ws = optimwarmstart(x,opts);    
   [sol,d,residual,exitflag,output,lambda] = lsqlin(eye(n),x,A,B,Aeq,beq,[],[],ws);
   Cp=sol.X;
   D=d;   
end
abcQ=clv.belongToAntiCoreQ(v,Cp,'rat',tol);
Cp=Cp';
resid=residual';
ADC=struct('ACp',Cp,'D',D,'acr_vaild',abcQ,'resid',resid,'ef',exitflag,'lambda',lambda,'x',x');
