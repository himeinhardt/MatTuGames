function COV=p_COV_propertyQ(clv,x,m,t,str)
% P_COV_PROPERTYQ verifies if the payoff x satisfies the covariance
% with strategic equivalence property w.r.t. (v,m,t) using Matlab's PCT.
% 
%  Usage: COV=clv.p_COV_propertyQ(x,m,t,str)
%
% Define variables:
%  output: Fields
%  covQ     -- Returns true (1), if x satisfies COV, otherwise
%              false (0).
%  sol_v2   -- Solution of game v2 w.r.t. the string variable
%              'str'. For permissible strings see below.
%  sgm      -- Strategic equivalent solution.
%  v2       -- Strategic equivalent game.
%  x        -- Input payoff.
%
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n). Must be efficient.
%  m        -- An arbitrary vector of size(1,n).
%  t        -- A scalar s.t. t>0.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'PRN' that is, checking covariance with strategic
%               equivalence in accordance with the pre-nucleolus.
%              'PRK' that is, checking covariance with strategic
%               equivalence in accordance with the pre-kernel.
%              'SHAP' that is, checking covariance with strategic
%               equivalence in accordance with the Shapley value.
%              'MODIC' that is, checking covariance with strategic
%               equivalence in accordance with the modiclus.
%              Default is 'PRK'.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%              (optional) 
%              

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/18/2015        0.7             hme
%   02/10/2018        0.9             hme
%


N=clv.tusize;
n=clv.tuplayers;
v=clv.tuvalues;
gt=clv.tutype;
stx=clv.tustpt;

if nargin < 2
    str='PRK' ;
    t=1;
    m=ones(1,n);
    if isempty(stx)
       x=clv.p_PreKernel;
    end
elseif nargin < 3
    m=ones(1,n);
    t=1;
    str='PRK';
    if isempty(x)
       x=clv.p_PreKernel;
    end
elseif nargin < 4
    t=1;
    str='PRK';
    if isempty(x)
       x=clv.p_PreKernel;
    end
    if isempty(m)
       m=ones(1,n);
    end
elseif nargin<5
    str='PRK';
    if isempty(x)
       x=clv.p_PreKernel;
    end
    if isempty(m)
       m=ones(1,n);
    end
    if isempty(t)
       t=1;
    end
else
    if isempty(x)
       x=clv.p_PreKernel;
    end
    if isempty(m)
       m=ones(1,n);
    end
    if isempty(t)
       t=1;
    end
end


em=additive_game(m);
v2=t*v+em;

if strcmp('PRK',str)
   sgm=t*x + m;
   if PrekernelQ(v2,sgm)
      sol_v2=sgm;
   else 
       sol_v2=p_PreKernel(v2,sgm);
   end
   covQ=all(abs(sgm-sol_v2)<10^6*eps);
elseif strcmp('PRN',str)
   try 
       sol_v2=cplex_prenucl_llp(v2);
   catch
       sol_v2=Prenucl_llp(v2);
   end
   sgm=t*x + m;
   covQ=all(abs(sgm-sol_v2)<10^6*eps);
elseif strcmp('SHAP',str)
   sol_v2=p_ShapleyValue(v2);
   sgm=t*x + m;
   covQ=all(abs(sgm-sol_v2)<10^6*eps);
elseif strcmp('MODIC',str)
   try
       sol_v2=cplex_modiclus(v2);
   catch
       sol_v2=Modiclus(v2);
   end
   sgm=t*x + m;
   covQ=all(abs(sgm-sol_v2)<10^6*eps);
else
   sgm=t*x + m;
   if PrekernelQ(v2,sgm)
      sol_v2=sgm;
   else 
       sol_v2=p_PreKernel(v2,sgm);
   end
   covQ=all(abs(sgm-sol_v2)<10^6*eps);
end    

COV=struct('covQ',covQ,'sol_v2',sol_v2,'sgm',sgm,'v2',v2,'x',x);
