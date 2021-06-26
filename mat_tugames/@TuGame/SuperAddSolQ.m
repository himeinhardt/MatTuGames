function SUPAQ=SuperAddSolQ(clv,x,str,tol)
% SUPER_ADDSOLQ checks if the vector x is an element of a super additive solution
% of the game v using MPT3. Note, it checks the property SUPA w.r.t. to a particular bi-decomposition of v
% and not w.r.t. all bi-decompositions of v.
%
%  Usage: SUPAQ=clv.SuperAddSolQ(x,str,tol) 
%
%
% Define variables:
%  output:
%  SUPAQ    -- Returns true (1) if x is an element of a super additive solution,
%              for a particular bi-decompostion of v, otherwise false (0).
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n)
%  str      -- A string that defines different Methods.
%              Permissible methods are:
%              'CORE' in accordance with the core solution. Must be satisfied for all bi-decompositions of v, whenever the core exists.
%                     Hence, a false indicates that the provided solution cannot be correct.       
%              'PRK' in accordance with the pre-kernel. Can be satisfied for a particular bi-decompoition of v.
%              'PRN' in accordance with the pre-nucleolus. Can be satisfied for a particular bi-decompoition of v.
%              'SHAP' in accordance with the Shapley value. Must be satisfied for all bi-decompositions of v.
%                     Hence, a false indicates that the provided solution cannot be correct.
%              Default is 'CORE'.
%  tol      -- Tolerance value, its default value is 10^6*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/21/2020        1.9             hme
%

if nargin < 2
   [~,x]=clv.LeastCore();
   str='CORE';
   tol=10^6*eps;
elseif nargin < 3
   str='CORE';
   tol=10^6*eps;
   if isempty(x)==1
      if strcmp(str,'CORE')
        [~,x]=clv.LeastCore();
      elseif strcmp(str,'PRK')
        x=clv.PreKernel();
      elseif strcmp(str,'PRN')
         try
           x=clv.cplex_prenucl(); 
         catch
           x=clv.PreNucl()
         end
      elseif strcmp(str,'SHAP')
        x=clv.ShapleyValue();
      end
   end
elseif nargin < 4
   tol=10^6*eps;
   if isempty(x)==1
      if strcmp(str,'CORE')
        [~,x]=clv.LeastCore();
      elseif strcmp(str,'PRK')
        x=clv.PreKernel();
      elseif strcmp(str,'PRN')
         try
           x=clv.cplex_prenucl();
         catch
           x=clv.PreNucl()
         end
      elseif strcmp(str,'SHAP')
        x=clv.ShapleyValue();
      end
   end
end

N=clv.tusize;
n=clv.tuplayers;
SUPAQ=false;

if strcmp(str,'CORE') 
   crQ_v=clv.CddCoreQ(tol);
   if crQ_v==0
        msg01='No core exists!';
        warning('SupaQ:No1',msg01);
        return;
   else
     DecG_v=clv.DecomposeGame();
     vert_v=clv.CddCoreVertices();
     Pv = Polyhedron(vert_v);
     x=x';
     bcrQ=Pv.contains(x);
     if bcrQ==1
        vert_w=CddCoreVertices(DecG_v.w);
        Pw = Polyhedron(vert_w);
        vert_z=CddCoreVertices(DecG_v.z);
        crv_zw=(vert_z+vert_w)';
        SUPAQ=all(Pv.contains( crv_zw ));
     else
      msg02='Solution vector is not a core element!';
      warning('SupaQ:No2',msg02);
     end
   end
elseif strcmp(str,'SHAP')
   DecG_v=clv.DecomposeGame(); 
   sh_w=ShapleyValue(DecG_v.w);
   sh_z=ShapleyValue(DecG_v.z);
   sh_wz=sh_w+sh_z;
   SUPAQ=clv.ShapleyQ(sh_wz);
elseif strcmp(str,'PRK')
   DecG_v=clv.DecomposeGame();
   pk_w=PreKernel(DecG_v.w);
   pk_z=PreKernel(DecG_v.z);
   pk_wz=pk_w+pk_z;
   SUPAQ=clv.PrekernelQ(pk_wz,tol);
elseif strcmp(str,'PRN')
   DecG_v=clv.DecomposeGame();
   try
     pn_w=cplex_prenucl(DecG_v.w);
     pn_z=cplex_prenucl(DecG_v.z);
   catch
     pn_w=PreNucl(DecG_v.w);
     pn_z=PreNucl(DecG_v.z);
   end
   pn_wz=pn_w+pn_z;
   SUPAQ=clv.balancedCollectionQ(pn_wz,tol);
end
