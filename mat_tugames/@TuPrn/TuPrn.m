classdef TuPrn < TuSol
% TUPRN creates the subclass object TuPrn to perform several computations for retrieving 
% and modifying game data. It stores relevant game information and
% the pre-nucleolus obtained by overloading functions from various solvers.
%
% Usage: clv = TuPrn(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- TuPrn class object (subclass of TuGame).
%
% input:
% v             -- A Tu-Game v of length 2^n-1.
% gtype         -- A string to define the game type.
%                    Permissible types are:
%                    -- 'cv' (convex/average-convex, semi-convex).
%                    -- 'cr' game with non-empty core.
%                    -- 'sv' simple game.
%                    -- 'acr' game with non-empty anti-core.
%                    -- ''   empty string (default)
% gformat       -- A string to define the game format.
%                    Permissible formats are:
%                    -- 'mattug' i.e., unique integer representation to perform computation
%                                under MatTuGames. (default)
%                    -- 'mama'   i.e. generic power set representation, i.e Mathematica.
%
% TuPrn properties:
%
%  tu_pn_cdd         -- stores the pre-nucleolus obtained with the CDD interface.
%  tu_pn_cplex       -- stores the pre-nucleolus obtained with the CPLEX interface.
%  tu_pn_glpk        -- stores the pre-nucleolus obtained with the GLPK interface.
%  tu_pn_gurobi      -- stores the pre-nucleolus obtained with the GUROBI interface.
%  tu_pn_lp          -- stores the pre-nucleolus obtained with MATLAB's linprog.
%  tu_pn_msk         -- stores the pre-nucleolus element obtained with the MOSEK interface.
%  tu_tol            -- stores tolerance value. Default is 10^6*eps.
%  cdd_pn_valid      -- returns 1 if tu_pn_cdd stores the pre-nucleolus, otherwise 0.
%  cplex_pn_valid    -- returns 1 if tu_pn_cplex stores the pre-nucleolus, otherwise 0.
%  glpk_pn_valid     -- returns 1 if tu_pn_glpk stores the pre-nucleolus, otherwise 0.
%  gurobi_pn_valid   -- returns 1 if tu_pn_gurobi stores the pre-nucleolus, otherwise 0.
%  lp_pn_valid       -- returns 1 if tu_pn_lp stores the pre-nucleolus, otherwise 0.
%  msk_pn_valid      -- returns 1 if tu_pn_msk stores the pre-nucleolus, otherwise 0.
%
%  Properties inherited from the superclass TuSol:
%
%  tu_prk       -- stores a pre-kernel element.
%  tu_prk2      -- stores a second pre-kernel element instead of the pre-nucleolus.
%  tu_prn       -- stores the pre-nucleolus.
%  tu_sh        -- stores the Shapley value.
%  tu_tauv      -- stores the Tau value.
%  tu_bzf       -- stores the Banzhaf value.
%  tu_aprk      -- stores the anti-pre-kernel.
%  prk_valid    -- returns 1 if tu_prk stores a pre-kernel element, otherwise 0.
%  prk2_valid   -- returns 1 if tu_prk2 stores a pre-kernel element, otherwise 0.
%  prn_valid    -- returns 1 if tu_prn stores the pre-nucleolus, otherwise 0.
%
%  Properties inherited from the superclass TuGame:
%
%  tuvalues     -- stores the characteristic values of a Tu-game.
%  tusize       -- stores the length of the game array/vector.
%  tuplayers    -- stores the number of players involved.
%  tutype       -- stores the game type information ('cv','cr', or 'sv').
%  tuessQ       -- stores if the game is essential.
%  tuformat     -- stores the format how the game is represented.
%  tumv         -- stores the value of the grand coalition.
%  tumnQ        -- stores the information whether a proper coalition has a higher value
%                  than the grand coalition
%  tuSi         -- stores the coalitions having size of n-1.
%  tuvi         -- stores the values of singleton coalitions.
%  tustpt       -- stores a starting point for doing computation. Has lower priority than
%                  providing a second input argument as an alternative starting point.
%                  Thus, if a starting point is provided as second input argument the
%                  information stored in tustpt will not be used.
%
% TuPrn methods:
%  TuPrn                       -- creates the class object TuPrn.
%  setAllVendorPrn             -- sets a set of vendor pre-nucleolus solutions to TuPrn.
%  setCddPrn                   -- sets the computed pre-nucleolus derived with the interface CDD to TuPrn. 
%  setCplexPrn                 -- sets the computed pre-nucleolus derived with the interface CPLEX to TuPrn. 
%  setGlpkPrn                  -- sets the computed pre-nucleolus derived with the interface GLPK to TuPrn. 
%  setGurobiPrn                -- sets the computed pre-nucleolus derived with the interface GUROBI to TuPrn. 
%  setLpPrn                    -- sets the computed pre-nucleolus derived with MATLAB's linprog to TuPrn. 
%  setMskPrn                   -- sets the computed pre-nucleolus derived with the interface MOSEK to TuPrn. 
%  copyTuSol                   -- copies all solutions from TuSol/p_TuSol to TuPrn.
%
%  Methods inherited from the superclass TuSol:
%
%  setAllSolutions  -- sets all solutions listed below to the class object TuSol.
%  setPreKernel     -- sets a pre-kernel element to the class object TuSol.
%  setPreNuc        -- sets the pre-nucleolus to the class object TuSol.
%  setShapley       -- sets the Shapley value to the class object TuSol.
%  setTauValue      -- sets the Tau value to the class object TuSol.
%  setBanzhaf       -- sets the Banzhaf value to the class object TuSol.
%  setAntiPreKernel -- sets an anti-pre-kernel element to the class object TuSol.
%
%  Methods inherited from the superclass TuGame:
%
%  startpt      -- sets a starting point for doing computation.
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/12/2013        0.5             hme
%   11/04/2016        0.9             hme
%   05/15/2019        1.1             hme
%





    properties(SetObservable = true)
       tu_pn_cdd
       tu_pn_cplex
       tu_pn_glpk
       tu_pn_gurobi
       tu_pn_lp  
       tu_pn_msk 
       tu_tol=10^6*eps; 
    end
  

    properties(GetAccess = 'public', SetAccess = 'private')
       cdd_pn_valid = false; 
       cplex_pn_valid = false;
       glpk_pn_valid = false;     
       gurobi_pn_valid = false;  
       lp_pn_valid = false;
       msk_pn_valid = false;     
     end


      
      
       methods
         function obj = TuPrn(w,gtype,gformat)
       % TUPRN creates the subclass object TuPrn to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information and
       % the pre-nucleolus needed by overloading functions from various solvers.
       %
       % Usage: clv = TuPrn(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- TuPrn class object (subclass of TuGame).
       %
       % input:
       % v             -- A Tu-Game v of length 2^n-1.
       % gtype         -- A string to define the game type.
       %                    Permissible types are:
       %                    -- 'cv' (convex/average-convex, semi-convex).
       %                    -- 'cr' game with non-empty core.
       %                    -- 'sv' simple game.
       %                    -- 'acr' game with non-empty anti-core.
       %                    -- ''   empty string (default)
       % gformat       -- A string to define the game format.
       %                    Permissible formats are:
       %                    -- 'mattug' i.e., unique integer representation to perform computation
       %                                 under MatTuGames. (default)
       %                    -- 'mama'   i.e. generic power set representation, i.e Mathematica.
       %
           if nargin > 3
              error('Too many input arguments');
           elseif nargin < 1
              error('Game information must be given as a 2^n-1 vector!');
           elseif nargin < 2 
              gtype = '';
              gformat = 'mattug';
           elseif nargin < 3
              gformat = 'mattug';
           else
           end
           obj = obj@TuSol(w,gtype,gformat); 
         end

         function obj = setCddPrn(obj,sol) 
         % SETCDDPRN sets the pre-nucleolus using CDD to the class object TuPrn.
         % 
         %  Usage: clv = setCddPrn(clv,sol)
         %
         %  output:
         %    clv       -- TuPrn class object.
         %
         %  input:
         %     clv      -- TuPrn class object.
         %     sol      -- the pre-nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_pn_cdd)
                   obj.tu_pn_cdd = [];   % Calls setter set.tu_pn_cdd.
               end
             else
               obj.tu_pn_cdd = sol;      % Calls setter set.tu_pn_cdd.
             end
         end


         function obj = set.tu_pn_cdd(obj,pt)
            if isempty(pt)
               try 
                 sol = CddPrenucl(obj);
                 obj.tu_pn_cdd = sol;
                 tol = obj.tu_tol;
                 obj.cdd_pn_valid=balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_pn_cdd = inf(1,n);
                 obj.cdd_pn_valid=false; 
               end    
            else
               obj.tu_pn_cdd = pt;
               tol = obj.tu_tol;
               obj.cdd_pn_valid=balancedCollectionQ(obj,pt,tol);
            end
         end         
         
         
         
         function obj = setCplexPrn(obj,sol) 
         % SETCPLEXPRN sets the pre-nucleolus using CPLEX to the class object TuPrn.
         % 
         %  Usage: clv = setCplexPrn(clv,sol)
         %
         %  output:
         %    clv       -- TuPrn class object.
         %
         %  input:
         %     clv      -- TuPrn class object.
         %     sol      -- the pre-nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_pn_cplex)
                   obj.tu_pn_cplex = [];   % Calls setter set.tu_pn_cplex.
               end
             else
               obj.tu_pn_cplex = sol;      % Calls setter set.tu_pn_cplex.
             end
         end


         function obj = set.tu_pn_cplex(obj,pt)
            if isempty(pt)
               try 
                 sol = cplex_prenucl(obj);
                 obj.tu_pn_cplex = sol;
                 tol = obj.tu_tol;
                 obj.cplex_pn_valid=balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_pn_cplex = inf(1,n);
                 obj.cplex_pn_valid=false; 
               end    
            else
               obj.tu_pn_cplex = pt;
               tol = obj.tu_tol;
               obj.cplex_pn_valid=balancedCollectionQ(obj,pt,tol);
            end
         end
         

         function obj = setGlpkPrn(obj,sol) 
         % SETGLPKPRN sets the pre-nucleolus using GLPK to the class object TuPrn.
         % 
         %  Usage: clv = setGlpkPrn(clv,sol)
         %
         %  output:
         %    clv       -- TuPrn class object.
         %
         %  input:
         %     clv      -- TuPrn class object.
         %     sol      -- the pre-nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_pn_glpk)
                   obj.tu_pn_glpk = [];   % Calls setter set.tu_pn_glpk.
               end
             else
               obj.tu_pn_glpk = sol;      % Calls setter set.tu_pn_glpk.
             end
         end


         function obj = set.tu_pn_glpk(obj,pt)
            if isempty(pt)
               try 
                 sol = glpk_prenucl(obj);
                 obj.tu_pn_glpk = sol;
                 tol = obj.tu_tol;
                 obj.glpk_pn_valid=balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_pn_glpk = inf(1,n);
                 obj.glpk_pn_valid=false; 
               end    
            else
               obj.tu_pn_glpk = pt;
               tol = obj.tu_tol;
               obj.glpk_pn_valid=balancedCollectionQ(obj,pt,tol);
            end
         end                  
         
         
         function obj = setGurobiPrn(obj,sol) 
         % SETGUROBIPRN sets the pre-nucleolus using GUROBI to the class object TuPrn.
         % 
         %  Usage: clv = setGurobiPrn(clv,sol)
         %
         %  output:
         %    clv       -- TuPrn class object.
         %
         %  input:
         %     clv      -- TuPrn class object.
         %     sol      -- the pre-nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_pn_gurobi)
                   obj.tu_pn_gurobi = [];   % Calls setter set.tu_pn_gurobi.
               end
             else
               obj.tu_pn_gurobi = sol;      % Calls setter set.tu_pn_gurobi.
             end
         end


         function obj = set.tu_pn_gurobi(obj,pt)
            if isempty(pt)
               try 
                 sol = gurobi_prenucl(obj);
                 obj.tu_pn_gurobi = sol;
                 tol = obj.tu_tol;
                 obj.gurobi_pn_valid=balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_pn_gurobi = inf(1,n);
                 obj.gurobi_pn_valid=false; 
               end    
            else
               obj.tu_pn_gurobi = pt;
               tol = obj.tu_tol;
               obj.gurobi_pn_valid=balancedCollectionQ(obj,pt,tol);
            end
         end                           
         

         function obj = setLpPrn(obj,sol) 
         % SETLPPRN sets a pre-nucleolus using MATLAB's linprog to the class object TuPrn.
         % 
         %  Usage: clv = setLpPrn(clv,sol)
         %
         %  output:
         %    clv       -- TuPrn class object.
         %
         %  input:
         %     clv      -- TuPrn class object.
         %     sol      -- the pre-nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_pn_lp)
                   obj.tu_pn_lp = [];   % Calls setter set.tu_pn_lp.
               end
             else
               obj.tu_pn_lp = sol;      % Calls setter set.tu_pn_lp.
             end
         end


         function obj = set.tu_pn_lp(obj,pt)
            if isempty(pt)
               try 
                 sol = PreNucl(obj);
                 obj.tu_pn_lp = sol;
                 tol = obj.tu_tol;
                 obj.lp_pn_valid=balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_pn_lp = inf(1,n);
                 obj.lp_pn_valid=false; 
               end    
            else
               obj.tu_pn_lp = pt;
               tol = obj.tu_tol;
               obj.lp_pn_valid=balancedCollectionQ(obj,pt,tol);
            end
         end                  
         

         function obj = setMskPrn(obj,sol) 
         % SETMSKPRN sets a pre-nucleolus using MOSEK to the class object TuPrn.
         % 
         %  Usage: clv = setMskPrn(clv,sol)
         %
         %  output:
         %    clv       -- TuPrn class object.
         %
         %  input:
         %     clv      -- TuPrn class object.
         %     sol      -- the pre-nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_pn_msk)
                   obj.tu_pn_msk = [];   % Calls setter set.tu_pn_msk.
               end
             else
               obj.tu_pn_msk = sol;      % Calls setter set.tu_pn_msk.
             end
         end


         function obj = set.tu_pn_msk(obj,pt)
            if isempty(pt)
               try 
                 sol = msk_prenucl(obj);
                 obj.tu_pn_msk = sol;
                 tol = obj.tu_tol;
                 obj.msk_pn_valid=balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_pn_msk = inf(1,n);
                 obj.msk_pn_valid=false; 
               end    
            else
               obj.tu_pn_msk = pt;
               tol = obj.tu_tol;
               obj.msk_pn_valid=balancedCollectionQ(obj,pt,tol);
            end
         end                  
         
         
         
         
         function obj = setAllVendorPrn(obj)
         % SETALLVENDORPRN sets a set of vendor pre-nucleolus solutions to the class object TuPrn.
         %
         %  Usage: clv = setAllVendorPrn(clv)
         %
         %  output:
         %    clv       -- TuPrn class object.
         %
         %  input:
         %     clv      -- TuPrn class object.
         %
                %
                % CDD Pre-Nucleolus
                if obj.tuplayers < 13
                  try 
                    solpn = CddPrenucl(obj);
                    obj.tu_pn_cdd = solpn;
                    tol = obj.tu_tol;
                    obj.cdd_pn_valid=balancedCollectionQ(obj,solpn,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_pn_cdd = inf(1,n);
                    obj.cdd_pn_valid=false; 
                  end
                end  
                %
                % CPLEX Pre-Nucleolus
                  try 
                    solpn = cplex_prenucl(obj);
                    obj.tu_pn_cplex = solpn;
                    tol = obj.tu_tol;
                    obj.cplex_pn_valid=balancedCollectionQ(obj,solpn,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_pn_cplex = inf(1,n);
                    obj.cplex_pn_valid=false;
                  end    
                %
                % GLPK Pre-Nucleolus
                if obj.tuplayers < 13
                  try 
                    solpn = glpk_prenucl(obj);
                    obj.tu_pn_glpk = solpn;
                    tol = obj.tu_tol;
                    obj.glpk_pn_valid=balancedCollectionQ(obj,solpn,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_pn_glpk = inf(1,n);
                    obj.glpk_pn_valid=false;
                  end
                end  
                %
                % GUROBI Pre-Nucleolus
                  try 
                    solpn = gurobi_prenucl(obj);
                    obj.tu_pn_gurobi = solpn;
                    tol = obj.tu_tol;
                    obj.gurobi_pn_valid=balancedCollectionQ(obj,solpn,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_pn_gurobi = inf(1,n);
                    obj.gurobi_pn_valid=false;
                  end
                %
                % MATLAB's linprog Pre-Nucleolus
                if obj.tuplayers < 13
                  try 
                    solpn = PreNucl(obj);
                    obj.tu_pn_lp = solpn;
                    tol = obj.tu_tol;
                    obj.lp_pn_valid=balancedCollectionQ(obj,solpn,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_pn_lp = inf(1,n);
                    obj.lp_pn_valid=false;
                  end
                end  
                %
                % MOSEK Pre-Nucleolus
                if obj.tuplayers < 24
                  try 
                    solpn = msk_prenucl(obj);
                    obj.tu_pn_msk = solpn;
                    tol = obj.tu_tol;
                    obj.msk_pn_valid=balancedCollectionQ(obj,solpn,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_pn_msk = inf(1,n);
                    obj.msk_pn_valid=false;
                  end
                end  
                %
                
        end
         
         


        function obj = copyTuSol(obj,clv)
        % copy constructor: class object clv to class object obj.
        %
        %  Usage: obj = copyTuSol(obj,clv)
        %
        %  output:
        %    obj       -- TuPrn class object.
        %
        %  input:
        %     obj      -- TuPrn class object.
        %     clv      -- TuSol/p_TuSol class object.
        %
          fns = properties(clv);
             for i=1:7
              try 
                 obj.(fns{i}) = clv.(fns{i});
              catch
                 if strcmp(obj.tutype,'sv') 
                    obj.(fns{i}) = clv.(fns{i});
                 end
              end
             end
        end


       end

 
end
