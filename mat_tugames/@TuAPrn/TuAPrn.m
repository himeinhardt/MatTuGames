classdef TuAPrn < TuSol
% TUAPRN creates the subclass object TuAPrn to perform several computations for retrieving 
% and modifying game data. It stores relevant game information and
% the anti pre-nucleolus obtained by overloading functions from various solvers.
%
% Usage: clv = TuAPrn(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- TuAPrn class object (subclass of TuGame).
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
%  tu_apn_cdd         -- stores the anti pre-nucleolus obtained with the CDD interface.
%  tu_apn_cplex       -- stores the anti pre-nucleolus obtained with the CPLEX interface.
%  tu_apn_glpk        -- stores the anti pre-nucleolus obtained with the GLPK interface.
%  tu_apn_gurobi      -- stores the anti pre-nucleolus obtained with the GUROBI interface.
%  tu_apn_lp          -- stores the anti pre-nucleolus obtained with MATLAB's linprog.
%  tu_apn_msk         -- stores the anti pre-nucleolus element obtained with the MOSEK interface.
%  tu_tol            -- stores tolerance value. Default is 10^6*eps.
%  cdd_apn_valid      -- returns 1 if tu_apn_cdd stores the anti pre-nucleolus, otherwise 0.
%  cplex_apn_valid    -- returns 1 if tu_apn_cplex stores the anti pre-nucleolus, otherwise 0.
%  glpk_apn_valid     -- returns 1 if tu_apn_glpk stores the anti pre-nucleolus, otherwise 0.
%  gurobi_apn_valid   -- returns 1 if tu_apn_gurobi stores the anti pre-nucleolus, otherwise 0.
%  lp_apn_valid       -- returns 1 if tu_apn_lp stores the anti pre-nucleolus, otherwise 0.
%  msk_apn_valid      -- returns 1 if tu_apn_msk stores the anti pre-nucleolus, otherwise 0.
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
% TuAPrn methods:
%  TuAPrn                       -- creates the class object TuAPrn.
%  setAllVendorAPrn             -- sets a set of vendor anti pre-nucleolus solutions to TuAPrn.
%  setCddAPrn                   -- sets the computed anti pre-nucleolus derived with the interface CDD to TuAPrn. 
%  setCplexAPrn                 -- sets the computed anti pre-nucleolus derived with the interface CPLEX to TuAPrn. 
%  setGlpkAPrn                  -- sets the computed anti pre-nucleolus derived with the interface GLPK to TuAPrn. 
%  setGurobiAPrn                -- sets the computed anti pre-nucleolus derived with the interface GUROBI to TuAPrn. 
%  setLpAPrn                    -- sets the computed anti pre-nucleolus derived with MATLAB's linprog to TuAPrn. 
%  setMskAPrn                   -- sets the computed anti pre-nucleolus derived with the interface MOSEK to TuAPrn. 
%  copyTuSol                    -- copies all solutions from TuSol/p_TuSol to TuAPrn.
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
%   08/09/2016        0.9             hme
%





    properties(SetObservable = true)
       tu_apn_cdd
       tu_apn_cplex
       tu_apn_glpk
       tu_apn_gurobi
       tu_apn_lp  
       tu_apn_msk 
       tu_tol=10^6*eps; 
    end
  

    properties(GetAccess = 'public', SetAccess = 'private')
       cdd_apn_valid = false; 
       cplex_apn_valid = false;
       glpk_apn_valid = false;     
       gurobi_apn_valid = false;  
       lp_apn_valid = false;
       msk_apn_valid = false;     
     end


      
      
       methods
         function obj = TuAPrn(w,gtype,gformat)
       % TUAPRN creates the subclass object TuAPrn to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information and
       % the pre-nucleolus needed by overloading functions from various solvers.
       %
       % Usage: clv = TuAPrn(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- TuAPrn class object (subclass of TuGame).
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

         function obj = setCddAPrn(obj,sol) 
         % SETCDDAPRN sets the anti pre-nucleolus using CDD to the class object TuAPrn.
         % 
         %  Usage: clv = setCddAPrn(clv,sol)
         %
         %  output:
         %    clv       -- TuAPrn class object.
         %
         %  input:
         %     clv      -- TuAPrn class object.
         %     sol      -- the anti pre-nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_apn_cdd)
                   obj.tu_apn_cdd = [];   % Calls setter set.tu_pn_cdd.
               end
             else
               obj.tu_apn_cdd = sol;      % Calls setter set.tu_pn_cdd.
             end
         end


         function obj = set.tu_apn_cdd(obj,pt)
            if isempty(pt)
               try 
                 sol = CddAntiPrenucl(obj);
                 obj.tu_apn_cdd = sol;
                 tol = obj.tu_tol;
                 obj.cdd_apn_valid=Anti_balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_apn_cdd = inf(1,n); 
               end    
            else
               obj.tu_apn_cdd = pt;
               tol = obj.tu_tol;
               obj.cdd_apn_valid=Anti_balancedCollectionQ(obj,pt,tol);
            end
         end         
         
         
         
         function obj = setCplexAPrn(obj,sol) 
         % SETCPLEXAPRN sets the anti pre-nucleolus using CPLEX to the class object TuAPrn.
         % 
         %  Usage: clv = setCplexAPrn(clv,sol)
         %
         %  output:
         %    clv       -- TuAPrn class object.
         %
         %  input:
         %     clv      -- TuAPrn class object.
         %     sol      -- the anti pre-nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_apn_cplex)
                   obj.tu_apn_cplex = [];   % Calls setter set.tu_pn_cplex.
               end
             else
               obj.tu_apn_cplex = sol;      % Calls setter set.tu_pn_cplex.
             end
         end


         function obj = set.tu_apn_cplex(obj,pt)
            if isempty(pt)
               try 
                 sol = cplex_AntiPreNucl(obj);
                 obj.tu_apn_cplex = sol;
                 tol = obj.tu_tol;
                 obj.cplex_apn_valid=Anti_balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_apn_cplex = inf(1,n); 
               end    
            else
               obj.tu_apn_cplex = pt;
               tol = obj.tu_tol;
               obj.cplex_apn_valid=Anti_balancedCollectionQ(obj,pt,tol);
            end
         end
         

         function obj = setGlpkAPrn(obj,sol) 
         % SETGLPKAPRN sets the anti pre-nucleolus using GLPK to the class object TuAPrn.
         % 
         %  Usage: clv = setGlpkAPrn(clv,sol)
         %
         %  output:
         %    clv       -- TuAPrn class object.
         %
         %  input:
         %     clv      -- TuAPrn class object.
         %     sol      -- the anti pre-nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_apn_glpk)
                   obj.tu_apn_glpk = [];   % Calls setter set.tu_pn_glpk.
               end
             else
               obj.tu_apn_glpk = sol;      % Calls setter set.tu_pn_glpk.
             end
         end


         function obj = set.tu_apn_glpk(obj,pt)
            if isempty(pt)
               try 
                 sol = glpk_AntiPreNucl(obj);
                 obj.tu_apn_glpk = sol;
                 tol = obj.tu_tol;
                 obj.glpk_apn_valid=Anti_balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_apn_glpk = inf(1,n); 
               end    
            else
               obj.tu_apn_glpk = pt;
               tol = obj.tu_tol;
               obj.glpk_apn_valid=Anti_balancedCollectionQ(obj,pt,tol);
            end
         end                  
         
         
         function obj = setGurobiAPrn(obj,sol) 
         % SETGUROBIAPRN sets the anti pre-nucleolus using GUROBI to the class object TuAPrn.
         % 
         %  Usage: clv = setGurobiAPrn(clv,sol)
         %
         %  output:
         %    clv       -- TuAPrn class object.
         %
         %  input:
         %     clv      -- TuAPrn class object.
         %     sol      -- the anti pre-nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_apn_gurobi)
                   obj.tu_apn_gurobi = [];   % Calls setter set.tu_pn_gurobi.
               end
             else
               obj.tu_apn_gurobi = sol;      % Calls setter set.tu_pn_gurobi.
             end
         end


         function obj = set.tu_apn_gurobi(obj,pt)
            if isempty(pt)
               try 
                 sol = gurobi_AntiPreNucl(obj);
                 obj.tu_apn_gurobi = sol;
                 tol = obj.tu_tol;
                 obj.gurobi_apn_valid=Anti_balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_apn_gurobi = inf(1,n); 
               end    
            else
               obj.tu_apn_gurobi = pt;
               tol = obj.tu_tol;
               obj.gurobi_apn_valid=Anti_balancedCollectionQ(obj,pt,tol);
            end
         end                           
         

         function obj = setLpAPrn(obj,sol) 
         % SETLPAPRN sets a anti pre-nucleolus using MATLAB's linprog to the class object TuAPrn.
         % 
         %  Usage: clv = setLpAPrn(clv,sol)
         %
         %  output:
         %    clv       -- TuAPrn class object.
         %
         %  input:
         %     clv      -- TuAPrn class object.
         %     sol      -- the anti pre-nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_apn_lp)
                   obj.tu_apn_lp = [];   % Calls setter set.tu_pn_lp.
               end
             else
               obj.tu_apn_lp = sol;      % Calls setter set.tu_pn_lp.
             end
         end


         function obj = set.tu_apn_lp(obj,pt)
            if isempty(pt)
               try 
                 sol = Anti_PreNucl(obj);
                 obj.tu_apn_lp = sol;
                 tol = obj.tu_tol;
                 obj.lp_apn_valid=Anti_balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_apn_lp = inf(1,n); 
               end    
            else
               obj.tu_apn_lp = pt;
               tol = obj.tu_tol;
               obj.lp_apn_valid=Anti_balancedCollectionQ(obj,pt,tol);
            end
         end                  
         

         function obj = setMskAPrn(obj,sol) 
         % SETMSKAPRN sets a anti pre-nucleolus using MOSEK to the class object TuAPrn.
         % 
         %  Usage: clv = setMskAPrn(clv,sol)
         %
         %  output:
         %    clv       -- TuAPrn class object.
         %
         %  input:
         %     clv      -- TuAPrn class object.
         %     sol      -- the anti pre-nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_apn_msk)
                   obj.tu_apn_msk = [];   % Calls setter set.tu_pn_msk.
               end
             else
               obj.tu_apn_msk = sol;      % Calls setter set.tu_pn_msk.
             end
         end


         function obj = set.tu_apn_msk(obj,pt)
            if isempty(pt)
               try 
                 sol = msk_AntiPreNucl(obj);
                 obj.tu_apn_msk = sol;
                 tol = obj.tu_tol;
                 obj.msk_apn_valid=Anti_balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_apn_msk = inf(1,n); 
               end    
            else
               obj.tu_apn_msk = pt;
               tol = obj.tu_tol;
               obj.msk_apn_valid=Anti_balancedCollectionQ(obj,pt,tol);
            end
         end                  
         
         
         
         
         function obj = setAllVendorAPrn(obj)
         % SETALLVENDORAPRN sets a set of vendor anti pre-nucleolus solutions to the class object TuAPrn.
         %
         %  Usage: clv = setAllVendorAPrn(clv)
         %
         %  output:
         %    clv       -- TuAPrn class object.
         %
         %  input:
         %     clv      -- TuAPrn class object.
         %
                %
                % CDD Anti Pre-Nucleolus
                if obj.tuplayers < 13
                  try 
                    solpn = CddAntiPrenucl(obj);
                    obj.tu_apn_cdd = solpn;
                    tol = obj.tu_tol;
                    obj.cdd_apn_valid=Anti_balancedCollectionQ(obj,solpn,tol);
                  catch
                    obj.tu_apn_cdd = inf(1,n); 
                  end
                end  
                %
                % CPLEX Anti Pre-Nucleolus
                  try 
                    solpn = cplex_AntiPreNucl(obj);
                    obj.tu_apn_cplex = solpn;
                    tol = obj.tu_tol;
                    obj.cplex_apn_valid=Anti_balancedCollectionQ(obj,solpn,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_apn_cplex = inf(1,n);
                  end    
                %
                % GLPK Anti Pre-Nucleolus
                if obj.tuplayers < 13
                  try 
                    solpn = glpk_AntiPreNucl(obj);
                    obj.tu_apn_glpk = solpn;
                    tol = obj.tu_tol;
                    obj.glpk_apn_valid=Anti_balancedCollectionQ(obj,solpn,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_apn_glpk = inf(1,n);
                  end
                end  
                %
                % GUROBI Anti Pre-Nucleolus
                  try 
                    solpn = gurobi_AntiPreNucl(obj);
                    obj.tu_apn_gurobi = solpn;
                    tol = obj.tu_tol;
                    obj.gurobi_apn_valid=Anti_balancedCollectionQ(obj,solpn,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_apn_gurobi = inf(1,n);
                  end
                %
                % MATLAB's linprog Anti Pre-Nucleolus
                if obj.tuplayers < 13
                  try 
                    solpn = Anti_PreNucl(obj);
                    obj.tu_apn_lp = solpn;
                    tol = obj.tu_tol;
                    obj.lp_apn_valid=Anti_balancedCollectionQ(obj,solpn,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_apn_lp = inf(1,n);
                  end
                end  
                %
                % MOSEK Anti Pre-Nucleolus
                if obj.tuplayers < 24
                  try 
                    solpn = msk_AntiPreNucl(obj);
                    obj.tu_apn_msk = solpn;
                    tol = obj.tu_tol;
                    obj.msk_apn_valid=Anti_balancedCollectionQ(obj,solpn,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_apn_msk = inf(1,n);
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
