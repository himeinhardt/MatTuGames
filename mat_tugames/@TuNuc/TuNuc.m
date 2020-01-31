classdef TuNuc < TuSol
% TUNUC creates the subclass object TuNuc to perform several computations for retrieving 
% and modifying game data. It stores relevant game information and
% the nucleolus obtained by overloading functions from various solvers.
%
% Usage: clv = TuNuc(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- TuNuc class object (subclass of TuGame).
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
% TuNuc properties:
%
%  tu_nc_cdd         -- stores the nucleolus obtained with the CDD interface.
%  tu_nc_cplex       -- stores the nucleolus obtained with the CPLEX interface.
%  tu_nc_glpk        -- stores the nucleolus obtained with the GLPK interface.
%  tu_nc_gurobi      -- stores the nucleolus obtained with the GUROBI interface.
%  tu_nc_lp          -- stores the nucleolus obtained with MATLAB's linprog.
%  tu_nc_msk         -- stores the nucleolus element obtained with the MOSEK interface.
%  tu_tol            -- stores tolerance value. Default is 10^6*eps.
%  cdd_nc_valid      -- returns 1 if tu_pn_cdd stores the nucleolus, otherwise 0.
%  cplex_nc_valid    -- returns 1 if tu_pn_cplex stores the nucleolus, otherwise 0.
%  glpk_nc_valid     -- returns 1 if tu_pn_glpk stores the nucleolus, otherwise 0.
%  gurobi_nc_valid   -- returns 1 if tu_pn_gurobi stores the nucleolus, otherwise 0.
%  lp_nc_valid       -- returns 1 if tu_pn_lp stores the nucleolus, otherwise 0.
%  msk_nc_valid      -- returns 1 if tu_pn_msk stores the nucleolus, otherwise 0.
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
% TuNuc methods:
%  TuNuc                       -- creates the class object TuNuc.
%  setAllVendorNuc             -- sets a set of vendor pre-nucleolus solutions to TuNuc.
%  setCddNuc                   -- sets the computed pre-nucleolus derived with the interface CDD to TuNuc. 
%  setCplexNuc                 -- sets the computed pre-nucleolus derived with the interface CPLEX to TuNuc. 
%  setGlpkNuc                  -- sets the computed pre-nucleolus derived with the interface GLPK to TuNuc. 
%  setGurobiNuc                -- sets the computed pre-nucleolus derived with the interface GUROBI to TuNuc. 
%  setLpNuc                    -- sets the computed pre-nucleolus derived with MATLAB's linprog to TuNuc. 
%  setMskNuc                   -- sets the computed pre-nucleolus derived with the interface MOSEK to TuNuc. 
%  copyTuSol                   -- copies all solutions from TuSol/p_TuSol to TuNuc.
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
%   02/07/2017        0.9             hme
%





    properties(SetObservable = true)
       tu_nc_cdd
       tu_nc_cplex
       tu_nc_glpk
       tu_nc_gurobi
       tu_nc_lp  
       tu_nc_msk 
       tu_tol=10^6*eps; 
    end
  

    properties(GetAccess = 'public', SetAccess = 'private')
       cdd_nc_valid = false; 
       cplex_nc_valid = false;
       glpk_nc_valid = false;     
       gurobi_nc_valid = false;  
       lp_nc_valid = false;
       msk_nc_valid = false;     
     end


      
      
       methods
         function obj = TuNuc(w,gtype,gformat)
       % TUNUC creates the subclass object TuNuc to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information and
       % the nucleolus needed by overloading functions from various solvers.
       %
       % Usage: clv = TuNuc(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- TuNuc class object (subclass of TuGame).
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
           eQ=EssentialQ(w);
           if eQ==0
              error('The game is not essential!')
           end
           obj = obj@TuSol(w,gtype,gformat);
         end

         function obj = setCddNuc(obj,sol) 
         % SETCDDNUC sets the nucleolus using CDD to the class object TuNuc.
         % 
         %  Usage: clv = setCddNuc(clv,sol)
         %
         %  output:
         %    clv       -- TuNuc class object.
         %
         %  input:
         %     clv      -- TuNuc class object.
         %     sol      -- the nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_nc_cdd)
                   obj.tu_nc_cdd = [];   % Calls setter set.tu_nc_cdd.
               end
             else
               obj.tu_nc_cdd = sol;      % Calls setter set.tu_nc_cdd.
             end
         end


         function obj = set.tu_nc_cdd(obj,pt)
            if isempty(pt)
               try 
                 sol = CddNucl(obj);
                 obj.tu_nc_cdd = sol;
                 tol = obj.tu_tol;
                 obj.cdd_nc_valid=B0_balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_nc_cdd = inf(1,n);
                 obj.cdd_nc_valid=false; 
               end    
            else
               obj.tu_nc_cdd = pt;
               tol = obj.tu_tol;
               obj.cdd_nc_valid=B0_balancedCollectionQ(obj,pt,tol);
            end
         end         
         
         
         
         function obj = setCplexNuc(obj,sol) 
         % SETCPLEXNUC sets the nucleolus using CPLEX to the class object TuNuc.
         % 
         %  Usage: clv = setCplexNuc(clv,sol)
         %
         %  output:
         %    clv       -- TuNuc class object.
         %
         %  input:
         %     clv      -- TuNuc class object.
         %     sol      -- the nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_nc_cplex)
                   obj.tu_nc_cplex = [];   % Calls setter set.tu_nc_cplex.
               end
             else
               obj.tu_nc_cplex = sol;      % Calls setter set.tu_nc_cplex.
             end
         end


         function obj = set.tu_nc_cplex(obj,pt)
            if isempty(pt)
               try 
                 sol = cplex_nucl(obj);
                 obj.tu_cn_cplex = sol;
                 tol = obj.tu_tol;
                 obj.cplex_nc_valid=B0_balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_nc_cplex = inf(1,n);
                 obj.cplex_nc_valid=false; 
               end    
            else
               obj.tu_nc_cplex = pt;
               tol = obj.tu_tol;
               obj.cplex_nc_valid=B0_balancedCollectionQ(obj,pt,tol);
            end
         end
         

         function obj = setGlpkNuc(obj,sol) 
         % SETGLPKNUC sets the nucleolus using GLPK to the class object TuNuc.
         % 
         %  Usage: clv = setGlpkNuc(clv,sol)
         %
         %  output:
         %    clv       -- TuNuc class object.
         %
         %  input:
         %     clv      -- TuNuc class object.
         %     sol      -- the nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_nc_glpk)
                   obj.tu_nc_glpk = [];   % Calls setter set.tu_nc_glpk.
               end
             else
               obj.tu_nc_glpk = sol;      % Calls setter set.tu_nc_glpk.
             end
         end


         function obj = set.tu_nc_glpk(obj,pt)
            if isempty(pt)
               try 
                 sol = glpk_nucl(obj);
                 obj.tu_nc_glpk = sol;
                 tol = obj.tu_tol;
                 obj.glpk_nc_valid=B0_balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_nc_glpk = inf(1,n);
                 obj.glpk_nc_valid=false; 
               end    
            else
               obj.tu_nc_glpk = pt;
               tol = obj.tu_tol;
               obj.glpk_nc_valid=B0_balancedCollectionQ(obj,pt,tol);
            end
         end                  
         
         
         function obj = setGurobiNuc(obj,sol) 
         % SETGUROBINUC sets the nucleolus using GUROBI to the class object TuNuc.
         % 
         %  Usage: clv = setGurobiNuc(clv,sol)
         %
         %  output:
         %    clv       -- TuNuc class object.
         %
         %  input:
         %     clv      -- TuNuc class object.
         %     sol      -- the nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_nc_gurobi)
                   obj.tu_nc_gurobi = [];   % Calls setter set.tu nc_gurobi.
               end
             else
               obj.tu_nc_gurobi = sol;      % Calls setter set.tu_nc_gurobi.
             end
         end


         function obj = set.tu_nc_gurobi(obj,pt)
            if isempty(pt)
               try 
                 sol = gurobi_nucl(obj);
                 obj.tu_nc_gurobi = sol;
                 tol = obj.tu_tol;
                 obj.gurobi_nc_valid=B0_balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_nc_gurobi = inf(1,n);
                 obj.gurobi_nc_valid=false; 
               end    
            else
               obj.tu_nc_gurobi = pt;
               tol = obj.tu_tol;
               obj.gurobi_nc_valid=B0_balancedCollectionQ(obj,pt,tol);
            end
         end                           
         

         function obj = setLpNuc(obj,sol) 
         % SETLPNUC sets a nucleolus using MATLAB's linprog to the class object TuNuc.
         % 
         %  Usage: clv = setLpNuc(clv,sol)
         %
         %  output:
         %    clv       -- TuNuc class object.
         %
         %  input:
         %     clv      -- TuNuc class object.
         %     sol      -- the nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_nc_lp)
                   obj.tu_nc_lp = [];   % Calls setter set.tu_nc_lp.
               end
             else
               obj.tu_nc_lp = sol;      % Calls setter set.tu_nc_lp.
             end
         end


         function obj = set.tu_nc_lp(obj,pt)
            if isempty(pt)
               try 
                 sol = nucl(obj);
                 obj.tu_nc_lp = sol;
                 tol = obj.tu_tol;
                 obj.lp_nc_valid=B0_balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_nc_lp = inf(1,n);
                 obj.lp_nc_valid=false; 
               end    
            else
               obj.tu_nc_lp = pt;
               tol = obj.tu_tol;
               obj.lp_nc_valid=B0_balancedCollectionQ(obj,pt,tol);
            end
         end                  
         

         function obj = setMskNuc(obj,sol) 
         % SETMSKPRN sets a nucleolus using MOSEK to the class object TuNuc.
         % 
         %  Usage: clv = setMskNuc(clv,sol)
         %
         %  output:
         %    clv       -- TuNuc class object.
         %
         %  input:
         %     clv      -- TuNuc class object.
         %     sol      -- the nucleolus (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_nc_msk)
                   obj.tu_nc_msk = [];   % Calls setter set.tu_nc_msk.
               end
             else
               obj.tu_nc_msk = sol;      % Calls setter set.tu_nc_msk.
             end
         end


         function obj = set.tu_nc_msk(obj,pt)
            if isempty(pt)
               try 
                 sol = msk_nucl(obj);
                 obj.tu_nc_msk = sol;
                 tol = obj.tu_tol;
                 obj.msk_nc_valid=B0_balancedCollectionQ(obj,sol,tol);
               catch
                 n=obj.tuplayers;
                 obj.tu_nc_msk = inf(1,n);
                 obj.msk_nc_valid=false; 
               end    
            else
               obj.tu_nc_msk = pt;
               tol = obj.tu_tol;
               obj.msk_nc_valid=B0_balancedCollectionQ(obj,pt,tol);
            end
         end                  
         
         
         
         
         function obj = setAllVendorNuc(obj)
         % SETALLVENDORNUc sets a set of vendor nucleolus solutions to the class object TuNuc.
         %
         %  Usage: clv = setAllVendorNuc(clv)
         %
         %  output:
         %    clv       -- TuNuc class object.
         %
         %  input:
         %     clv      -- TuNuc class object.
         %
                %
                % CDD Nucleolus
                if obj.tuplayers < 13
                  try 
                    solnc = CddNucl(obj);
                    obj.tu_nc_cdd = solnc;
                    tol = obj.tu_tol;
                    obj.cdd_nc_valid=B0_balancedCollectionQ(obj,solnc,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_nc_cdd = inf(1,n);
                    obj.cdd_nc_valid=false; 
                  end
                end  
                %
                % CPLEX Nucleolus
                  try 
                    solnc = cplex_nucl(obj);
                    obj.tu_nc_cplex = solnc;
                    tol = obj.tu_tol;
                    obj.cplex_nc_valid=B0_balancedCollectionQ(obj,solnc,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_nc_cplex = inf(1,n);
                    obj.cplex_nc_valid=false;
                  end    
                %
                % GLPK Nucleolus
                if obj.tuplayers < 13
                  try 
                    solnc = glpk_nucl(obj);
                    obj.tu_nc_glpk = solnc;
                    tol = obj.tu_tol;
                    obj.glpk_nc_valid=B0_balancedCollectionQ(obj,solnc,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_nc_glpk = inf(1,n);
                    obj.glpk_nc_valid=false;
                  end
                end  
                %
                % GUROBI Nucleolus
                  try 
                    solnc = gurobi_nucl(obj);
                    obj.tu_nc_gurobi = solnc;
                    tol = obj.tu_tol;
                    obj.gurobi_nc_valid=B0_balancedCollectionQ(obj,solnc,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_nc_gurobi = inf(1,n);
                    obj.gurobi_nc_valid=false;
                  end
                %
                % MATLAB's linprog Nucleolus
                if obj.tuplayers < 13
                  try 
                    solnc = nucl(obj);
                    obj.tu_nc_lp = solnc;
                    tol = obj.tu_tol;
                    obj.lp_nc_valid=B0_balancedCollectionQ(obj,solnc,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_nc_lp = inf(1,n);
                    obj.lp_nc_valid=false;
                  end
                end  
                %
                % MOSEK Nucleolus
                if obj.tuplayers < 24
                  try 
                    solnc = msk_nucl(obj);
                    obj.tu_nc_msk = solnc;
                    tol = obj.tu_tol;
                    obj.msk_nc_valid=B0_balancedCollectionQ(obj,solnc,tol);
                  catch
                    n=obj.tuplayers;
                    obj.tu_nc_msk = inf(1,n);
                    obj.msk_nc_valid=false;
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
        %    obj       -- TuNuc class object.
        %
        %  input:
        %     obj      -- TuNuc class object.
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
