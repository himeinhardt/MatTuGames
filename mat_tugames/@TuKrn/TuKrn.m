classdef TuKrn < TuSol
% TUPRK creates the subclass object TuKrn to perform several computations for retrieving 
% and modifying game data. It stores relevant game information and
% kernel elements obtained by overloading functions from various solvers.
%
% Usage: clv = TuKrn(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- TuKrn class object (subclass of TuGame).
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
% TuKrn properties:
%
%  tu_kr_cplex       -- stores a kernel element obtained with the CPLEX interface.
%  tu_kr_cvx         -- stores a kernel element obtained with the CVX interface.
%  tu_kr_gurobi      -- stores a kernel element obtained with the GUROBI interface.
%  tu_kr_ipopt       -- stores a kernel element obtained with the IPOPT interface.
%  tu_kr_msk         -- stores a kernel element obtained with the MOSEK interface.
%  tu_kr_oases       -- stores a kernel element obtained with the OASES interface.
%  tu_kr_qpc         -- stores a kernel element obtained with the QPC interface.
%  tu_tol            -- stores the tolerance value. Default is set to 10^6*eps.
%  cplex_kr_valid    -- returns 1 if tu_kr_cplex stores a kernel element, otherwise 0.
%  cvx_kr_valid      -- returns 1 if tu_kr_cvx stores a kernel element, otherwise 0.
%  gurobi_kr_valid   -- returns 1 if tu_kr_gurobi stores a kernel element, otherwise 0.
%  ipopt_kr_valid    -- returns 1 if tu_kr_ipopt stores a kernel element, otherwise 0.
%  msk_kr_valid      -- returns 1 if tu_kr_msk stores a kernel element, otherwise 0.
%  oases_kr_valid    -- returns 1 if tu_kr_oases stores a kernel element, otherwise 0.
%  qpc_kr_valid      -- returns 1 if tu_kr_qpc stores a kernel element, otherwise 0.
%
%  Properties inherited from the superclass TuSol:
%
%  tu_prk       -- stores a pre-kernel element.
%  tu_prk2      -- stores a second pre-kernel element instead of the pre-nucleolus.
%  tu_prn       -- stores the pre-nucleolus.
%  tu_sh        -- stores the Shapley value.
%  tu_tauv      -- stores the Tau value.
%  tu_bzf       -- stores the Banzhaf value.
%  tu_aprk      -- stores the anti-kernel.
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
% TuKrn methods:
%  TuKrn                       -- creates the class object TuKrn.
%  setAllVendorKrn             -- sets all computed kernel elements to TuKrn.
%  setCplexKrn                 -- sets the computed kernel element dervied with the interface CPLEX to TuKrn. 
%  setCvxKrn                   -- sets the computed kernel element dervied with the interface CVX to TuKrn. 
%  setGurobiKrn                -- sets the computed kernel element dervied with the interface GUROBI to TuKrn. 
%  setHslKrn                   -- sets the computed kernel element dervied with the interface HSL to TuKrn. 
%  setIpoptKrn                 -- sets the computed kernel element dervied with the interface IPOPT to TuKrn. 
%  setLinKrn                   -- sets the computed kernel element dervied with MATLAB's linsolve to TuKrn. 
%  setMskKrn                   -- sets the computed kernel element dervied with the interface MOSEK to TuKrn. 
%  setOasesKrn                 -- sets the computed kernel element dervied with the interface OASES to TuKrn. 
%  setOlsKrn                   -- sets the computed kernel element dervied with MATLAB's lscov to TuKrn. 
%  copyTuSol                   -- copies all solutions from TuSol/p_TuSol to TuKrn.
%  copy_p_TuKrn                -- copies consistency properties from p_TuKrn to TuKrn.
%
%  Methods inherited from the superclass TuSol:
%
%  setAllSolutions  -- sets all solutions listed below to the class object TuSol.
%  setPreKernel     -- sets a pre-kernel element to the class object TuSol.
%  setPreNuc        -- sets the pre-nucleolus to the class object TuSol.
%  setShapley       -- sets the Shapley value to the class object TuSol.
%  setTauValue      -- sets the Tau value to the class object TuSol.
%  setBanzhaf       -- sets the Banzhaf value to the class object TuSol.
%  setAntiPreKernel -- sets an anti-kernel element to the class object TuSol.
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
       tu_kr_cplex
       tu_kr_cvx  
       tu_kr_gurobi
       tu_kr_ipopt
       tu_kr_msk  
       tu_kr_oases
       tu_kr_qpc
       tu_tol=10^6*eps; 
    end
  
    properties(GetAccess = 'public', SetAccess = 'private')
       cplex_kr_valid = false; 
       cvx_kr_valid = false;     
       gurobi_kr_valid = false;  
       ipopt_kr_valid = false;   
       msk_kr_valid = false;     
       oases_kr_valid = false;   
       qpc_kr_valid = false;
     end
      
      
       methods
         function obj = TuKrn(w,gtype,gformat)
       % TUPRK creates the subclass object TuKrn to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information and
       % kernel elements needed by overloading functions from various solvers.
       %
       % Usage: clv = TuKrn(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- TuKrn class object (subclass of TuGame).
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

         function obj = setCplexKrn(obj,sol) 
         % SETCPLEXKRN sets a kernel element using CPLEX to the class object TuKrn.
         % 
         %  Usage: clv = setCplexKrn(clv,sol)
         %
         %  output:
         %    clv       -- TuKrn class object.
         %
         %  input:
         %     clv      -- TuKrn class object.
         %     sol      -- a kernel element (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_kr_cplex)
                   obj.tu_kr_cplex = [];   % Calls setter set.tu_kr_cplex.
               end
             else
               obj.tu_kr_cplex = sol;      % Calls setter set.tu_kr_cplex.
             end
         end


         function obj = set.tu_kr_cplex(obj,pt)
            if isempty(pt)
               try 
                 sol = cplex_kernel(obj);
                 tol = obj.tu_tol;
                 pkQ = kernelQ(obj,sol,tol);
                 obj.tu_kr_cplex = sol;
                 obj.cplex_kr_valid = pkQ;
               catch
                 n=obj.tuplayers;
                 obj.tu_kr_cplex = inf(1,n);
                 obj.cplex_kr_valid = false;                 
               end    
            else
               obj.tu_kr_cplex = pt;
               tol = obj.tu_tol;
               obj.cplex_kr_valid = kernelQ(obj,pt,tol);
            end
         end
         

         function obj = setCvxKrn(obj,sol) 
         % SETCVXKRN sets a kernel element using CVX to the class object TuKrn.
         % 
         %  Usage: clv = setCvxKrn(clv,sol)
         %
         %  output:
         %    clv       -- TuKrn class object.
         %
         %  input:
         %     clv      -- TuKrn class object.
         %     sol      -- a kernel element (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_kr_cvx)
                   obj.tu_kr_cvx = [];   % Calls setter set.tu_kr_cvx.
               end
             else
               obj.tu_kr_cvx = sol;      % Calls setter set.tu_kr_cvx.
             end
         end


         function obj = set.tu_kr_cvx(obj,pt)
            if isempty(pt)
               try 
                 sol = cvx_kernel(obj);
                 tol = obj.tu_tol;
                 pkQ = kernelQ(obj,sol,tol);
                 obj.tu_kr_cvx = sol;
                 obj.cvx_kr_valid = pkQ;
               catch
                 n=obj.tuplayers;
                 obj.tu_kr_cvx = inf(1,n);
                 obj.cvx_kr_valid = false;                 
               end    
            else
               obj.tu_kr_cvx = pt;
               tol = obj.tu_tol;
               obj.cvx_kr_valid = kernelQ(obj,pt,tol);
            end
         end         
         

         function obj = setGurobiKrn(obj,sol) 
         % SETGUROBIKRN sets a kernel element using GUROBI to the class object TuKrn.
         % 
         %  Usage: clv = setGurobiKrn(clv,sol)
         %
         %  output:
         %    clv       -- TuKrn class object.
         %
         %  input:
         %     clv      -- TuKrn class object.
         %     sol      -- a kernel element (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_kr_gurobi)
                   obj.tu_kr_gurobi = [];   % Calls setter set.tu_kr_gurobi.
               end
             else
               obj.tu_kr_gurobi = sol;      % Calls setter set.tu_kr_gurobi.
             end
         end


         function obj = set.tu_kr_gurobi(obj,pt)
            if isempty(pt)
               try 
                 sol = gurobi_kernel(obj);
                 tol = obj.tu_tol;
                 pkQ = kernelQ(obj,sol,tol);
                 obj.tu_kr_gurobi = sol;
                 obj.gurobi_kr_valid = pkQ;
               catch
                 n=obj.tuplayers;
                 obj.tu_kr_gurobi = inf(1,n);
                 obj.gurobi_kr_valid = false;                 
               end    
            else
               obj.tu_kr_gurobi = pt;
               tol = obj.tu_tol;
               obj.gurobi_kr_valid = kernelQ(obj,pt,tol);
            end
         end                  
         
         
         
         function obj = setIpoptKrn(obj,sol) 
         % SETIPOPTKRN sets a kernel element using IPOPT to the class object TuKrn.
         % 
         %  Usage: clv = setIpoptKrn(clv,sol)
         %
         %  output:
         %    clv       -- TuKrn class object.
         %
         %  input:
         %     clv      -- TuKrn class object.
         %     sol      -- a kernel element (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_kr_ipopt)
                   obj.tu_kr_ipopt = [];   % Calls setter set.tu_kr_ipopt.
               end
             else
               obj.tu_kr_ipopt = sol;      % Calls setter set.tu_kr_ipopt.
             end
         end


         function obj = set.tu_kr_ipopt(obj,pt)
            if isempty(pt)
               try 
                 sol = ipopt_kernel(obj);
                 tol = obj.tu_tol;
                 pkQ = kernelQ(obj,sol,tol);
                 obj.tu_kr_ipopt = sol;
                 obj.ipopt_kr_valid = pkQ;
               catch
                 n=obj.tuplayers;
                 obj.tu_kr_ipopt = inf(1,n);
                 obj.ipopt_kr_valid = false;                 
               end    
            else
               obj.tu_kr_ipopt = pt;
               tol=obj.tu_tol;
               obj.ipopt_kr_valid = kernelQ(obj,pt,tol);
            end
         end                  

         

         function obj = setMskKrn(obj,sol) 
         % SETMSKKRN sets a kernel element using MOSEK to the class object TuKrn.
         % 
         %  Usage: clv = setMskKrn(clv,sol)
         %
         %  output:
         %    clv       -- TuKrn class object.
         %
         %  input:
         %     clv      -- TuKrn class object.
         %     sol      -- a kernel element (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_kr_msk)
                   obj.tu_kr_msk = [];   % Calls setter set.tu_kr_msk.
               end
             else
               obj.tu_kr_msk = sol;      % Calls setter set.tu_kr_msk.
             end
         end


         function obj = set.tu_kr_msk(obj,pt)
            if isempty(pt)
               try 
                 sol = msk_kernel(obj);
                 tol = obj.tu_tol;
                 pkQ = kernelQ(obj,sol,tol);
                 obj.tu_kr_msk = sol;
                 obj.msk_kr_valid = pkQ;
               catch
                 n=obj.tuplayers;
                 obj.tu_kr_msk = inf(1,n);
                 obj.msk_kr_valid = false;                 
               end    
            else
               obj.tu_kr_msk = pt;
               tol = obj.tu_tol;
               obj.msk_kr_valid = kernelQ(obj,pt,tol);
            end
         end                           
         

         function obj = setOasesKrn(obj,sol) 
         % SETOASESKRN sets a kernel element using OASES to the class object TuKrn.
         % 
         %  Usage: clv = setOasesKrn(clv,sol)
         %
         %  output:
         %    clv       -- TuKrn class object.
         %
         %  input:
         %     clv      -- TuKrn class object.
         %     sol      -- a kernel element (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_kr_oases)
                   obj.tu_kr_oases = [];   % Calls setter set.tu_kr_oases.
               end
             else
               obj.tu_kr_oases = sol;      % Calls setter set.tu_kr_oases.
             end
         end


         function obj = set.tu_kr_oases(obj,pt)
            if isempty(pt)
               try 
                 sol = oases_kernel(obj);
                 tol = obj.tu_tol;
                 pkQ = kernelQ(obj,sol,tol);
                 obj.tu_kr_oases = sol;
                 obj.oases_kr_valid = pkQ;
               catch
                 n=obj.tuplayers;
                 obj.tu_kr_oases = inf(1,n);
                 obj.oases_kr_valid = false;                 
               end    
            else
               obj.tu_kr_oases = pt;
               tol = obj.tu_tol;
               obj.oases_kr_valid = kernelQ(obj,pt,tol);
            end
         end                           
         
         
         function obj = setQpcKrn(obj,sol) 
         % SETQPCKRN sets a kernel element using QPC to the class object TuKrn.
         % 
         %  Usage: clv = setQpcKrn(clv,sol)
         %
         %  output:
         %    clv       -- TuKrn class object.
         %
         %  input:
         %     clv      -- TuKrn class object.
         %     sol      -- a kernel element (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_kr_qpc)
                   obj.tu_kr_qpc = [];   % Calls setter set.tu_kr_qpc.
               end
             else
               obj.tu_kr_qpc = sol;      % Calls setter set.tu_kr_qpc.
             end
         end


         function obj = set.tu_kr_qpc(obj,pt)
            if isempty(pt)
               try 
                 sol = qpc_kernel(obj);
                 tol = obj.tu_tol;
                 pkQ = kernelQ(obj,sol);
                 obj.tu_kr_qpc = sol;
                 obj.qpc_kr_valid = pkQ;
               catch
                 n=obj.tuplayers;
                 obj.tu_kr_qpc = inf(1,n);
                 obj.qpc_kr_valid = false;                 
               end    
            else
               obj.tu_kr_qpc = pt;
               tol = obj.tu_tol;
               obj.qpc_kr_valid = kernelQ(obj,pt,tol);
            end
         end                  
         
         
         
         
         
         function obj = setAllVendorKrn(obj)
         % SETALLVENDORKRN sets a set of kernel solutions to the class object TuKrn.
         %
         %  Usage: clv = setAllVendorKrn(clv)
         %
         %  output:
         %    clv       -- TuKrn class object.
         %
         %  input:
         %     clv      -- TuKrn class object.
         %
                %
                % CPLEX kernel
                  try 
                    solpk = cplex_kernel(obj);
                    obj.tu_kr_cplex = solpk;
                    tol = obj.tu_tol;
                    pkQ = kernelQ(obj,solpk,tol);
                    obj.cplex_kr_valid = pkQ;                  
                  catch
                    n=obj.tuplayers;
                    obj.tu_kr_cplex = inf(1,n);
                    obj.cplex_kr_valid = false;
                  end    
                %
                % CVX kernel
                  try
                    solpk = cvx_kernel(obj);
                    obj.tu_kr_cvx = solpk;
                    tol = obj.tu_tol;
                    pkQ = kernelQ(obj,solpk,tol);
                    obj.cvx_kr_valid = pkQ;
                  catch
                    n=obj.tuplayers;
                    obj.tu_kr_cvx = inf(1,n);
                    obj.cvx_kr_valid = false;                      
                  end    
                %
                % GUROBI kernel
                  try
                    solpk = gurobi_kernel(obj);
                    obj.tu_kr_gurobi = solpk;
                    tol = obj.tu_tol;
                    pkQ = kernelQ(obj,solpk,tol);
                    obj.gurobi_kr_valid = pkQ;
                  catch
                    n=obj.tuplayers;
                    obj.tu_kr_gurobi = inf(1,n);
                    obj.gurobi_kr_valid = false;                      
                  end    
                %
                % IPOPT kernel
                  try
                    solpk = ipopt_kernel(obj);
                    obj.tu_kr_ipopt = solpk;
                    tol = obj.tu_tol;
                    pkQ = kernelQ(obj,solpk,tol);
                    obj.ipopt_kr_valid = pkQ;
                  catch
                    n=obj.tuplayers;
                    obj.tu_kr_ipopt = inf(1,n);
                    obj.ipopt_kr_valid = false;                      
                  end    
                %
                % MOSEK kernel
                  try
                    solpk = msk_kernel(obj);
                    obj.tu_kr_msk = solpk;
                    tol = obj.tu_tol;
                    pkQ = kernelQ(obj,solpk,tol);
                    obj.msk_kr_valid = pkQ;
                  catch
                    n=obj.tuplayers;
                    obj.tu_kr_msk = inf(1,n);
                    obj.msk_kr_valid = false;                      
                  end                    
                %                
                % OASES kernel
                  try
                    solpk = oases_kernel(obj);
                    obj.tu_kr_oases = solpk;
                    tol = obj.tu_tol;
                    pkQ = kernelQ(obj,solpk,tol);
                    obj.oases_kr_valid = pkQ;
                  catch
                    n=obj.tuplayers;
                    obj.tu_kr_oases = inf(1,n);
                    obj.oases_kr_valid = false;                      
                  end                    
                %
                % QPC kernel
                  try
                    solpk = qpc_kernel(obj);
                    obj.tu_kr_qpc = solpk;
                    tol = obj.tu_tol;
                    pkQ = kernelQ(obj,solpk,tol);
                    obj.qpc_kr_valid = pkQ;
                  catch
                    n=obj.tuplayers;
                    obj.tu_kr_qpc = inf(1,n);
                    obj.qpc_kr_valid = false;                      
                  end                    
                %
                
        end
         
         


        function obj = copyTuSol(obj,clv)
        % copy constructor: class object clv to class object obj.
        %
        %  Usage: obj = copyTuSol(obj,clv)
        %
        %  output:
        %    obj       -- TuKrn class object.
        %
        %  input:
        %     obj      -- TuKrn class object.
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


        function obj = copy_p_TuKrn(obj,clv)
        % copy constructor: class object clv to class object obj.
        %
        %  Usage: obj = copy_p_TuKrn(obj,clv)
        %
        %  output:
        %    obj       -- TuKrn class object.
        %
        %  input:
        %     obj      -- TuKrn class object.
        %     clv      -- p_TuKrn class object.
        %
          fns = properties(clv);
             for i=1:10
              try
                 obj.(fns{i}) = clv.(fns{i});
              catch
              end
             end
        end



       end

 
end
