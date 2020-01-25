classdef TuShRep < TuSol
% TUSHREP creates the subclass object TuShRep to perform several computations for retrieving 
% and modifying game data. It stores relevant game information needed to replicate the 
% Shapley value by overloading functions.
%
% Usage: clv = TuShRep(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- TuShRep class object (subclass of TuGame).
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
%
% TuShRep properties:
%  RepShap;     -- stores the result of a Shapley value replication. For details see replicate_Shapley().
%  tu_x;        -- stores the Shapley value that has been replicated. 
%  x_sh_valid   -- returns 1 if tu_x is the Shapley value.
%  scl;         -- stores the scaling factor. 
%  tol;         -- stores the tolerance value.
%
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
% TuShRep methods:
%  TuShRep                -- creates the class object TuShRep.
%  setReplicate_Shapley   -- replicates the Shapley value of game v
%                            as the Shapley vlaue of the game space v_sp.
%  copyTuSol              -- copies all solutions from TuSol/p_TuSol to TuShRep.
%  copy_p_TuShRep         -- copies replication results from p_TuShRep to TuShRep.
%
%
%
%  Methods inherited from the superclass TuSol:
%
%  TuSol            -- creates the class object TuSol.
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
%   06/23/2013        0.4             hme
%


     properties(SetObservable = true)
       RepShap;
       tu_x;
       scl = 1;
       tol=10^6*eps;
     end

    properties(GetAccess = 'public', SetAccess = 'private')
       x_sh_valid = false;
    end


       methods
         function obj = TuShRep(w,gtype,gformat)
       % TUSHREP creates the subclass object TuShRep to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information needed to replicate the 
       % Shapley value by overloading functions.
       %
       % Usage: clv = TuShRep(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- TuShRep class object (subclass of TuGame).
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


         function obj = setReplicate_Shapley(obj,x,scl,tol)
         % SETREPLICATE_SHAPLEY sets results of the Shapley value replication to the class object TuRep.
         %
         %  Usage: clv = setReplicate_Shapley(clv,scl,tol)
         %
         %  output:
         %    clv       -- TuShRep class object.
         %
         %  input:
         %     clv      -- TuShRep class object.
         %     x        -- Shapley value of game v.    
         %     scl      -- scaling factor (default is 1)
         %     tol      -- tolerance value, default is 10^6*eps;
         %
          if nargin > 4 
             error('Too many input arguments');
          elseif nargin < 2 
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
                obj.tu_x = obj.tu_sh;
                obj.scl = 1;
                obj.tol = 10^6*eps;
                shQ = checkShapley(obj);
                obj.x_sh_valid = shQ;
             end
          elseif nargin < 3
             if isempty(obj.tu_sh)
                obj = setAllSolutions(obj);
             end
             obj.tu_x = x;
             if obj.x_sh_valid == 0
                obj.tu_sh = x;
                shQ = checkShapley(obj);
                if shQ == 1
                   obj.x_sh_valid = shQ;
                else
                   obj.tu_sh = ShapleyValue(obj);
                end
             else
                obj.tu_sh = obj.tu_x;
             end
                obj.scl = 1;
                obj.tol = 10^6*eps;
          elseif nargin < 4
             obj.scl = scl;
             if isempty(obj.tu_sh)
                obj = setAllSolutions(obj);
             end
             obj.tu_x = x;
             if obj.x_sh_valid == 0
                obj.tu_sh = ShapleyValue(obj);
                shQ = checkShapley(obj);
                if shQ == 1
                   obj.x_sh_valid = shQ;
                else 
                   obj.tu_sh = ShapleyValue(obj);
                end
             else
                obj.tu_x = obj.tu_sh;
             end
          else
             if isempty(obj.tu_sh)
                obj = setAllSolutions(obj);
             end              
             obj.tu_x = x;
             obj.scl = scl;
             obj.tol = tol;
             if obj.x_sh_valid == 0
                obj.tu_sh = ShapleyValue(obj);
                shQ = checkShapley(obj);
                if shQ == 1
                   obj.x_sh_valid = shQ;
                else 
                   obj.tu_sh = ShapleyValue(obj);
                end
             else
                obj.tu_x = obj.tu_sh;
             end
          end       
             obj.RepShap=replicate_Shapley(obj,obj.scl,obj.tol);
         end

         function obj = copyTuSol(obj,clv)
         % copy constructor: class object clv to class object obj.
         %
         %  Usage: obj = copyTuSol(obj,clv)
         %
         %  output:
         %    obj       -- TuShRep class object.
         %
         %  input:
         %     obj      -- TuShRep class object.
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


         function obj = copy_p_TuShRep(obj,clv)
         % copy constructor: class object clv to class object obj.
         %
         %  Usage: obj = copy_p_TuShRep(obj,clv)
         %
         %  output:
         %    obj       -- TuShRep class object.
         %
         %  input:
         %     obj      -- TuShRep class object.
         %     clv      -- p_TuShRep class object.
         %
           fns = properties(clv);
              for i=1:16
               try
                  obj.(fns{i}) = clv.(fns{i});
               catch
               end
              end
         end


      
       end

       methods(Access = 'private')

         function shQ = checkShapley(obj)
           if isempty(obj.tu_sh)
              return;
           elseif isempty(obj.tu_x)
               obj.tu_x=ShapleyValue(obj);
           end              
           shQ=all(abs(ShapleyValue(obj)-obj.tu_x)<obj.tol);
         end 

       end


end
