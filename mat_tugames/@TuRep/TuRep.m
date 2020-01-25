classdef TuRep < TuSol
% TUREP creates the subclass object TuRep to perform several computations for retrieving 
% and modifying game data. It stores relevant game information needed to replicate a 
% pre-kernel element by overloading functions.
%
% Usage: clv = TuRep(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- TuRep class object (subclass of TuGame).
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
% TuRep properties:
%  RepSol;      -- stores the result of a pre-kernel element replication. For details see replicate_prk().
%  tu_x;        -- stores the pre-kernel vector that has been replicated. 
%  x_prk_valid  -- returns 1 if tu_x is a pre-kernel element.
%  scl;         -- stores the scaling factor. 
%  smc;         -- stores the cardinality of most effective coalitions, that is, smallest/largest (1/0).
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
% TuRep methods:
%  TuRep                  -- creates the class object TuRep.
%  setReplicate_Prk       -- replicates a pre-kernel solution x as a pre-kernel of
%                            the game space v_sp.
%  copyTuSol              -- copies all solutions from TuSol/p_TuSol to TuRep.
%  copy_p_TuRep           -- copies replication results from p_TuRep to TuRep.
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
%   11/17/2012        0.3             hme
%   06/23/2013        0.4             hme
%


     properties(SetObservable = true)
       RepSol;
       tu_x;
       scl = 1;
       smc = 1;
     end

    properties(GetAccess = 'public', SetAccess = 'private')
      x_prk_valid = false;
    end


       methods
         function obj = TuRep(w,gtype,gformat)
       % TUREP creates the subclass object TuRep to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information needed to replicate a 
       % pre-kernel element by overloading functions.
       %
       % Usage: clv = TuRep(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- TuRep class object (subclass of TuGame).
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


         function obj = setReplicate_Prk(obj,x,scl,smc)
         % SETREPLICATE_PRK sets results of pre-kernel replication to the class object TuRep.
         %
         %  Usage: clv = setReplicate_Prk(clv,x,scl,smc)
         %
         %  output:
         %    clv       -- TuRep class object.
         %
         %  input:
         %     clv      -- TuRep class object.
         %     x        -- pre-kernel vector.
         %     scl      -- scaling factor (default is 1)
         %     smc      -- selecting from effc the smallest/largest (1/0)
         %                 cardinality. default is 1. 
         %
          if nargin > 4 
             error('Too many input arguments');
          elseif nargin < 2 
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
                obj.x_prk_valid = obj.prk_valid;
             end
          elseif nargin < 3
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
             end
             obj.tu_prk = x;
             if obj.prk_valid == 0
                obj.tu_prk = PreKernel(obj);
                pkQ = checkPreKernel(obj);
                obj.x_prk_valid = pkQ;
                obj.tu_x = obj.tu_prk;
             else
                obj.tu_x = obj.tu_prk;
                obj.x_prk_valid = obj.prk_valid;
             end
          elseif nargin < 4
             obj.scl = scl;
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
             end
             obj.tu_prk = x;
             if obj.prk_valid == 0
                obj.tu_prk = PreKernel(obj);
                pkQ = checkPreKernel(obj);
                obj.x_prk_valid = pkQ;
                obj.tu_x = obj.tu_prk;
             else
                obj.tu_x = obj.tu_prk;
                obj.x_prk_valid = obj.prk_valid;
             end
          else
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
             end
             obj.scl = scl;
             obj.tu_prk = x;
             if obj.prk_valid == 0
                obj.tu_prk = PreKernel(obj);
                pkQ = checkPreKernel(obj);
                obj.x_prk_valid = pkQ;
                obj.tu_x = obj.tu_prk;
             else
                obj.tu_x = obj.tu_prk;
                obj.x_prk_valid = obj.prk_valid;
             end
             if smc > 1
                obj.smc = 1;
             elseif smc < 0;
                obj.smc = 1;
             else
                obj.smc=round(smc);
             end
          end       
             obj.RepSol=replicate_prk(obj,obj.tu_prk,obj.scl,obj.smc);
         end

         function obj = copyTuSol(obj,clv)
         % copy constructor: class object clv to class object obj.
         %
         %  Usage: obj = copyTuSol(obj,clv)
         %
         %  output:
         %    obj       -- TuRep class object.
         %
         %  input:
         %     obj      -- TuRep class object.
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


         function obj = copy_p_TuRep(obj,clv)
         % copy constructor: class object clv to class object obj.
         %
         %  Usage: obj = copy_p_TuRep(obj,clv)
         %
         %  output:
         %    obj       -- TuRep class object.
         %
         %  input:
         %     obj      -- TuRep class object.
         %     clv      -- p_TuRep class object.
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

         function pkQ = checkPreKernel(obj)
           if isempty(obj.tu_prk)
              return;
           end
           pkQ=PrekernelQ(obj,obj.tu_prk);
         end 

       end


end
