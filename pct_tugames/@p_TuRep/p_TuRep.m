classdef p_TuRep < p_TuSol
% P_TUREP creates the subclass object TuRep to perform several computations for retrieving 
% and modifying game data. It stores relevant game information needed to replicate a 
% pre-kernel element by overloading functions using Matlab's PCT.
%
% Usage: clv = p_TuRep(v,'gtype','gformat',K)
%
% Define variables:
% output:
% clv           -- p_TuRep class object (subclass of TuGame).
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
% K             -- A positive integer between 2 and n.
%
%
% p_TuRep properties:
%  RepSol;      -- stores the result of a pre-kernel element replication. For details see replicate_prk().
%  tu_x;        -- stores the pre-kernel vector that has been replicated. 
%  x_prk_valid  -- returns 1 if tu_x is a pre-kernel element.
%  scl;         -- stores the scaling factor. 
%  smc;         -- stores the cardinality of most effective coalitions, that is, smallest/largest (1/0).
%
%
%  Properties inherited from the superclass p_TuSol:
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
%  tuvi         -- stores the values of singelton coalitions.
%  tustpt       -- stores a starting point for doing computation. Has lower priority than
%                  providing a second input argument as an alternative starting point.
%                  Thus, if a starting point is provided as second input argument the
%                  information stored in tustpt will not be used.
%
% p_TuRep methods:
%  p_TuRep                -- creates the class object p_TuRep.
%  p_setReplicate_Prk     -- replicates a pre-kernel solution x as a pre-kernel of
%                            the game space v_sp.
%  p_copyTuSol            -- copies all solutions from TuSol/p_TuSol to p_TuRep.
%  copyTuRep              -- copies replication results from TuRep to p_TuRep.
%
%
%
%  Methods inherited from the superclass p_TuSol:
%
%  p_TuSol              -- creates the class object p_TuSol.
%  p_setAllSolutions    -- sets all solutions listed below to the class object p_TuSol.
%  p_setPreKernel       -- sets a pre-kernel element to the class object p_TuSol.
%  p_setPreNuc          -- sets the pre-nucleolus to the class object p_TuSol.
%  p_setShapley         -- sets the Shapley value to the class object p_TuSol.
%  p_setTauValue        -- sets the Tau value to the class object p_TuSol.
%  p_setBanzhaf         -- sets the Banzhaf value to the class object p_TuSol.
%  p_setAntiPreKernel   -- sets an anti-pre-kernel element to the class object p_TuSol.
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
         function obj = p_TuRep(w,gtype,gformat)
       % P_TUREP creates the subclass object p_TuRep to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information needed to replicate a 
       % pre-kernel element by overloading functions using Matlab's PCT. 
       %
       % Usage: clv = p_TuRep(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- p_TuRep class object (subclass of TuGame).
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
           obj = obj@p_TuSol(w,gtype,gformat);
         end


         function obj = p_setReplicate_Prk(obj,x,scl,smc)
         % P_SETREPLICATE_PRK sets results of pre-kernel replication to the class object p_TuRep.
         %
         %  Usage: clv = p_setReplicate_Prk(clv,x,scl,smc)
         %
         %  output:
         %    clv       -- p_TuRep class object.
         %
         %  input:
         %     clv      -- p_TuRep class object.
         %     x        -- pre-kernel vector.
         %     scl      -- scaling factor (default is 1)
         %     smc      -- selecting from effc the smallest/largest (1/0)
         %                 cardinality. default is 1. 
         %
          if nargin > 4 
             error('Too many input arguments');
          elseif nargin < 2
             if isempty(obj.tu_prk)
                obj = p_setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
                obj.x_prk_valid = obj.prk_valid;
             end
          elseif nargin < 3
             if isempty(obj.tu_prk)
                obj = p_setAllSolutions(obj);
             end
             obj.tu_prk = x;
             if obj.prk_valid == 0
                obj.tu_prk = p_PreKernel(obj);
                pkQ = p_checkPreKernel(obj);
                obj.x_prk_valid = pkQ;
                obj.tu_x = obj.tu_prk;
             else
                obj.tu_x = obj.tu_prk;
                obj.x_prk_valid = obj.prk_valid;
             end
          elseif nargin < 4
             obj.scl = scl; 
             if isempty(obj.tu_prk)
                obj = p_setAllSolutions(obj);
             end
             obj.tu_prk = x;
             if obj.prk_valid == 0
                obj.tu_prk = p_PreKernel(obj);
                pkQ = p_checkPreKernel(obj);
                obj.x_prk_valid = pkQ;
                obj.tu_x = obj.tu_prk;
             else
                obj.tu_x = obj.tu_prk;
                obj.x_prk_valid = obj.prk_valid;
             end
          else
             if isempty(obj.tu_prk)
                obj = p_setAllSolutions(obj);
             end
             obj.scl = scl;
             obj.tu_prk = x;
             if obj.prk_valid == 0
                obj.tu_prk = p_PreKernel(obj);
                pkQ = p_checkPreKernel(obj);
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
             obj.RepSol=p_replicate_prk(obj,obj.tu_prk,obj.scl,obj.smc);
         end

         function obj = p_copyTuSol(obj,clv)
         % copy constructor: class object clv to class object obj.
         %
         %  Usage: obj = p_copyTuSol(obj,clv)
         %
         %  output:
         %    obj       -- p_TuRep class object.
         %
         %  input:
         %     obj      -- p_TuRep class object.
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

         function obj = copyTuRep(obj,clv)
         % copy constructor: class object clv to class object obj.
         %
         %  Usage: obj = p_copyTuSol(obj,clv)
         %
         %  output:
         %    obj       -- p_TuRep class object.
         %
         %  input:
         %     obj      -- p_TuRep class object.
         %     clv      -- TuRep class object.
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

         function pkQ = p_checkPreKernel(obj)
           if isempty(obj.tu_prk)
              return;
           end
           pkQ=p_PrekernelQ(obj,obj.tu_prk);
         end 

       end


end
