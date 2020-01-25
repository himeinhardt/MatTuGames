classdef TuCons < TuSol
% TUCONS creates the subclass object TuCons to perform several computations for retrieving 
% and modifying game data. It stores relevant game information needed to apply various
% consistency investigations on a proposed solution by overloading functions. 
%
% Usage: clv = TuCons(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- TuCons class object (subclass of TuGame).
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
% TuCons properties:
%
%  tu_x             -- stores the payoff vector to be checked on consistency.
%  tu_CRGP_PRN      -- stores the information whether tu_x satisfies CRGP w.r.t. DM-reduced game in accordance with the pre-nucleolus.
%  tu_CRGPC_PRN     -- stores the information w.r.t. the DM-reduced games restricted to tu_x. 
%  tu_CRGP_PRK      -- stores the information whether tu_x satisfies CRGP w.r.t. DM-reduced game in accordance with the pre-kernel.
%  tu_CRGPC_PRK     -- stores the information w.r.t. the DM-reduced games restricted to tu_x.
%  tu_CRGP_SHAP     -- stores the information whether tu_x satisfies CRGP w.r.t. HMS-reduced game in accordance with the Shapley value.
%  tu_CRGPC_SHAP    -- stores the information w.r.t. the HMS-reduced games restricted to tu_x.
%  tu_RCP_PRN       -- stores the information whether tu_x satisfies RCP w.r.t. DM-reduced game in accordance with the pre-nucleolus.
%  tu_RCPC_PRN      -- stores the information w.r.t. the DM-reduced games restricted to tu_x.
%  tu_RCP_PRK       -- stores the information whether tu_x satisfies RCP w.r.t. DM-reduced game in accordance with the pre-kernel.
%  tu_RCPC_PRK      -- stores the information w.r.t. the DM-reduced games restricted to tu_x.
%  tu_RCP_SHAP      -- stores the information whether tu_x satisfies RCP w.r.t. HMS-reduced game in accordance with the Shapley value.
%  tu_RCPC_SHAP     -- stores the information w.r.t. the HMS-reduced games restricted to tu_x.
%  tu_RCP_HMS_PN    -- stores the information whether tu_x satisfies RCP w.r.t. HMS-reduced game in accordance with the pre-nucleolus.
%  tu_RCPC_HMS_PN   -- stores the information w.r.t. the HMS-reduced games restricted to tu_x.
%  tu_RCP_HMS_PK    -- stores the information whether tu_x satisfies RCP w.r.t. HMS-reduced game in accordance with the pre-kernel.
%  tu_RCPC_HMS_PK   -- stores the information w.r.t. the HMS-reduced games restricted to tu_x.
%  tu_RGP_PRN       -- stores the information whether tu_x satisfies RGP w.r.t. DM-reduced game in accordance with the pre-nucleolus.
%  tu_RGPC_PRN      -- stores the information w.r.t. the DM-reduced games restricted to tu_x.
%  tu_RGP_PRK       -- stores the information whether tu_x satisfies RGP w.r.t. DM-reduced game in accordance with the pre-kernel.
%  tu_RGPC_PRK      -- stores the information w.r.t. the DM-reduced games restricted to tu_x.
%  tu_RGP_SHAP      -- stores the information whether tu_x satisfies RGP w.r.t. HMS-reduced game in accordance with the Shapley value 
%  tu_RGPC_SHAP     -- stores the information w.r.t. the HMS-reduced games restricted to tu_x.
%  tu_RGP_HMS_PN    -- stores the information whether tu_x satisfies RGP w.r.t. HMS-reduced game in accordance with the pre-nucleolus.
%  tu_RGPC_HMS_PN   -- stores the information w.r.t. the HMS-reduced games restricted to tu_x.
%  tu_RGP_HMS_PK    -- stores the information whether tu_x satisfies RGP w.r.t. HMS-reduced game in accordance with the pre-kernel.
%  tu_RGPC_HMS_PK   -- stores the information w.r.t. the HMS-reduced games restricted to tu_x.
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
% TuCons methods:
%  TuCons                       -- creates the class object TuCons.
%  setConverse_RGP              -- sets results of CRGP to the class object TuCons.
%  setReconfirmation_property   -- sets results of RCP to the class object TuCons.
%  setReduced_game_property     -- sets resutls of RGP to the class object TuCons.
%  copyTuSol                    -- copies all solutions from TuSol/p_TuSol to TuCons.
%  copy_p_TuCons                -- copies consistency properties from p_TuCons to TuCons.
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
%   11/11/2012        0.3             hme
%





      properties(SetObservable = true)
      tu_x
      tu_CRGP_PRN
      tu_CRGPC_PRN
      tu_CRGP_PRK
      tu_CRGPC_PRK
      tu_CRGP_SHAP
      tu_CRGPC_SHAP
      tu_RCP_PRN
      tu_RCPC_PRN
      tu_RCP_PRK
      tu_RCPC_PRK
      tu_RCP_SHAP
      tu_RCPC_SHAP
      tu_RCP_HMS_PN
      tu_RCPC_HMS_PN
      tu_RCP_HMS_PK
      tu_RCPC_HMS_PK
      tu_RGP_PRN
      tu_RGPC_PRN
      tu_RGP_PRK
      tu_RGPC_PRK
      tu_RGP_SHAP
      tu_RGPC_SHAP
      tu_RGP_HMS_PN
      tu_RGPC_HMS_PN
      tu_RGP_HMS_PK
      tu_RGPC_HMS_PK
      end
  

       methods
         function obj = TuCons(w,gtype,gformat)
       % TUCONS creates the subclass object TuCons to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information needed to apply various
       % consistency investigations on a proposed solution by overloading functions.
       %
       % Usage: clv = TuCons(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- TuCons class object (subclass of TuGame).
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

         function obj = setConverse_RGP(obj,x)
         % SETCONVERSE_RGP sets results of CRGP to the class object TuCons.
         %
         %  Usage: clv = setConverse_RGP(clv,x) 
         %
         %  output:
         %    clv       -- TuCons class object.
         %
         %  input:
         %     clv      -- TuCons class object.
         %     x        -- payoff vector (efficient).
         %        
          if nargin > 2
             error('Too many input arguments');
          elseif nargin < 2
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             end
             [obj.tu_CRGP_PRN, obj.tu_CRGPC_PRN]=Converse_RGP_Q(obj,obj.tu_x,'PRN');
             [obj.tu_CRGP_PRK, obj.tu_CRGPC_PRK]=Converse_RGP_Q(obj,obj.tu_x,'PRK');
             [obj.tu_CRGP_SHAP, obj.tu_CRGPC_SHAP]=Converse_RGP_Q(obj,obj.tu_x,'SHAP');          
          else
            if all(size(x)==[1, obj.tuplayers ])
             obj.tu_x = x;
            else
              error('Payoff vector has not the correct size! It has to be (1,n).')
            end
             [obj.tu_CRGP_PRN, obj.tu_CRGPC_PRN]=Converse_RGP_Q(obj,obj.tu_x,'PRN');
             [obj.tu_CRGP_PRK, obj.tu_CRGPC_PRK]=Converse_RGP_Q(obj,obj.tu_x,'PRK');
             [obj.tu_CRGP_SHAP, obj.tu_CRGPC_SHAP]=Converse_RGP_Q(obj,obj.tu_x,'SHAP');
          end
         end


         function obj = setReconfirmation_property(obj,x)
         % SETRECONFIRMATION_PROPERTY sets results of RCP to the class object TuCons.
         %
         %  Usage: clv = setReconfirmation_property(clv,x)
         %
         %  output:
         %    clv       -- TuCons class object.
         %
         %  input:
         %     clv      -- TuCons class object.
         %     x        -- payoff vector (efficient).
         %
           if nargin > 2
              error('Too many input arguments');
           elseif nargin < 2
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             end
                [obj.tu_RCP_PRN, obj.tu_RCPC_PRN]=Reconfirmation_propertyQ(obj,obj.tu_x,'PRN');
                [obj.tu_RCP_PRK, obj.tu_RCPC_PRK]=Reconfirmation_propertyQ(obj,obj.tu_x,'PRK');
                [obj.tu_RCP_SHAP, obj.tu_RCPC_SHAP]=Reconfirmation_propertyQ(obj,obj.tu_x,'SHAP');
                [obj.tu_RCP_HMS_PN, obj.tu_RCPC_HMS_PN]=Reconfirmation_propertyQ(obj,obj.tu_x,'HMS_PN');
                [obj.tu_RCP_HMS_PK, obj.tu_RCPC_HMS_PK]=Reconfirmation_propertyQ(obj,obj.tu_x,'HMS_PK');
           else
            if all(size(x)==[1, obj.tuplayers ])
             obj.tu_x = x;
            else
              error('Payoff vector has not the correct size! It has to be (1,n).')
            end
                [obj.tu_RCP_PRN, obj.tu_RCPC_PRN]=Reconfirmation_propertyQ(obj,obj.tu_x,'PRN');
                [obj.tu_RCP_PRK, obj.tu_RCPC_PRK]=Reconfirmation_propertyQ(obj,obj.tu_x,'PRK');
                [obj.tu_RCP_SHAP, obj.tu_RCPC_SHAP]=Reconfirmation_propertyQ(obj,obj.tu_x,'SHAP');
                [obj.tu_RCP_HMS_PN, obj.tu_RCPC_HMS_PN]=Reconfirmation_propertyQ(obj,obj.tu_x,'HMS_PN');
                [obj.tu_RCP_HMS_PK, obj.tu_RCPC_HMS_PK]=Reconfirmation_propertyQ(obj,obj.tu_x,'HMS_PK');
           end
         end

         function obj = setReduced_game_property(obj,x)
         % SETREDUCED_GAME_PROPERTY sets results of RCP to the class object TuCons.
         %
         %  Usage: clv = setReduced_game_property(clv,x)
         %
         %  output:
         %    clv       -- TuCons class object.
         %
         %  input:
         %     clv      -- TuCons class object.
         %     x        -- payoff vector (efficient).
         %
           if nargin > 2 
             error('Too many input arguments');
           elseif nargin < 2
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             end
                [obj.tu_RGP_PRN, obj.tu_RGPC_PRN]=Reduced_game_propertyQ(obj,obj.tu_x,'PRN');
                [obj.tu_RGP_PRK, obj.tu_RGPC_PRK]=Reduced_game_propertyQ(obj,obj.tu_x,'PRK');
                [obj.tu_RGP_SHAP, obj.tu_RGPC_SHAP]=Reduced_game_propertyQ(obj,obj.tu_x,'SHAP');
                [obj.tu_RGP_HMS_PN, obj.tu_RGPC_HMS_PN]=Reduced_game_propertyQ(obj,obj.tu_x,'HMS_PN');
                [obj.tu_RGP_HMS_PK, obj.tu_RGPC_HMS_PK]=Reduced_game_propertyQ(obj,obj.tu_x,'HMS_PK');
           else
            if all(size(x)==[1, obj.tuplayers ])
             obj.tu_x = x;
            else
              error('Payoff vector has not the correct size! It has to be (1,n).')
            end
                [obj.tu_RGP_PRN, obj.tu_RGPC_PRN]=Reduced_game_propertyQ(obj,obj.tu_x,'PRN');
                [obj.tu_RGP_PRK, obj.tu_RGPC_PRK]=Reduced_game_propertyQ(obj,obj.tu_x,'PRK');
                [obj.tu_RGP_SHAP, obj.tu_RGPC_SHAP]=Reduced_game_propertyQ(obj,obj.tu_x,'SHAP');
                [obj.tu_RGP_HMS_PN, obj.tu_RGPC_HMS_PN]=Reduced_game_propertyQ(obj,obj.tu_x,'HMS_PN');
                [obj.tu_RGP_HMS_PK, obj.tu_RGPC_HMS_PK]=Reduced_game_propertyQ(obj,obj.tu_x,'HMS_PK');
           end
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


        function obj = copy_p_TuCons(obj,clv)
        % copy constructor: class object clv to class object obj.
        %
        %  Usage: obj = copy_p_TuCons(obj,clv)
        %
        %  output:
        %    obj       -- TuCons class object.
        %
        %  input:
        %     obj      -- TuCons class object.
        %     clv      -- p_TuCons class object.
        %
          fns = properties(clv);
             for i=1:26
              try
                 obj.(fns{i}) = clv.(fns{i});
              catch
              end
             end
        end



       end

 
end
