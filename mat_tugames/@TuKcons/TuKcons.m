classdef TuKcons < TuSol 
% TUKCONS creates the subclass object TuKcons to perform several computations for retrieving 
% and modifying game data. It stores relevant game information needed to appy various
% generalized consistency investigations on a proposed solution by overloading functions.
%
% Usage: clv = TuKcons(v,'gtype','gformat',K)
%
% Define variables:
% output:
% clv           -- TuKcons class object (subclass of TuGame).
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
% TuKcons properties:
%
%  tu_x              -- stores the payoff vector to be checked on consistency.
%  tu_K              -- A positive integer between 2 and n.
%  tu_kCRGP_PRN      -- stores the information whether tu_x satisfies K-CRGP w.r.t. DM-reduced game in accordance with the pre-nucleolus.
%  tu_kCRGP_PRK      -- stores the information whether tu_x satisfies K-CRGP w.r.t. DM-reduced game in accordance with the pre-kernel.
%  tu_kCRGP_SHAP     -- stores the information whether tu_x satisfies K-CRGP w.r.t. HMS-reduced game in accordance with the Shapley value.
%  tu_kRCP_PRN       -- stores the information whether tu_x satisfies K-RCP w.r.t. DM-reduced game in accordance with the pre-nucleolus.
%  tu_kRCP_PRK       -- stores the information whether tu_x satisfies K-RCP w.r.t. DM-reduced game in accordance with the pre-kernel.
%  tu_kRCP_SHAP      -- stores the information whether tu_x satisfies K-RCP w.r.t. HMS-reduced game in accordance with the Shapley value.
%  tu_kRCP_HMS_PN    -- stores the information whether tu_x satisfies K-RCP w.r.t. HMS-reduced game in accordance with the pre-nucleolus.
%  tu_kRCP_HMS_PK    -- stores the information whether tu_x satisfies K-RCP w.r.t. HMS-reduced game in accordance with the pre-kernel.
%  tu_kSCRGP_PRN     -- stores the information whether tu_x satisfies K-SCRGP w.r.t. DM-reduced game in accordance with the pre-nucleolus.
%  tu_kSCRGP_PRK     -- stores the information whether tu_x satisfies K-SCRGP w.r.t. DM-reduced game in accordance with the pre-kernel.
%  tu_kSCRGP_SHAP    -- stores the information whether tu_x satisfies K-SCRGP w.r.t. HMS-reduced game in accordance with the Shapley value 
%  tu_kRGP_PRN       -- stores the information whether tu_x satisfies KRGP w.r.t. DM-reduced game in accordance with the pre-nucleolus.
%  tu_kRGP_PRK       -- stores the information whether tu_x satisfies K-RGP w.r.t. DM-reduced game in accordance with the pre-kernel.
%  tu_kRGP_SHAP      -- stores the information whether tu_x satisfies K-RGP w.r.t. HMS-reduced game in accordance with the Shapley value.
%  tu_kRGP_HMS_PN    -- stores the information whether tu_x satisfies K-RGP w.r.t. HMS-reduced game in accordance with the pre-nucleolus.
%  tu_kRGP_HMS_PN    -- stores the information whether tu_x satisfies K-RGP w.r.t. HMS-reduced game in accordance with the pre-kernel.
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
% TuKcons methods:
%  TuKcons                       -- creates the class object TuKcons.
%  setKConverse_RGP              -- sets results of k-CRGP to the class object TuKcons.
%  setKReconfirmation_property   -- sets results of k-RCP to the class object TukCons.
%  setKStrConverse_RGP           -- sets results of k-SCRGP to the class object TuKcons.
%  setKReducedGameProperty       -- sets results of k-RGP to the class object TuKcons.  
%                                   k is a positive integer between 2 and n. 
%  copyTuSol                     -- copies all solutions from TuSol/p_TuSol to TuKcons.
%  copy_p_TuKcons                -- copies k-consistency properties from p_TuKcons to TuKcons.
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
      tu_K = 2;
      tu_kCRGP_PRN
      tu_kCRGP_PRK
      tu_kCRGP_SHAP
      tu_kRCP_PRN
      tu_kRCP_PRK
      tu_kRCP_SHAP
      tu_kRCP_HMS_PN
      tu_kRCP_HMS_PK
      tu_kSCRGP_PRN
      tu_kSCRGP_PRK
      tu_kSCRGP_SHAP
      tu_kRGP_PRN
      tu_kRGP_PRK
      tu_kRGP_SHAP
      tu_kRGP_HMS_PN
      tu_kRGP_HMS_PK
      end
  

       methods
         function obj = TuKcons(w,gtype,gformat,K)
       % TUCKONS creates the subclass object TuKcons to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information needed to appy various
       % generalized consistency investigations on a proposed solution by overloading functions.
       %
       % Usage: clv = TuKcons(v,'gtype','gformat',K)
       %
       % Define variables:
       % output:
       % clv           -- TuKcons class object (subclass of TuGame).
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
       % K             -- A positive integer between 2 and n.
       %
           if nargin > 4
              error('Too many input arguments');
           elseif nargin < 1
              error('Game information must be given as a 2^n-1 vector!');
           elseif nargin < 2 
              gtype = '';
              gformat = 'mattug';
              K =2;
           elseif nargin < 3
              gformat = 'mattug';
              K = 2;
           elseif nargin < 4
              K = 2;
           else
           end
           obj = obj@TuSol(w,gtype,gformat);
           obj.tu_K = K; 
         end

         function obj = setKConverse_RGP(obj,x,K)
         % SETCONVERSE_RGP sets results of CRGP to the class object TuCons.
         %
         %  Usage: clv = setKConverse_RGP(clv,x,K) 
         %
         %  output:
         %    clv       -- TuKcons class object.
         %
         %  input:
         %     clv      -- TuKcons class object.
         %     x        -- payoff vector (efficient).
         %     K        -- a positve integer between 2 and n.
         %                 default tu_K=2.
         %        
          if nargin > 3 
             error('Too many input arguments');
          elseif nargin < 2
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             end
             obj.tu_kCRGP_PRN=k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRN');
             obj.tu_kCRGP_PRK=k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRK');
             obj.tu_kCRGP_SHAP=k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'SHAP');
          elseif nargin < 3
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             else
                if all(size(x)==[1, obj.tuplayers ])
                  obj.tu_x = x;
                else
                  error('Payoff vector has not the correct size! It has to be (1,n).')
                end
             end
             obj.tu_kCRGP_PRN=k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRN');
             obj.tu_kCRGP_PRK=k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRK');
             obj.tu_kCRGP_SHAP=k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'SHAP');
          else
            if(K < 2 || K > obj.tuplayers) 
              error('K must be an integer between 2 and n!')
            end
            obj.tu_K = K;
            if all(size(x)==[1, obj.tuplayers ])
             obj.tu_x = x;
            else
              error('Payoff vector has not the correct size! It has to be (1,n).')
            end
             obj.tu_kCRGP_PRN=k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRN');
             obj.tu_kCRGP_PRK=k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRK');
             obj.tu_kCRGP_SHAP=k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'SHAP');
          end
         end


         function obj = setKReconfirmation_property(obj,x,K)
         % SETRECONFIRMATION_PROPERTY sets results of RCP to the class object TuCons.
         %
         %  Usage: clv = setKReconfirmation_property(clv,x,K)
         %
         %  output:
         %    clv       -- TuKcons class object.
         %
         %  input:
         %     clv      -- TuKcons class object.
         %     x        -- payoff vector (efficient).
         %     K        -- a positve integer between 2 and n.
         %                 default K=2.
         %
           if nargin > 3 
              error('Too many input arguments');
           elseif nargin < 2
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             end
                obj.tu_kRCP_PRN=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'PRN');
                obj.tu_kRCP_PRK=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'PRK');
                obj.tu_kRCP_SHAP=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'SHAP');
                obj.tu_kRCP_HMS_PN=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PN');
                obj.tu_kRCP_HMS_PK=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PK');
           elseif nargin < 3
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             else
                if all(size(x)==[1, obj.tuplayers ])
                  obj.tu_x = x;
                else
                  error('Payoff vector has not the correct size! It has to be (1,n).')
                end
             end
                obj.tu_kRCP_PRN=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'PRN');
                obj.tu_kRCP_PRK=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'PRK');
                obj.tu_kRCP_SHAP=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'SHAP');
                obj.tu_kRCP_HMS_PN=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PN');
                obj.tu_kRCP_HMS_PK=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PK');
           else
            if(K < 2 || K > obj.tuplayers)
              error('K must be an integer between 2 and n!')
            end
            obj.tu_K = K;
            if all(size(x)==[1, obj.tuplayers ])
             obj.tu_x = x;
            else
              error('Payoff vector has not the correct size! It has to be (1,n).')
            end
                obj.tu_kRCP_PRN=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'PRN');
                obj.tu_kRCP_PRK=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'PRK');
                obj.tu_kRCP_SHAP=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'SHAP');
                obj.tu_kRCP_HMS_PN=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PN');
                obj.tu_kRCP_HMS_PK=k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PK');
           end
         end

         function obj = setKStrConverse_RGP(obj,x,K)
         % SETREDUCED_GAME_PROPERTY sets results of RCP to the class object TuCons.
         %
         %  Usage: clv = setKStrConverse_RGP(clv,x,K)
         %
         %  output:
         %    clv       -- TuKcons class object.
         %
         %  input:
         %     clv      -- TuKcons class object.
         %     x        -- payoff vector (efficient).
         %     K        -- a positve integer between 2 and n.
         %                 default K=2. 
         %
           if nargin > 3 
             error('Too many input arguments');
           elseif nargin < 2
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             end
             obj.tu_kSCRGP_PRN=k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRN');
             obj.tu_kSCRGP_PRK=k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRK');
             obj.tu_kSCRGP_SHAP=k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'SHAP');
           elseif nargin < 3
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             else
                if all(size(x)==[1, obj.tuplayers ])
                  obj.tu_x = x;
                else
                  error('Payoff vector has not the correct size! It has to be (1,n).')
                end
             end
             obj.tu_kSCRGP_PRN=k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRN');
             obj.tu_kSCRGP_PRK=k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRK');
             obj.tu_kSCRGP_SHAP=k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'SHAP');
           else
            if(K < 2 || K > obj.tuplayers)
              error('K must be an integer between 2 and n!')
            end
            obj.tu_K = K;
            if all(size(x)==[1, obj.tuplayers ])
             obj.tu_x = x;
            else
              error('Payoff vector has not the correct size! It has to be (1,n).')
            end
             obj.tu_kSCRGP_PRN=k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRN');
             obj.tu_kSCRGP_PRK=k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRK');
             obj.tu_kSCRGP_SHAP=k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'SHAP');
           end
         end


         function obj = setKReducedGameProperty(obj,x,K)
         % SETKREDUCEDGAMEPROPERTY sets results of RGP to the class object TuCons.
         %
         %  Usage: clv = setKReducedGameProperty(clv,x,K)
         %
         %  output:
         %    clv       -- TuKcons class object.
         %
         %  input:
         %     clv      -- TuKcons class object.
         %     x        -- payoff vector (efficient).
         %     K        -- a positve integer between 2 and n.
         %                 default K=2.
         %
           if nargin > 3
              error('Too many input arguments');
           elseif nargin < 2
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             end
                obj.tu_kRGP_PRN=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'PRN');
                obj.tu_kRGP_PRK=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'PRK');
                obj.tu_kRGP_SHAP=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'SHAP');
                obj.tu_kRGP_HMS_PN=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PN');
                obj.tu_kRGP_HMS_PK=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PK');
           elseif nargin < 3
             if isempty(obj.tu_prk)
                obj = setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             else
                if all(size(x)==[1, obj.tuplayers ])
                  obj.tu_x = x;
                else
                  error('Payoff vector has not the correct size! It has to be (1,n).')
                end
             end
                obj.tu_kRGP_PRN=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'PRN');
                obj.tu_kRGP_PRK=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'PRK');
                obj.tu_kRGP_SHAP=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'SHAP');
                obj.tu_kRGP_HMS_PN=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PN');
                obj.tu_kRGP_HMS_PK=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PK');
           else
            if(K < 2 || K > obj.tuplayers)
              error('K must be an integer between 2 and n!')
            end
            obj.tu_K = K;
            if all(size(x)==[1, obj.tuplayers ])
             obj.tu_x = x;
            else
              error('Payoff vector has not the correct size! It has to be (1,n).')
            end
                obj.tu_kRGP_PRN=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'PRN');
                obj.tu_kRGP_PRK=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'PRK');
                obj.tu_kRGP_SHAP=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'SHAP');
                obj.tu_kRGP_HMS_PN=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PN');
                obj.tu_kRGP_HMS_PK=k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PK');
           end
         end


        function obj = copyTuSol(obj,clv)
        % copy constructor: class object clv to class object obj.
        %
        %  Usage: obj = copyTuSol(obj,clv)
        %
        %  output:
        %    obj       -- TuKcons class object.
        %
        %  input:
        %     obj      -- TuKcons class object.
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


         function obj = copy_p_TuKcons(obj,clv)
         % copy constructor: class object clv to class object obj.
         %
         %  Usage: obj = copy_p_TuKcons(obj,clv)
         %
         %  output:
         %    obj       -- TuKcons class object.
         %
         %  input:
         %     obj      -- TuKcons class object.
         %     clv      -- p_TuKcons class object.
         %
           fns = properties(clv);
              for i=1:13
               try
                  obj.(fns{i}) = clv.(fns{i});
               catch
               end
              end
         end



         function Cobj = copyConstructor(obj)
         % copy constructor: obj to struct
          if isa(obj,'TuKcons')
             fns = properties(obj);
             for i=1:length(fns)
              Cobj.(fns{i}) = obj.(fns{i});
             end
          else
             error('Wrong class object! It must be TuKcons!');
          end
         end



       end




 
end
