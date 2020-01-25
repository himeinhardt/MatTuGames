classdef p_TuKcons < p_TuSol 
% P_TUKCONS creates the subclass object p_TuKcons to perform several computations for retrieving 
% and modifying game data. It stores relevant game information needed to appy various
% generalized consistency investigations on a proposed solution by overloading functions using Matlab's PCT.
%
% Usage: clv = p_TuKcons(v,'gtype','gformat',K)
%
% Define variables:
% output:
% clv           -- p_TuKcons class object (subclass of TuGame).
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
% p_TuKcons properties:
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
% p_TuKcons methods:
%  p_TuKcons                       -- creates the class object p_TuKcons.
%  p_setKConverse_RGP              -- sets results of k-CRGP to the class object p_TuKcons.
%  p_setKReconfirmation_property   -- sets results of k-RCP to the class object TukCons.
%  p_setKStrConverse_RGP           -- sets results of k-SCRGP to the class object p_TuKcons. 
%  p_setKReducedGameProperty       -- sets results of k-RGP to the class object p_TuKcons.
%                                     k is a positive integer between 2 and n. 
%  p_copyTuSol                     -- copies all solutions from TuSol/p_TuSol to TuRep.
%  p_copyTuKcons                   -- copies k-consistency properties from TuKcons to p_TuKcons.
%
%  Methods inherited from the superclass p_TuSol:
%
%  p_TuSol            -- creates the class object p_TuSol.
%  p_setAllSolutions  -- sets all solutions listed below to the class object p_TuSol.
%  p_setPreKernel     -- sets a pre-kernel element to the class object p_TuSol.
%  p_setPreNuc        -- sets the pre-nucleolus to the class object p_TuSol.
%  p_setShapley       -- sets the Shapley value to the class object p_TuSol.
%  p_setTauValue      -- sets the Tau value to the class object p_TuSol.
%  p_setBanzhaf       -- sets the Banzhaf value to the class object p_TuSol.
%  p_setAntiPreKernel -- sets an anti-pre-kernel element to the class object p_TuSol.
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
         function obj = p_TuKcons(w,gtype,gformat,K)
       % P_TUCKONS creates the subclass object p_TuKcons to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information needed to appy various
       % generalized consistency investigations on a proposed solution by overloading functions using Matlab's PCT. 
       %
       % Usage: clv = p_TuKcons(v,'gtype','gformat',K)
       %
       % Define variables:
       % output:
       % clv           -- p_TuKcons class object (subclass of TuGame).
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
           obj = obj@p_TuSol(w,gtype,gformat);
           obj.tu_K = K; 
         end

         function obj = p_setKConverse_RGP(obj,x,K)
         % P_SETCONVERSE_RGP sets results of CRGP to the class object p_TuCons.
         %
         %  Usage: clv = p_setKConverse_RGP(clv,tu_x,K) 
         %
         %  output:
         %    clv       -- p_TuKcons class object.
         %
         %  input:
         %     clv      -- p_TuKcons class object.
         %     x        -- payoff vector (efficient).
         %     K        -- a positve integer between 2 and n.
         %                 default tu_K=2.
         %        
          if nargin > 3 
             error('Too many input arguments');
          elseif nargin < 2
             if isempty(obj.tu_prk)
                obj = p_setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             end
             obj.tu_kCRGP_PRN=p_k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRN');
             obj.tu_kCRGP_PRK=p_k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRK');
             obj.tu_kCRGP_SHAP=p_k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'SHAP');
          elseif nargin < 3
             if isempty(obj.tu_prk)
                obj = p_setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             else
                if all(size(x)==[1, obj.tuplayers ])
                  obj.tu_x = x;
                else
                  error('Payoff vector has not the correct size! It has to be (1,n).')
                end
             end
             obj.tu_kCRGP_PRN=p_k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRN');
             obj.tu_kCRGP_PRK=p_k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRK');
             obj.tu_kCRGP_SHAP=p_k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'SHAP');
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
             obj.tu_kCRGP_PRN=p_k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRN');
             obj.tu_kCRGP_PRK=p_k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRK');
             obj.tu_kCRGP_SHAP=p_k_Converse_RGP_Q(obj,obj.tu_x,obj.tu_K,'SHAP');
          end
         end


         function obj = p_setKReconfirmation_property(obj,x,K)
         % P_SETRECONFIRMATION_PROPERTY sets results of RCP to the class object p_TuCons.
         %
         %  Usage: clv = p_setKReconfirmation_property(clv,x,K)
         %
         %  output:
         %    clv       -- p_TuKcons class object.
         %
         %  input:
         %     clv      -- p_TuKcons class object.
         %     x        -- payoff vector (efficient).
         %     K        -- a positve integer between 2 and n.
         %                 default K=2.
         %
           if nargin > 3 
              error('Too many input arguments');
           elseif nargin < 2
             if isempty(obj.tu_prk)
                obj = p_setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             end
                obj.tu_kRCP_PRN=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'PRN');
                obj.tu_kRCP_PRK=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'PRK');
                obj.tu_kRCP_SHAP=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'SHAP');
                obj.tu_kRCP_HMS_PN=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PN');
                obj.tu_kRCP_HMS_PK=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PK');
           elseif nargin < 3
             if isempty(obj.tu_prk)
                obj = p_setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             else
                if all(size(x)==[1, obj.tuplayers ])
                  obj.tu_x = x;
                else
                  error('Payoff vector has not the correct size! It has to be (1,n).')
                end
             end
                obj.tu_kRCP_PRN=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'PRN');
                obj.tu_kRCP_PRK=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'PRK');
                obj.tu_kRCP_SHAP=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'SHAP');
                obj.tu_kRCP_HMS_PN=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PN');
                obj.tu_kRCP_HMS_PK=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PK');
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
                obj.tu_kRCP_PRN=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'PRN');
                obj.tu_kRCP_PRK=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'PRK');
                obj.tu_kRCP_SHAP=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'SHAP');
                obj.tu_kRCP_HMS_PN=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PN');
                obj.tu_kRCP_HMS_PK=p_k_Reconfirmation_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PK');
           end
         end

         function obj = p_setKStrConverse_RGP(obj,x,K)
         % P_SETREDUCED_GAME_PROPERTY sets results of RGP to the class object p_TuCons.
         %
         %  Usage: clv = p_setKStrConverse_RGP(clv,x,K)
         %
         %  output:
         %    clv       -- p_TuKcons class object.
         %
         %  input:
         %     clv      -- p_TuKcons class object.
         %     x        -- payoff vector (efficient).
         %     K        -- a positve integer between 2 and n.
         %                 default K=2. 
         %
           if nargin > 3 
             error('Too many input arguments');
           elseif nargin < 2
             if isempty(obj.tu_prk)
                obj = p_setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             end
             obj.tu_kSCRGP_PRN=p_k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRN');
             obj.tu_kSCRGP_PRK=p_k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRK');
             obj.tu_kSCRGP_SHAP=p_k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'SHAP');
           elseif nargin < 3
             if isempty(obj.tu_prk)
                obj = p_setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             else
                if all(size(x)==[1, obj.tuplayers ])
                  obj.tu_x = x;
                else
                  error('Payoff vector has not the correct size! It has to be (1,n).')
                end
             end
             obj.tu_kSCRGP_PRN=p_k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRN');
             obj.tu_kSCRGP_PRK=p_k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRK');
             obj.tu_kSCRGP_SHAP=p_k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'SHAP');
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
             obj.tu_kSCRGP_PRN=p_k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRN');
             obj.tu_kSCRGP_PRK=p_k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'PRK');
             obj.tu_kSCRGP_SHAP=p_k_StrConverse_RGP_Q(obj,obj.tu_x,obj.tu_K,'SHAP');
           end
         end


         function obj = p_setKReducedGameProperty(obj,x,K)
         % P_SETKREDUCEDGAMEPROPERTY sets results of RGP to the class object p_TuCons.
         %
         %  Usage: clv = p_setKReducedGameProperty(clv,x,K)
         %
         %  output:
         %    clv       -- p_TuKcons class object.
         %
         %  input:
         %     clv      -- p_TuKcons class object.
         %     x        -- payoff vector (efficient).
         %     K        -- a positve integer between 2 and n.
         %                 default K=2.
         %
           if nargin > 3
              error('Too many input arguments');
           elseif nargin < 2
             if isempty(obj.tu_prk)
                obj = p_setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             end
                obj.tu_kRGP_PRN=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'PRN');
                obj.tu_kRGP_PRK=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'PRK');
                obj.tu_kRGP_SHAP=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'SHAP');
                obj.tu_kRGP_HMS_PN=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PN');
                obj.tu_kRGP_HMS_PK=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PK');
           elseif nargin < 3
             if isempty(obj.tu_prk)
                obj = p_setAllSolutions(obj);
                obj.tu_x = obj.tu_prk;
             else
                if all(size(x)==[1, obj.tuplayers ])
                  obj.tu_x = x;
                else
                  error('Payoff vector has not the correct size! It has to be (1,n).')
                end
             end
                obj.tu_kRGP_PRN=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'PRN');
                obj.tu_kRGP_PRK=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'PRK');
                obj.tu_kRGP_SHAP=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'SHAP');
                obj.tu_kRGP_HMS_PN=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PN');
                obj.tu_kRGP_HMS_PK=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PK');
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
                obj.tu_kRGP_PRN=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'PRN');
                obj.tu_kRGP_PRK=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'PRK');
                obj.tu_kRGP_SHAP=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'SHAP');
                obj.tu_kRGP_HMS_PN=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PN');
                obj.tu_kRGP_HMS_PK=p_k_Reduced_game_propertyQ(obj,obj.tu_x,obj.tu_K,'HMS_PK');
           end
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
       

         function obj = p_copyTuKcons(obj,clv)
         % copy constructor: class object clv to class object obj.
         %
         %  Usage: obj = p_copyTuKcons(obj,clv)
         %
         %  output:
         %    obj       -- p_TuKcons class object.
         %
         %  input:
         %     obj      -- p_TuKcons class object.
         %     clv      -- TuKcons class object.
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
          if isa(obj,'p_TuKcons')
             fns = properties(obj);
             for i=1:length(fns)
              Cobj.(fns{i}) = obj.(fns{i});
             end
          else
             error('Wrong class object! It must be p_TuKcons!');
          end
         end





       end

 
end
