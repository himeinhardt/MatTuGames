classdef TuAVert < TuSol
% TUAVERT is a subclass object of TUSOL to perform several computations for retrieving 
% and modifying game data. It stores relevant game information needed to compute all 
% vertices of an anti-core/imputation set by overloading functions.
%
% Usage: clv = TuAVert(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- TuAVert class object (subclass of TuSol).
%
% input:
% v             -- A Tu-Game v of length 2^n-1.
% gtype         -- A string to define the game type.
%                    Permissible types are:
%                    -- 'cv' (convex/average-convex, semi-convex).
%                    -- 'cr' game with non-empty core.
%                    -- 'sv' simple game.
%                    -- 'acr' game with non-empty anti-core.
%                    -- ''   empty string. (default)
% gformat       -- A string to define the game format.
%                    Permissible formats are:
%                    -- 'mattug' i.e., unique integer representation to perform computation
%                                under MatTuGames. (default)
%                    -- 'mama'   i.e. generic power set representation, i.e Mathematica.
%
% TuVert properties:
%  tu_acrv       -- stores the anti-core vertices (cdd).
%  tu_acddv      -- stores the anti-core vertices (cddmex).
%  tu_aimp       -- stores the anti-imputation vertices (cdd).
%  tu_acddi      -- stores the anti-imputation vertices (cddmex).
%  tu_accv       -- stores the core cover vertices of the game.
%
%  
% Properties inherihed from the superclass TuSol
%  tu_prk       -- stores a pre-kernel element.
%  tu_prk2      -- stores a second pre-kernel element instead of the pre-nucleolus.
%  tu_prn       -- stores the pre-nucleolus.
%  tu_sh        -- stores the Shapley value.
%  tu_tauv      -- stores the Tau value.
%  tu_bzf       -- stores the Banzhaf value.
%  tu_aprk      -- stores the anti-prek-kernel.
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
% TuAVert methods:
%  TuAVert                      -- creates the class object TuVert.
%  setAllAntiVertices          -- sets all anti-vertices listed below.
%  setAntiCoreVertices         -- sets all anti-core vertices (cdd).
%  setCddAntiVertices          -- sets all anti-core vertices (cddmex).
%  setAntiImpVertices          -- sets the anti-imputation vertices (cdd).
%  setCddAntiImpVertices       -- sets the anti-imputation vertices (cddmex).
%  setCddAntiCoverVertices     -- sets the anti-core cover vertices. 
%
%
%  Methods inherited from the superclass TuSol:
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

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/17/2015        0.7             hme
%



    properties(SetObservable = true, SetAccess = public, Dependent = false)
      tu_acrv;
      tu_acddv;
      tu_aimp;
      tu_acddi;
      tu_accv;
    end

    methods
       function obj = TuAVert(w,gtype,gformat)
       % TUAVERT creates the subclass object TuSol to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information needed to compute all 
       % vertices of an anti-core/imputation set by overloading functions.
       %
       % Usage: clv = TuAVert(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- TuAVert class object (subclass of TuSol).
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
           if nargin > 3
              error('Too many input arguments');
           elseif nargin < 1
              error('Game information must be given as a 2^n-1 vector!');
           elseif nargin < 2
              [r1, N]=size(w);
              [~, n]=log2(N);
              if N~=2^n-1;
                 error('Game has not the correct size!');
              elseif r1~=1
                 error('Game has to be a row vector of length 2^n-1.');
              end 
              gtype = 'cr';
              gformat = 'mattug';
           elseif nargin < 3
              [r1, N]=size(w);
              [~, n]=log2(N);
              if N~=2^n-1;
                 error('Game has not the correct size!');
              elseif r1~=1
                 error('Game has to be a row vector of length 2^n-1.');
              end
              gformat = 'mattug';
           else
              [r1, N]=size(w);
              [~, n]=log2(N);
              if N~=2^n-1;
                 error('Game has not the correct size!');
              elseif r1~=1
                 error('Game has to be a row vector of length 2^n-1.');
              end
           end
           obj = obj@TuSol(w,gtype,gformat);
           obj = setAllSolutions(obj); 
       end
     

         function obj = setAllAntiVertices(obj)
         % SETAllANTIVERTICES sets the anti-core/imputation vertices to the class object TuAVert.
         %
         %  Usage: clv = setAllAntiVertices(obj)
         %
         %  output:
         %    clv       -- TuAVert class object.
         %
         %  input:
         %     clv      -- TuAVert class object.
         %
             if nargin < 2
               if isempty(obj.tu_acrv)
                  obj.tu_acrv=AntiCoreVertices(obj);
                  obj.tu_acddv=CddAntiCoreVertices(obj);
                  obj.tu_aimp=AntiImputationVertices(obj);
                  obj.tu_acddi=CddAntiImputationVertices(obj);
                  obj.tu_accv=CddAntiCoreCoverVertices(obj);
               end
             end
         end


         function obj = setAntiCoreVertices(obj,crv)
         % SETVANTIERTICES sets the anti-core vertices to the class object TuAVert. (cdd)
         %
         %  Usage: clv = setAntiVertices(obj,crv)
         %
         %  output:
         %    clv       -- TuAVert class object.
         %
         %  input:
         %     clv      -- TuAVert class object.
         %     crv      -- anti-core vertices.
         %
             if nargin < 2
               if isempty(obj.tu_acrv)
                  obj.tu_acrv=AntiCoreVertices(obj);
               end
             elseif nargin < 3
               if isempty(obj.tu_acrv)
                  obj.tu_acrv = crv;
               end   
             end
         end         



         function obj = setCddAntiVertices(obj,crv)
         % SETCDDANTIVERTICES sets the anti-core vertices to the class object TuAVert. (cddmex)
         %
         %  Usage: clv = setCddAntiVertices(obj,crv)
         %
         %  output:
         %    clv       -- TuAVert class object.
         %
         %  input:
         %     clv      -- TuAVert class object.
         %     crv      -- anti-core vertices.
         %
             if nargin < 2
               if isempty(obj.tu_acddv)
                  obj.tu_acddv=CddAntiCoreVertices(obj);
               end
             elseif nargin < 3
               if isempty(obj.tu_acddv)
                  obj.tu_acddv = crv;
               end
             end
         end


         function obj = setAntiImpVertices(obj,irv)
         % SETANTIIMPVERTICES sets the anti-imputation vertices to the class object TuAVert. (cdd)
         %
         %  Usage: clv = setAntiImpVertices(obj,irv)
         %
         %  output:
         %    clv       -- TuAVert class object.
         %
         %  input:
         %     clv      -- TuAVert class object.
         %     irv      -- anti imputation vertices (cdd).
         %
             if nargin < 2
               if isempty(obj.tu_aimp)
                  obj.tu_aimp=ImputationVertices(obj);
               end
             elseif nargin < 3
               if isempty(obj.tu_aimp)
                  obj.tu_aimp = irv;
               end
             end
         end


         function obj = setCddAntiImpVertices(obj,irv)
         % SETCDDANTIIMPVERTICES sets the anti imputation vertices to the class object TuAVert. (cddmex)
         %
         %  Usage: clv = setCddAntiImpVertices(obj,irv)
         %
         %  output:
         %    clv       -- TuAVert class object.
         %
         %  input:
         %     clv      -- TuAVert class object.
         %     irv      -- anti-imputation vertices (cddmex).
         %
             if nargin < 2
               if isempty(obj.tu_acddi)
                  obj.tu_acddi=CddImputationVertices(obj);
               end
             elseif nargin < 3
               if isempty(obj.tu_acddi)
                  obj.tu_aimp = irv;
               end
             end
         end


         function obj = setCddCoverVertices(obj,ccv)
         % SETCDDCOVERVERTICES sets the core cover vertices to the class object TuVert. (cddmex)
         %
         %  Usage: clv = setCddCoverVertices(obj,crv)
         %
         %  output:
         %    clv       -- TuVert class object.
         %
         %  input:
         %     clv      -- TuVert class object.
         %     crv      -- core vertices.
         %
             if nargin < 2
               if isempty(obj.tu_ccv)
                  obj.tu_ccv=CddCoreCoverVertices(obj);
               end
             elseif nargin < 3
               if isempty(obj.tu_ccv)
                  obj.tu_ccv = ccv;
               end
             end
         end



    end
end
