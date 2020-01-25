classdef TuACore < TuSol
% TUACORE is a subclass object of TUSOL to perform several computations for retrieving 
% and modifying game data. It stores relevant game information needed to draw the 
% anti-core by overloading functions.
%
% Usage: clv = TuACore(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- TuACore class object (subclass of TuSol).
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
% TuCore properties:
%  tu_acrv      -- stores the anti-core vertices of the game (cdd).
%  tu_acddv     -- stores the anti-core vertices of the game (cddmex).
%  tu_imp       -- stores the anti-imputation vertices of the dual game (cdd).
%  tu_cddi      -- stores the anti-imputation vertices of the dual game (cddmex).
%  tu_zov       -- stores the zero-one normalized game.
%
%  
% Properties inherihed from the superclass TuSol
%  tu_prk       -- stores a pre-kernel element of the zero-one normalized game.
%  tu_prk2      -- stores a second pre-kernel element instead of the pre-nucleolus.
%  tu_prn       -- stores the pre-nucleolus of the zero-one normalized game.
%  tu_sh        -- stores the Shapley value of the zero-one normalized game.
%  tu_tauv      -- stores the Tau value of the zero-one normalized game.
%  tu_bzf       -- stores the Banzhaf value of the zero-one normalized game.
%  tu_aprk      -- stores the anti-prek-kernel of the zero-one normalized game.
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
% TuACore methods:
%  TuACore          -- creates the class object TuACore.
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
%   11/11/2012        0.3             hme
%   07/15/2015        0.7             hme
%



    properties(SetObservable = true, SetAccess = public, Dependent = false)
      tu_acrv;
      tu_acddv;
      tu_imp;
      tu_cddi;
      tu_zov;
    end

    methods
       function obj = TuACore(w,gtype,gformat)
       % TUACORE creates the subclass object TuACore to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information needed to draw the 
       % anti-core by overloading functions. 
       %
       % Usage: clv = TuACore(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- TuACore class object (subclass of TuSol).
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
                 error('Game has not the correct size! Must be a three or four person game only!');
              elseif N < 7  
                 error('Game has not the correct size! Must be a three or four person game only!');
              elseif N > 15
                 error('Game has not the correct size! Must be a three or four person game only!');
              elseif r1~=1
                 error('Game has to be a row vector of length 7 or 15.');
              end 
              zov=ZeroOne_Normalization(w);
              gtype = 'acr';
              gformat = 'mattug';
           elseif nargin < 3
              [r1, N]=size(w);
              [~, n]=log2(N);
              if N~=2^n-1;
                 error('Game has not the correct size! Must be a three or four person game only!');
              elseif N < 7
                 error('Game has not the correct size! Must be a three or four person game only!');
              elseif N > 15
                 error('Game has not the correct size! Must be a three or four person game only!');
              elseif r1~=1
                 error('Game has to be a row vector of length 7 or 15.');
              end
              zov=ZeroOne_Normalization(w);
              if isempty(gtype)
                 gtype = 'acr'; % default
              end
              gformat = 'mattug';
           else
              [r1, N]=size(w);
              [~, n]=log2(N);
              if N~=2^n-1;
                 error('Game has not the correct size! Must be a three or four person game only!');
              elseif N < 7
                 error('Game has not the correct size! Must be a three or four person game only!');
              elseif N > 15
                 error('Game has not the correct size! Must be a three or four person game only!');
              elseif r1~=1
                 error('Game has to be a row vector of length 7 or 15.');
              end
              zov=ZeroOne_Normalization(w);
              if isempty(gtype)
                 gtype = 'acr'; % default
              end
           end
           obj = obj@TuSol(w,gtype,gformat);
           obj = setAllSolutions(obj); 
           obj.tu_acrv=AntiCoreVertices(obj);
           obj.tu_acddv=CddAntiCoreVertices(obj);
           obj.tu_imp=AntiImputationVertices(obj);
           obj.tu_cddi=CddAntiImputationVertices(obj);
           obj.tu_zov=zov;
       end
      


    end
end
