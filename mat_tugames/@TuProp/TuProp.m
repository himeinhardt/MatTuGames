classdef TuProp < TuGame
% TUPROP creates the subclass object TuProp to perform several computations for retrieving 
% and modifying game data. It stores relevant game information, and properties needed by
% overloading functions.
%
% Usage: clv = TuProp(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- TuProp class object (subclass of TuGame).
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
% TuProp properties:
%
%  cv_valid          -- convex
%  acv_valid         -- average-convex
%  scv_valid         -- semi-convex
%  kcv_valid         -- k-convex
%  sad_valid         -- super-additive
%  wsad_valid        -- weakly super-additive
%  mon_valid         -- monotone
%  zmon_valid        -- zero-monotone
%  cr_valid          -- core
%  acr_valid         -- anti-core
%
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
% TuProp methods:
%
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
%   05/27/2013        0.3               hme
%   11/11/2012        0.3  beta         hme
%


      properties (GetAccess = 'public', SetAccess = 'private')
      cv_valid = false;    
      acv_valid = false;  
      scv_valid = false 
      kcv_valid = false; 
      sad_valid = false;
      wsad_valid = false; 
      mon_valid = false;  
      zmon_valid = false; 
      cr_valid = false;
      acr_valid = false; 
      end




     methods
       function obj = TuProp(w,gtype,gformat)
       % TUPROP creates the subclass object TuProp to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information and properties needed by
       % overloading functions and to store game properties like convexity.
       %
       % Usage: clv = TuProp(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- TuProp class object (subclass of TuGame).
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
           obj = obj@TuGame(w,gtype,gformat);
           obj = checkAllProperties(obj);
        end



    end



    methods(Access = 'private')


        function obj = checkAllProperties(obj)
        % CHECKALLPROPERTIES sets all derived game properties to the class object TuProp.
        %
        %  Usage: clv = checkAllProperties(clv)
        %
        %  output:
        %    clv       -- TuProp class object.
        %
        %  input:
        %     clv      -- TuProp class object.
        %
        %
             global CTOLR
             if isempty(CTOLR), CTOLR=10^6*eps; end  % tolerance value for coreQ.
             obj.cv_valid = convex_gameQ(obj);
             obj.acv_valid = average_convexQ(obj);
             obj.scv_valid = semi_convexQ(obj);
             obj.kcv_valid = k_convexQ(obj);
             obj.sad_valid = super_additiveQ(obj);
             obj.wsad_valid = weakly_super_additiveQ(obj);
             obj.mon_valid = monotone_gameQ(obj);
             obj.zmon_valid = zero_monotonicQ(obj);
             try
                obj.cr_valid = CddCoreQ(obj,CTOLR);  % faster than coreQ().
             catch
                obj.cr_valid = coreQ(obj,CTOLR); 
             end
             try
                obj.acr_valid = CddAntiCoreQ(obj,CTOLR);
             catch
                obj.acr_valid = anti_coreQ(obj,CTOLR);
             end
        end



    end


end
