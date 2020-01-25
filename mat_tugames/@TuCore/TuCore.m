classdef TuCore < TuSol
% TUCORE is a subclass object of TUSOL to perform several computations for retrieving 
% and modifying game data. It stores relevant game information needed to draw the 
% core/imputation set by overloading functions.
%
% Usage: clv = TuCore(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- TuCore class object (subclass of TuSol).
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
%  tu_crv       -- stores the core vertices of the game (cdd).
%  tu_cddv      -- stores the core vertices of the game (cddmex).
%  tu_imp       -- stores the imputation vertices of the game (cdd).
%  tu_cddi      -- stores the imputation vertices of the game (cddmex).
%  tu_mgc       -- stores the marginal contributions of players (Weber set).
%  tu_ccv       -- stores the core cover vertices of the game.
%  tu_vals      -- stores the zero-one normalized game.
%  eps_vls      -- stores the epsilon values of the various epsilon games.
%  StrCr        -- stores for the various epsilon games
%      crv_eps    -- the vertices 
%      crst_eps   -- the constraints
%      vol_eps    -- the volumes
%      P          -- the Polyhedrons of all strong epsilon cores.                                                  
%      Mov        -- the frames of the movie.
%
%  
% Properties inherited from the superclass TuSol
%  tu_prk       -- stores a pre-kernel element of the game.
%  tu_prk2      -- stores a second pre-kernel element instead of the pre-nucleolus.
%  tu_prn       -- stores the pre-nucleolus of the game.
%  tu_sh        -- stores the Shapley value of the game.
%  tu_tauv      -- stores the Tau value of the game.
%  tu_bzf       -- stores the Banzhaf value of the game.
%  tu_aprk      -- stores the anti-pre-kernel of the game.
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
% TuCore methods:
%  TuCore           -- creates the class object TuCore.
%  setStrongCores   -- sets the vertices of the strong epsilon cores to class object TuCore. 
%  setAllSolutions  -- sets all solutions listed below to the class object TuCore. 
%  setPreKernel     -- sets a pre-kernel element to the class object TuCore.
%  setPreNuc        -- sets the pre-nucleolus to the class object TuCore.
%  setShapley       -- sets the Shapley value to the class object TuCore.
%  setTauValue      -- sets the Tau value to the class object TuCore.
%  setBanzhaf       -- sets the Banzhaf value to the class object TuCore.
%  setAntiPreKernel -- sets an anti-pre-kernel element to the class object TuCore.
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
%   08/02/2014        0.5             hme
%



    properties(SetObservable = true, SetAccess = public, Dependent = false)
      tu_crv;
      tu_cddv;
      tu_imp;
      tu_cddi;
      tu_mgc;
      tu_ccv;
      tu_vals;
      eps_vls;
      StrCr;
    end

    methods
       function obj = TuCore(w,gtype,gformat)
       % TUCORE creates the subclass object TuCore to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information needed to draw the 
       % core/imputation set by overloading functions. 
       %
       % Usage: clv = TuCore(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- TuCore class object (subclass of TuSol).
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
              gtype = 'cr';
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
                 gtype = 'cr'; % default
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
                 gtype = 'cr'; % default
              end
           end
           obj = obj@TuSol(w,gtype,gformat);
           obj = setAllSolutions(obj); 
           obj.tu_crv=CoreVertices(obj);
           obj.tu_cddv=CddCoreVertices(obj);
           obj.tu_imp=ImputationVertices(obj);
           obj.tu_cddi=CddImputationVertices(obj);
           obj.tu_mgc=CddWeberSet(obj);
           obj.tu_ccv=CddCoreCoverVertices(obj);
           obj.tu_vals=zov;
       end
    

        function obj = setStrongCores(obj,it_stps,crit_val)
        % SETSTRONGCORES sets a family of vertices to the class object TuSol.
        % 
        %  Usage: clv = setStrongCores(clv,it_stps,crit_val)
        %
        %  output:
        %    clv    -- TuCore class object with structure elements:
        %                     
        %       crv_eps  -- set of all strong epsilon core vertices.
        %      crst_eps  -- set of all strong eplison core constraints. 
        %       vol_eps  -- the corresponding core volumes.
        %             P  -- class object Polyhedron.
        %           Mov  -- contains the frames of the movie.
        %
        %  input:
        %     clv       -- TuCore class object.
        %     it_stps   -- number of iteration steps (optional), default is 5.
        %     crit_val  -- a critical strong epsilon value (optional), 
        %                  see functions critical_values*.m   
        %

          if nargin < 2
             it_stps=5;
             crit_val='';
          elseif nargin < 3
             crit_val='';
          end

          fmin=obj.CddLeastCore();
          crQ=obj.CddCoreQ();

          if isempty(crit_val)
             ctv1=obj.critical_value1();
             ctv2=obj.critical_value2();
             ctv3=obj.critical_value_star();
             vc=[ctv1,ctv2,ctv3];
             ctv=max(vc);
             if ctv<=0
               ctv=2;
             end
          else
             ctv=crit_val;
          end
          div=it_stps*25;
          sz=abs(fmin-ctv)/div;
          t=ctv:-sz:fmin;
          t(end+1)=fmin;         
          y=range(obj.tu_cddv);
          [~, idx]=min(y);

          for k=1:length(t)
              obj.eps_vls(k)=t(k);
              v_eps(k,:)=obj.streps_value(t(k));
              [obj.StrCr.crv_eps{k},obj.StrCr.crst_eps{k},obj.StrCr.vol_eps(k),obj.StrCr.P(k)]=CddCoreVertices(v_eps(k,:),idx);
          end

         
        end  


    end
end
