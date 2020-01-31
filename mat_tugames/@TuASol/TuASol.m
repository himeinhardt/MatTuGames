classdef TuASol < TuGame
% TUASOL is a subclass object of TUGAME to perform several computations for retrieving 
% and modifying game data. It stores relevant game information and solutions needed by
% overloading functions.
%
% Usage: clv = TuASol(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- TuASol class object (subclass of TuGame).
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
% TuASol properties:
%  tu_aprk      -- stores an anti pre-kernel element.
%  tu_aprk2     -- stores a second anti pre-kernel element instead of the anti pre-nucleolus.
%  tu_aprn      -- stores the anti pre-nucleolus.
%  tu_sh        -- stores the Shapley value.
%  tu_tauv      -- stores the Tau value.
%  tu_bzf       -- stores the Banzhaf value.
%  tu_prk       -- stores the pre-kernel. 
%  tu_tol       -- stores tolerance value. Default is 10^6*eps.
%  aprk_valid   -- returns 1 if tu_aprk stores an anti pre-kernel element, otherwise 0.
%  aprk2_valid  -- returns 1 if tu_aprk2 stores an anti pre-kernel element, otherwise 0.
%  aprn_valid   -- returns 1 if tu_aprn stores the anti pre-nucleolus, otherwise 0.
%  prk_valid    -- returns 1 if tu_prk stores a pre-kernel element, otherwise 0.
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
% TuASol methods:
%  TuASol               -- creates the class object TuASol.
%  setAllAntiSolutions  -- sets all anti solutions listed below to the class object TuASol. 
%  setAPreKernel        -- sets an anti pre-kernel element to the class object TuASol.
%  setAPreNuc           -- sets the anti pre-nucleolus to the class object TuASol.
%  setShapley           -- sets the Shapley value to the class object TuASol.
%  setTauValue          -- sets the Tau value to the class object TuASol.
%  setBanzhaf           -- sets the Banzhaf value to the class object TuASol.
%  setPreKernel         -- sets a pre-kernel element to the class object TuASol.
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
%   08/10/2016        0.9               hme
%




    properties(SetObservable = true, SetAccess = public, Dependent = false)
      tu_aprk;
      tu_aprn;
      tu_aprk2;
      tu_sh;
      tu_tauv;
      tu_bzf;
      tu_prk;
      tu_tol=10^6*eps; 
    end

    properties(GetAccess = 'public', SetAccess = 'private')
      aprk_valid = false;
      aprn_valid = false;
      aprk2_valid = false;
      prk_valid = false;
    end

    methods
       function obj = TuASol(w,gtype,gformat)
       % TUASOL creates the subclass object TuSol to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information and solutions needed by
       % overloading functions.
       %
       % Usage: clv = TuASol(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- TuASol class object (subclass of TuGame).
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
              gtype = '';
              gformat = 'mattug';
           elseif nargin < 3
              gformat = 'mattug';
           else
           end
           obj = obj@TuGame(w,gtype,gformat); 
       end
      
         
         function obj = setAPreKernel(obj,sol) 
         % SETAPREKERNEL sets an anti pre-kernel element to the class object TuASol.
         % 
         %  Usage: clv = setAPreKernel(clv,sol)
         %
         %  output:
         %    clv       -- TuASol class object.
         %
         %  input:
         %     clv      -- TuASol class object.
         %     sol      -- an anti pre-kernel element (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_aprk)
                   obj.tu_aprk = [];   % Calls setter set.tu_aprk.
               end
             else
               obj.tu_aprk = sol;      % Calls setter set.tu_aprk.
             end
         end


        function obj = set.tu_aprk(obj,pt)
          if isempty(pt)
             sol = Anti_PreKernel(obj);
             tol = obj.tu_tol;
             pkQ = Anti_PrekernelQ(obj,sol,tol);
             if pkQ == 1
               obj.tu_aprk = sol;
               obj.aprk_valid = checkAPreKernel(obj);
             else
                obj.tu_aprk = Anti_PreKernel(obj);
                obj.aprk_valid = checkAPreKernel(obj);
             end
          else
             obj.tu_aprk = pt;
             obj.aprk_valid = checkAPreKernel(obj);
          end
        end


        function obj = set.tu_aprk2(obj,pt)
          if isempty(pt)
             sol = StrategicEquivalentAPrK(obj);
             tol = obj.tu_tol;
             pkQ = Anti_PrekernelQ(obj,sol,tol);
             if pkQ == 1
               obj.tu_aprk2 = sol;
               obj.aprk2_valid = checkAPreKernel2(obj);
             else
                obj.tu_aprk2 = StrategicEquivalentAPrK(obj);
                obj.aprk2_valid = checkAPreKernel2(obj);
             end
          else
             obj.tu_aprk2 = pt;
             obj.aprk2_valid = checkAPreKernel2(obj);
          end
        end


         function obj = setAPreNuc(obj,sol)
         % SETAPRENUC sets the anti pre-nucleolus to the class object TuASol.
         %
         %  Usage: clv = setAPreNuc(clv,sol)
         %
         %  output:
         %    clv       -- TuASol class object.
         %
         %  input:
         %     clv      -- TuASol class object.
         %     sol      -- the anti pre-nucleolus (optional).
         %
              if nargin < 2
                 if isempty(obj.tu_aprn)
                    obj.tu_aprn = [];   % Calls setter set.tu_aprn.
                 end
              else
                  obj.tu_aprn = sol;    % Calls setter set.tu_aprn.
              end
         end 


        function obj = set.tu_aprn(obj,pt)
          if isempty(pt)
             try
               sol = cplex_AntiPreNucl(obj);
               obj.tu_aprn = sol;
               tol = obj.tu_tol;
               obj.aprn_valid = Anti_balancedCollectionQ(obj,sol,tol);
             catch
               sol= Anti_PreNucl(obj);
               obj.tu_aprn = sol;
               tol = obj.tu_tol;
               obj.aprn_valid = Anti_balancedCollectionQ(obj,sol,tol);
             end
          else
             obj.tu_aprn = pt;
             obj.aprn_valid = Anti_balancedCollectionQ(obj,sol,tol);
          end
        end


         function obj = setShapley(obj,sol)
         % SETSHAPLEY sets the Shapley value to the class object TuASol.
         %
         %  Usage: clv = setShapley(clv,sol)
         %
         %  output:
         %    clv       -- TuASol class object.
         %
         %  input:
         %     clv      -- TuASol class object.
         %     sol      -- the Shapley value (optional).
         %
                if nargin < 2
                  if isempty(obj.tu_sh)
                     obj.tu_sh = [];
                  end
                else
                  obj.tu_sh = sol;
                end
         end


        function obj = set.tu_sh(obj,pt)
          if isempty(pt)
             obj.tu_sh = ShapleyValue(obj);
          else
             obj.tu_sh = pt;
          end
        end


         function obj = setTauValue(obj,sol)
         % SETTAUVALUE sets the Tau value to the class object TuASol.
         %
         %  Usage: clv = setTauValue(clv,sol)
         %
         %  output:
         %    clv       -- TuASol class object.
         %
         %  input:
         %     clv      -- TuASol class object.
         %     sol      -- the Tau value (optional).
         %
                if nargin < 2
                  if isempty(obj.tu_tauv)
                     obj.tu_tauv = [];
                  end
                else
                  obj.tu_tauv = sol;
                end
         end


        function obj = set.tu_tauv(obj,pt)
          if isempty(pt)
             try
               obj.tu_tauv = TauValue(obj);
             catch
               obj.tu_tauv = [];
             end
          else
             obj.tu_tauv = pt;
          end
        end



         function obj = setBanzhaf(obj,sol)
         % SETBANZHAF sets the Banzhaf value to the class object TuASol.
         %
         %  Usage: clv = setBanzhaf(clv,sol)
         %
         %  output:
         %    clv       -- TuASol class object.
         %
         %  input:
         %     clv      -- TuASol class object.
         %     sol      -- the Banzhaf value (optional). 
         %

                if nargin < 2 
                  if isempty(obj.tu_bzf)
                     obj.tu_bzf = [];
                  end
                else
                  obj.tu_bzf = sol;
                end

         end


        function obj = set.tu_bzf(obj,pt)
          if isempty(pt)
             obj.tu_bzf = banzhaf(obj);
          else
             obj.tu_bzf = pt;
          end
        end



         function obj = setPreKernel(obj,sol)
         % SETPREKERNEL sets a pre-kernel element to the class object TuASol.
         %
         %  Usage: clv = setPreKernel(clv,sol)
         %
         %  output:
         %    clv       -- TuASol class object.
         %
         %  input:
         %     clv      -- TuASol class object.
         %     sol      -- a pre-kernel element (optional).
         %
               if nargin < 2
                  if isempty(obj.tu_prk)
                       obj.tu_prk = [];
                  end
               else isempty(sol)
                  obj.tu_prk = sol;
               end         
         end


        function obj = set.tu_prk(obj,pt)
          if isempty(pt)
             sol = PreKernel(obj);
             tol = obj.tu_tol;
             pkQ = PrekernelQ(obj,sol,tol);
             if pkQ == 1
               obj.tu_prk = sol;
               obj.prk_valid = checkPreKernel(obj);
             else
                obj.tu_prk = PreKernel(obj);
                obj.prk_valid = checkPreKernel(obj);
             end
          else
             obj.tu_prk = pt;
             obj.prk_valid = checkPreKernel(obj);
          end
        end


         function obj = setAllAntiSolutions(obj)
         % SETALLANTISOLUTIONS sets a set of game solutions to the class object TuASol.
         %
         %  Usage: clv = setAllAntiSolutions(clv)
         %
         %  output:
         %    clv       -- TuASol class object.
         %
         %  input:
         %     clv      -- TuASol class object.
         %
                %
                % Anti Pre-Kernel
                  solpk = Anti_PreKernel(obj);
                  tol = obj.tu_tol;
                  pkQ = Anti_PrekernelQ(obj,solpk,tol);
                  if pkQ == 1
                     obj.tu_aprk = solpk;
                  else
                     obj.tu_aprk = solpk;
                     obj.aprk_valid = checkAPreKernel(obj);
                  end
                %
                % Anti Pre-Nucleolus
                  if obj.tuplayers <= 15
                     obj.tu_aprn = []; 
                  end
                %
                % Shapley value
                  obj.tu_sh = ShapleyValue(obj);
                %
                % PreKernel
                  asolpk = PreKernel(obj);
                  tol = obj.tu_tol;
                  apkQ = PrekernelQ(obj,asolpk,tol);
                  if apkQ == 1
                     obj.tu_prk = asolpk;
                  else
                     obj.tu_prk = asolpk;
                     obj.prk_valid = checkPreKernel(obj);
                  end
                %
                % Banzhaf value
                  if strcmp(obj.tutype,'sv')
                   obj.tu_bzf = banzhaf(obj);
                  else
                % Tau value
                   try
                      obj.tu_tauv = TauValue(obj);
                   catch
                      obj.tu_tauv = []; 
                   end
                  end

         end


      end


    methods(Access = 'private')

       function pkQ = checkAPreKernel(obj)
           if isempty(obj.tu_aprk)
              return;
           end
           tol = obj.tu_tol;
           pkQ=Anti_PrekernelQ(obj,obj.tu_aprk,tol);


       end


       function pkQ = checkAPreKernel2(obj)
           if isempty(obj.tu_aprk2)
              return;
           end
           tol = obj.tu_tol;
           pkQ=Anti_PrekernelQ(obj,obj.tu_aprk2,tol);


       end



       function apkQ = checkPreKernel(obj)
           if isempty(obj.tu_prk)
              return;
           end
           tol = obj.tu_tol;
           apkQ=PrekernelQ(obj,obj.tu_prk,tol);


       end

    end







end
