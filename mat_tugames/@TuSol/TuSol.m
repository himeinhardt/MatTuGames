classdef TuSol < TuGame
% TUSOL is a subclass object of TUGAME to perform several computations for retrieving 
% and modifying game data. It stores relevant game information and solutions needed by
% overloading functions.
%
% Usage: clv = TuSol(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- TuSol class object (subclass of TuGame).
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
% TuSol properties:
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
% TuSol methods:
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


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/27/2013        0.3               hme
%   08/05/2013        0.4               hme
%   08/03/2016        0.9               hme
%




    properties(SetObservable = true, SetAccess = public, Dependent = false)
      tu_prk;
      tu_prn;
      tu_prk2;
      tu_sh;
      tu_tauv;
      tu_bzf;
      tu_aprk;
    end

    properties(GetAccess = 'public', SetAccess = 'private')
      prk_valid = false;
      prn_valid = false;
      prk2_valid = false;
      aprk_valid = false;
    end

    methods
       function obj = TuSol(w,gtype,gformat)
       % TUSOL creates the subclass object TuSol to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information and solutions needed by
       % overloading functions.
       %
       % Usage: clv = TuSol(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- TuSol class object (subclass of TuGame).
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
      
         
         function obj = setPreKernel(obj,sol) 
         % SETPREKERNEL sets a pre-kernel element to the class object TuSol.
         % 
         %  Usage: clv = setPreKernel(clv,sol)
         %
         %  output:
         %    clv       -- TuSol class object.
         %
         %  input:
         %     clv      -- TuSol class object.
         %     sol      -- a pre-kernel element (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_prk)
                   obj.tu_prk = [];   % Calls setter set.tu_prk.
               end
             else
               obj.tu_prk = sol;      % Calls setter set.tu_prk.
             end
         end


        function obj = set.tu_prk(obj,pt)
          if isempty(pt)
             sol = PreKernel(obj);
             pkQ = PrekernelQ(obj,sol);
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


        function obj = set.tu_prk2(obj,pt)
          if isempty(pt)
             sol = ols_prekernel(obj);
             pkQ = PrekernelQ(obj,sol);
             if pkQ == 1
               obj.tu_prk2 = sol;
               obj.prk2_valid = checkPreKernel2(obj);
             else
                obj.tu_prk2 = ols_prekernel(obj);
                obj.prk2_valid = checkPreKernel2(obj);
             end
          else
             obj.tu_prk2 = pt;
             obj.prk2_valid = checkPreKernel2(obj);
          end
        end


         function obj = setPreNuc(obj,sol)
         % SETPRENUC sets the pre-nucleolus to the class object TuSol.
         %
         %  Usage: clv = setPreNuc(clv,sol)
         %
         %  output:
         %    clv       -- TuSol class object.
         %
         %  input:
         %     clv      -- TuSol class object.
         %     sol      -- the pre-nucleolus (optional).
         %
              if nargin < 2
                 if isempty(obj.tu_prn)
                    obj.tu_prn = [];   % Calls setter set.tu_prn.
                 end
              else
                  obj.tu_prn = sol;    % Calls setter set.tu_prn.
              end
         end 


        function obj = set.tu_prn(obj,pt)
          if isempty(pt)
             try
               solvQ = PreNucl_llp(obj);
               y=0;
             catch
               [solvQ,qq,pp,y]= prenucl(obj);
             end
              if y==Inf
                 pnQ = PrekernelQ(obj,solvQ);
                 if pnQ == 1
                    obj.tu_prn = solvQ;
                    obj.prn_valid = balancedCollectionQ(obj);;
                    obj.tu_prk2 = solvQ;
                    obj.prk2_valid = true;
                 else
                    obj.tu_prn = solvQ;
                    obj.prn_valid = false;
                 end
              else
                   obj.tu_prn = solvQ;
                   pnQ = PrekernelQ(obj,solvQ);
                   if pnQ == 1
                       obj.prn_valid = balancedCollectionQ(obj);
                   else
                       obj.prn_valid = false;
                   end
              end
          else
             obj.tu_prn = pt;
             obj.prn_valid = checkPreNuc(obj);
          end
        end


         function obj = setShapley(obj,sol)
         % SETSHAPLEY sets the Shapley value to the class object TuSol.
         %
         %  Usage: clv = setShapley(clv,sol)
         %
         %  output:
         %    clv       -- TuSol class object.
         %
         %  input:
         %     clv      -- TuSol class object.
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
         % SETTAUVALUE sets the Tau value to the class object TuSol.
         %
         %  Usage: clv = setTauValue(clv,sol)
         %
         %  output:
         %    clv       -- TuSol class object.
         %
         %  input:
         %     clv      -- TuSol class object.
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
         % SETBANZHAF sets the Banzhaf value to the class object TuSol.
         %
         %  Usage: clv = setBanzhaf(clv,sol)
         %
         %  output:
         %    clv       -- TuSol class object.
         %
         %  input:
         %     clv      -- TuSol class object.
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



         function obj = setAntiPreKernel(obj,sol)
         % SETANTIPREKERNEL sets an anti-pre-kernel element to the class object TuSol.
         %
         %  Usage: clv = setAntiPreKernel(clv,sol)
         %
         %  output:
         %    clv       -- TuSol class object.
         %
         %  input:
         %     clv      -- TuSol class object.
         %     sol      -- an anti-pre-kernel element (optional).
         %
               if nargin < 2
                  if isempty(obj.tu_aprk)
                       obj.tu_aprk = [];
                  end
               else isempty(sol)
                  obj.tu_aprk = sol;
               end         
         end


        function obj = set.tu_aprk(obj,pt)
          if isempty(pt)
             sol = Anti_PreKernel(obj);
             apkQ = Anti_PrekernelQ(obj,sol);
             if apkQ == 1
               obj.tu_aprk = sol;
               obj.aprk_valid = checkAntiPreKernel(obj);
             else
                obj.tu_aprk = Anit_PreKernel(obj);
                obj.aprk_valid = checkAntiPreKernel(obj);
             end
          else
             obj.tu_aprk = pt;
             obj.aprk_valid = checkAntiPreKernel(obj);
          end
        end


         function obj = setAllSolutions(obj)
         % SETALLSOLUTIONS sets a set of game solutions to the class object TuSol.
         %
         %  Usage: clv = setAllSolutions(clv)
         %
         %  output:
         %    clv       -- TuSol class object.
         %
         %  input:
         %     clv      -- TuSol class object.
         %
                %
                % Pre-Kernel
                  solpk = PreKernel(obj);
                  pkQ = PrekernelQ(obj,solpk);
                  if pkQ == 1
                     obj.tu_prk = solpk;
                  else
                     obj.tu_prk = solpk;
                     obj.prk_valid = checkPreKernel(obj);
                  end
                %
                % Pre-Nucleolus
                  if obj.tuplayers <= 15
                     obj.tu_prn = []; 
                  end
                %
                % Shapley value
                  obj.tu_sh = ShapleyValue(obj);
                %
                % Anti-PreKernel
                  asolpk = Anti_PreKernel(obj);
                  apkQ = Anti_PrekernelQ(obj,asolpk);
                  if apkQ == 1
                     obj.tu_aprk = asolpk;
                  else
                     obj.tu_aprk = asolpk;
                     obj.aprk_valid = checkAntiPreKernel(obj);
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

       function pkQ = checkPreKernel(obj)
           if isempty(obj.tu_prk)
              return;
           end
           pkQ=PrekernelQ(obj,obj.tu_prk);


       end


       function pkQ = checkPreKernel2(obj)
           if isempty(obj.tu_prk2)
              return;
           end
           pkQ=PrekernelQ(obj,obj.tu_prk2);


       end


       function pnQ = checkPreNuc(obj)
              try
                 [solvQ,qq,pp,y]= prenucl2(obj,obj.tu_prn);
              catch
                 y=0;
              end
              if y==Inf
                 pkQ = PrekernelQ(obj,solvQ);
                 if pkQ == 1
                    pnQ = balancedCollectionQ(obj);
                    obj.tu_prk2 = solvQ;
                    obj.prk2_valid = true;
                 else
                    pnQ = false;
                 end
              else 
                pkQ = PrekernelQ(obj,obj.tu_prn);
                if pkQ == 1
                   pnQ = balancedCollectionQ(obj);
                else
                   pnQ = false;
                 end
              end
       end

       function apkQ = checkAntiPreKernel(obj)
           if isempty(obj.tu_aprk)
              return;
           end
           apkQ=Anti_PrekernelQ(obj,obj.tu_aprk);


       end

    end







end
