classdef p_TuVal < TuGame
% P_TUVAL is a subclass object of TUGAME to perform several computations for retrieving 
% and modifying game data. It stores relevant game information and
% fairness values needed by overloading functions.
%
% Usage: clv = p_TuVal(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- p_TuVal class object (subclass of TuGame).
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
% p_TuVal properties:
%  tu_sh        -- stores the Shapley value.
%  tu_ad        -- stores the Aumann-Dreze value.
%  tu_ow        -- stores the Owen value.
%  tu_my        -- stores the Myerson value.
%  tu_myus      -- stores the Myerson value w.r.t. a union stable system.
%  tu_psus      -- stores the position value w.r.t. a union stable system.
%  tu_pshs      -- stores the position value w.r.t. a hypergraph system.
%  tu_sl        -- stores the solidarity value.
%  tu_csl       -- stores the coalition solidarity value.
%  tu_slsh      -- stores the solidarity Shapley value.
%  tu_asl      -- stores the solidarity value w.r.t. a priori unions.    
%  tu_cs        -- stores a communication structure a la Myerson (mandatory).
%  tu_ptn       -- stores a partition of the grand coalition.
%  tu_us        -- stores a union stable system.
%  tu_hs        -- stores a hypergraph system.    
%  us_valid     -- stores a true (1), if tu_us is a union stable
%                  system, otherwise false (0).
%  hs_valid     -- stores a true (1), if tu_hs is a hypergraph 
%                  communication system, otherwise false (0).
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
% p_TuVal methods:
%  p_TuVal                   -- creates the class object p_TuVal.
%  p_setAllValues            -- sets all values listed below to the class object p_TuVal. 
%  p_setCoalitionStructures  -- sets the required coalition structures
%                               to the class object p_TuVal.
%  p_setShapley              -- sets the Shapley value to the class object p_TuVal.
%  p_setADvalue              -- sets the Aumann-Dreze value to the class object p_TuVal.
%  p_setOwen                 -- sets the Owen value to the class object p_TuVal.
%  p_setMyerson              -- sets the Myerson value to the class object p_TuVal.
%  p_setMyersonUS            -- sets the Myerson value w.r.t. a union
%                               stable system to the class object p_TuVal.
%  p_setPosition             -- sets the position value w.r.t. a union
%                               stable system to the class object p_TuVal.
%  p_setPositionHS           -- sets the position value w.r.t. a hypergraph
%                               system to the class object p_TuVal.
%  p_setSolidarity           -- sets the solidarity value to the class object p_TuVal.
%  p_setCoalitionSolidarity  -- sets the coalition solidarity value to the class object p_TuVal.
%  p_setSolidarityShapley    -- sets the solidarity Shapely value to the class object p_TuVal.    
%  p_setApuSolidarity        -- sets the solidarity value w.r.t. a
%                             priori unions to the class object p_TuVal.    
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
%   08/02/2013        0.4               hme
%




    properties(SetObservable = true, SetAccess = public, Dependent = false)
      tu_sh;
      tu_ad;
      tu_ow;
      tu_my;
      tu_myus;
      tu_psus;
      tu_pshs;
      tu_sl;
      tu_csl;
      tu_slsh;
      tu_asl;      
      tu_cs;
      tu_ptn;
      tu_us;
      tu_hs;    
    end

    properties(GetAccess = 'public', SetAccess = 'private')
      us_valid = false;
      hs_valid = false;
    end

    methods
       function obj = p_TuVal(w,gtype,gformat)
       % P_TUVAL creates the subclass object p_TuVal to perform several computations for retrieving 
       % and modifying game data. It stores relevant game information and fairness values needed 
       % by overloading functions.
       %
       % Usage: clv = p_TuVal(v,'gtype','gformat')
       %
       % Define variables:
       % output:
       % clv           -- p_TuVal class object (subclass of TuGame).
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
      
         
         function obj = p_setShapley(obj,sol)
         % P_SETSHAPLEY sets the Shapley value to the class object p_TuVal.
         %
         %  Usage: clv = p_setShapley(clv,sol)
         %
         %  output:
         %    clv       -- p_TuVal class object.
         %
         %  input:
         %     clv      -- p_TuVal class object.
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
       
       


         function obj = p_setADvalue(obj,sol) 
         % P_SETADVALUE sets the Aumann-Dreze value to the class object p_TuVal.
         % 
         %  Usage: clv = p_setADvalue(clv,sol)
         %
         %  output:
         %    clv       -- p_TuVal class object.
         %
         %  input:
         %     clv      -- p_TuVal class object.
         %     sol      -- the Aumann-Dreze value (optional).   
         %
             if nargin < 2
               if isempty(obj.tu_ad)
                   obj.tu_ad = [];    % Calls setter set.tu_sh.
               end
             else
               obj.tu_ad = sol;      % Calls setter set.tu_sh.
             end
         end


        function obj = set.tu_ad(obj,pt)
          if isempty(pt)
             obj.tu_ad = p_ADvalue(obj,obj.tu_ptn);
          else
             obj.tu_ad = pt;
          end
        end



         function obj = p_setOwen(obj,sol)
         % P_SETOWEN sets the Owen value to the class object p_TuVal.
         %
         %  Usage: clv = p_setOwen(clv,sol)
         %
         %  output:
         %    clv       -- p_TuVal class object.
         %
         %  input:
         %     clv      -- p_TuVal class object.
         %     sol      -- the Owen value (optional).
         %
                if nargin < 2
                  if isempty(obj.tu_ow)
                     obj.tu_ow = [];
                  end
                else
                  obj.tu_ow = sol;
                end
         end


         function obj = set.tu_ow(obj,pt)
           if isempty(pt)
              obj.tu_ow = p_OwenValue(obj,obj.tu_ptn);
           else
              obj.tu_ow = pt;
           end
         end

         function obj = p_setMyerson(obj,sol)
         % P_SETMYERSON sets the Myerson value to the class object p_TuVal.
         %
         %  Usage: clv = p_setMyerson(clv,sol)
         %
         %  output:
         %    clv       -- p_TuVal class object.
         %
         %  input:
         %     clv      -- p_TuVal class object.
         %     sol      -- the Myerson value (optional).
         %
                if nargin < 2
                  if isempty(obj.tu_my)
                     obj.tu_my = [];
                  end
                else
                  obj.tu_my = sol;
                end
         end


         function obj = set.tu_my(obj,pt)
           if isempty(pt)
              obj.tu_my = p_MyersonValue(obj,obj.tu_cs,'cs');
           else
              obj.tu_my = pt;
           end
         end

         
         function obj = p_setMyersonUS(obj,sol)
         % P_SETMYERSONUS sets the Myerson value w.r.t. a union stable
         % system to the class object p_TuVal.
         %
         %  Usage: clv = p_setMyersonUS(clv,sol)
         %
         %  output:
         %    clv       -- p_TuVal class object.
         %
         %  input:
         %     clv      -- p_TuVal class object.
         %     sol      -- the Myerson value (optional).
         %
                if nargin < 2
                  if isempty(obj.tu_myus)
                     obj.tu_myus = [];
                  end
                else
                  obj.tu_myus = sol;
                end
         end


         function obj = set.tu_myus(obj,pt)
           if isempty(pt)
              obj.tu_myus = p_MyersonValue(obj,obj.tu_us,'us');
           else
              obj.tu_myus = pt;
           end
         end         

         
         function obj = p_setPosition(obj,sol)
         % P_SETPOSITION sets the position value w.r.t. a union
         % stable system to the class object p_TuVal.
         %
         %  Usage: clv = p_setPosition(clv,sol)
         %
         %  output:
         %    clv       -- p_TuVal class object.
         %
         %  input:
         %     clv      -- p_TuVal class object.
         %     sol      -- the position value (optional).
         %
                if nargin < 2
                  if isempty(obj.tu_psus)
                     obj.tu_psus = [];
                  end
                else
                  obj.tu_psus = sol;
                end
         end


         function obj = set.tu_psus(obj,pt)
           if isempty(pt)
              obj.tu_psus = p_PositionValue(obj,obj.tu_us,'us');
           else
              obj.tu_psus = pt;
           end
         end         
         
         function obj = p_setPositionHS(obj,sol)
         % P_SETPOSITIONHS sets the position value w.r.t. a hypergraph
         % system to the class object p_TuVal.
         %
         %  Usage: clv = p_setPositionHS(clv,sol)
         %
         %  output:
         %    clv       -- p_TuVal class object.
         %
         %  input:
         %     clv      -- p_TuVal class object.
         %     sol      -- the position value (optional).
         %
                if nargin < 2
                  if isempty(obj.tu_pshs)
                     obj.tu_pshs = [];
                  end
                else
                  obj.tu_pshs = sol;
                end
         end
         
         function obj = set.tu_pshs(obj,pt)
           if isempty(pt)
              obj.tu_pshs = p_PositionValue(obj,obj.tu_hs,'hs');
           else
              obj.tu_pshs = pt;
           end
         end         
         
         
         
         
         function obj = p_setSolidarity(obj,sol)
         % P_SETSOLIDARITY sets the solidarity value to the class object p_TuVal.
         %
         %  Usage: clv = p_setSolidarity(clv,sol)
         %
         %  output:
         %    clv       -- p_TuVal class object.
         %
         %  input:
         %     clv      -- p_TuVal class object.
         %     sol      -- the solidarity value (optional).
         %
                if nargin < 2
                  if isempty(obj.tu_sl)
                     obj.tu_sl = [];
                  end
                else
                  obj.tu_sl = sol;
                end
         end


         function obj = set.tu_sl(obj,pt)
           if isempty(pt)
              obj.tu_sl = p_SolidarityValue(obj);
           else
              obj.tu_sl = pt;
           end
         end

         function obj = p_setCoalitionSolidarity(obj,sol)
         % P_SETCOALITIONSOLIDARITY sets the coalition solidarity value to the class object p_TuVal.
         %
         %  Usage: clv = p_setCoalitionSolidarity(clv,sol)
         %
         %  output:
         %    clv       -- p_TuVal class object.
         %
         %  input:
         %     clv      -- p_TuVal class object.
         %     sol      -- the coalition solidarity value (optional).
         %
                if nargin < 2
                  if isempty(obj.tu_csl)
                     obj.tu_csl = [];
                  end
                else
                  obj.tu_csl = sol;
                end
         end


         function obj = set.tu_csl(obj,pt)
           if isempty(pt)
              obj.tu_csl = p_CoalitionSolidarity(obj,obj.tu_ptn);
           else
              obj.tu_csl = pt;
           end
         end
         
         function obj = p_setSolidarityShapley(obj,sol)
         % P_SETSOLIDARITYSHAPLEY sets the solidarity shapley value to the class object p_TuVal.
         %
         %  Usage: clv = p_setSolidarityShapley(clv,sol)
         %
         %  output:
         %    clv       -- p_TuVal class object.
         %
         %  input:
         %     clv      -- p_TuVal class object.
         %     sol      -- the solidarity Shapley value (optional).
         %
                if nargin < 2
                  if isempty(obj.tu_slsh)
                     obj.tu_slsh = [];
                  end
                else
                  obj.tu_slsh = sol;
                end
         end
         
         
         function obj = set.tu_slsh(obj,pt)
           if isempty(pt)
              obj.tu_slsh = p_SolidarityShapleyValue(obj,obj.tu_ptn);
           else
              obj.tu_slsh = pt;
           end
         end         
         
         function obj = p_setApuSolidarity(obj,sol)
         % P_SETAPUSOLIDARITY sets the solidarity value w.r.t. a
         % priori unions to the class object p_TuVal.
         %
         %  Usage: clv = p_setApuSolidarity(clv,sol)
         %
         %  output:
         %    clv       -- p_TuVal class object.
         %
         %  input:
         %     clv      -- p_TuVal class object.
         %     sol      -- the solidarity value w.r.t. a priori
         %                  unions (optional).
         %
                if nargin < 2
                  if isempty(obj.tu_asl)
                     obj.tu_asl = [];
                  end
                else
                  obj.tu_asl = sol;
                end
         end


         function obj = set.tu_asl(obj,pt)
           if isempty(pt)
              obj.tu_asl = p_apu_SolidarityValue(obj,obj.tu_ptn);
           else
              obj.tu_asl = pt;
           end
         end

         
         
         function obj = p_setCoalitionStructures(obj,cs,ptn,us,hs)
         % P_SETCOALITIONSTRUCTURES sets the required coalition structures to the class object p_TuVal.
         %
         %
         %  Usage: clv = p_setCoalitionStructures(obj,cs,ptn,us,hs)
         %
         %  output:
         %    clv       -- p_TuVal class object.
         %
         %  input:
         %    clv       -- p_TuVal class object.
         %    cs        -- a communication structure a la Myerson (mandatory).
         %    ptn       -- a partition of the grand coalition.
         %    us        -- a union stable system.
         %    hs        -- a hypergraph system.    
         %
         %
                if nargin < 2
                  if isempty(obj.tu_cs)
                     error('A least a coalition structure must be given!') 
                     obj.tu_cs = [];
                     obj.tu_ptn = [];
                     obj.tu_us = [];                     
                     obj.tu_hs = [];                                       
                  end
                elseif nargin < 3
                  obj.tu_cs=cs; 
                  if isempty(obj.tu_ptn)
                     obj.tu_ptn = [];
                     obj.tu_us = [];                                          
                     obj.tu_hs = [];                                       
                  end
                elseif nargin < 4
                  obj.tu_cs=cs; 
                  obj.tu_ptn=ptn;                      
                  if isempty(obj.tu_us)
                     obj.tu_us = [];
                     obj.tu_hs = [];                                       
                  end
                elseif nargin < 5
                  obj.tu_cs=cs; 
                  obj.tu_ptn=ptn; 
                  obj.tu_us=us;                     
                  if isempty(obj.tu_hs)
                     obj.tu_hs = [];                                       
                  end                
                else
                   obj.tu_cs = cs;
                   obj.tu_ptn = ptn;
                   obj.tu_us = us;                                         
                   obj.tu_hs = hs;                                                         
                end
         end


         function obj = set.tu_cs(obj,cs)
           if isempty(cs)
              error('A least a coalition structure must be given!')
           else
              obj.tu_cs = cs;
           end
         end             
       
         function obj = set.tu_ptn(obj,ptn)
           if isempty(ptn)
              gc=obj.tusize;
              n=obj.tuplayers;
              obj.tu_ptn = PartitionSL(gc,obj.tu_cs,n);
           else
              obj.tu_ptn = ptn;
           end
         end             

         function obj = set.tu_us(obj,us)
           if isempty(us)
              obj.tu_us=p_genUnionStable(obj.tu_cs);
              obj.us_valid = true;
           else
              usQ=union_stableQ(us); 
              if usQ==0
                 obj.tu_us=p_genUnionStable(obj.tu_cs);
                 obj.us_valid = p_union_stableQ(obj.tu_us);
              else    
                obj.tu_us = us;
                obj.us_valid = usQ;
              end  
           end
         end                      
         
         function obj = set.tu_hs(obj,hcs)
           if isempty(hcs)
             obj.hs_valid=false;  
           else
             obj.tu_hs = hcs;  
             obj.hs_valid=hypergraphQ(hcs,obj.tuplayers);
           end  
         end     
         
         
         
         function obj = p_setAllValues(obj)
         % P_SETALLVALUESS sets a set of fairness values to the class object p_TuVal.
         %
         %  Usage: clv = p_setAllValues(clv)
         %
         %  output:
         %    clv       -- p_TuVal class object.
         %
         %  input:
         %     clv      -- p_TuVal class object.
         %
                %
                % Shapley value
                  obj.tu_sh = ShapleyValue(obj);
                % Aumann-Dreze value
                  obj.tu_ad = p_ADvalue(obj,obj.tu_ptn);                  
                % Owen value
                  obj.tu_ow = p_OwenValue(obj,obj.tu_ptn);                  
                % Myerson Value
                  obj.tu_my = p_MyersonValue(obj,obj.tu_cs);                  
                % Myerson Value w.r.t a union stable system
                  obj.tu_myus = p_MyersonValue(obj,obj.tu_us,'us');                                    
                % Position Value w.r.t a union stable system
                  obj.tu_psus = p_PositionValue(obj,obj.tu_us,'us');  
                % Position Value w.r.t a hypergraph system
                  obj.tu_pshs = p_PositionValue(obj,obj.tu_hs,'hs');  
                % Solidarity value
                  obj.tu_sl = p_SolidarityValue(obj);                                    
                % Coalition solidarity value
                  obj.tu_csl = p_CoalitionSolidarity(obj,obj.tu_ptn);                  
                % Solidarity Shapley value
                  obj.tu_slsh = p_SolidarityShapleyValue(obj,obj.tu_ptn);                                    
                % APU solidarity value
                  obj.tu_asl = p_apu_SolidarityValue(obj,obj.tu_ptn);
         end


      end










end
