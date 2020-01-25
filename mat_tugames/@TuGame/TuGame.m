classdef TuGame < hgsetget & matlab.mixin.Copyable
% TUGAME is a class object to perform several computations for retrieving and modifying game data.
% It stores relevant game information needed by overloading functions.
%
% A data array v containing the information of a Tu-game will be converted
% to a class object. By doing so, it will be checked that the data array
% has the correct size and format. These information among others will be 
% stored in some class instances which will be readout by overloading functions.
%
% Usage: clv = TuGame(v,'gtype','gformat')
%
% Define variables:
% output:
% clv           -- TuGame class object.
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
% TuGame properties:
%  tuvalues     -- stores the characteristic values of a Tu-game. 
%  tusize       -- stores the length of the game array/vector.
%  tuplayers    -- stores the number of players involved.
%  tutype       -- stores the game type information ('cv','cr', or 'sv').
%  tuessQ       -- stores if the game is essential.
%  tuformat     -- stores the format how the game is represented.
%  tumv         -- stores the largest value of the game.
%  tumnQ        -- stores the information whether a proper coalition has a higher value 
%                  than the grand coalition 
%  tuSi         -- stores the coalitions having size of n-1.
%  tuvi         -- stores the values of singleton coalitions.
%  tustpt       -- stores a starting point for doing computation. Has lower priority than 
%                  providing a second input argument as an alternative starting point.
%                  Thus, if a starting point is provided as second input argument the 
%                  information stored in tustpt will not be used. 
%
% TuGame methods:
%  TuGame       -- creates the class object TuGame.
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
   tuvalues;
   tusize;
   tuplayers;
   tutype;
   tuessQ;
   tuformat;
   tumv;
   tumnQ;
   tuSi;
   tuvi;
   tustpt;
  end


  methods
     function obj = TuGame(w,gtype,gformat)
     % TUGAME creates the class object TuGame.
     %
     % Usage: clv = TuGame(v,'gtype','gformat').
     %
     % Define variables:
     % output:
     % clv           -- TuGame class object.
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
     %                                 under MatTuGames. (default)
     %                    -- 'mama'   i.e. generic power set representation, i.e Mathematica.
     %
     %
       if nargin > 3 
         error('Too many input arguments');
       elseif nargin < 1
         error('Game information must be given as a 2^n-1 vector!');
       elseif nargin < 2 
         [r1, N]=size(w);
         [~, n]=log2(N);
         if N~=2^n-1
           error('Game has not the correct size! Must be a 2^n-1 vector.');
         elseif r1~=1
           error('Game has to be a row vector of length 2^n-1.');
         end

         mv=max(w);
         mnQ=mv>w(N);
         k=1:n;
         if N==1,
            Si=N;
         else
            Si=bitset(N,k,0);
         end
         obj.tuvalues = w;
         obj.tusize = N;
         obj.tuplayers = n;
         obj.tuvi=w(bitset(0,k));
         obj.tuessQ=sum(obj.tuvi)<=w(N);
         if islogical(w)
            obj.tutype = 'sv';
           else
            obj.tutype = '';
         end
         obj.tuformat = 'mattug';
         obj.tumnQ = mnQ;
         obj.tumv = mv;
         obj.tuSi = Si;
       elseif nargin < 3 
         [r1, N]=size(w);
         [~, n]=log2(N);
         if N~=2^n-1
           error('Game has not the correct size! Must be a 2^n-1 vector.');
         elseif r1~=1
           error('Game has to be a row vector of length 2^n-1.');
         end

         mv=max(w);
         mnQ=mv>w(N);
         k=1:n;
         if N==1,
            Si=N;
         else
            Si=bitset(N,k,0);
         end
         obj.tuvalues = w;
         obj.tusize = N;
         obj.tuplayers = n;
         obj.tuvi=w(bitset(0,k));
         obj.tuessQ=sum(obj.tuvi)<=w(N);
         if islogical(w)
            obj.tutype = 'sv';
           else
            obj.tutype = gtype;
         end
         obj.tuformat = 'mattug';
         obj.tumnQ = mnQ;
         obj.tumv = mv;
         obj.tuSi = Si;
       else
         [r1, N]=size(w);
         [~, n]=log2(N);
         if N~=2^n-1
           error('Game has not the correct size! Must be a 2^n-1 vector.');
         elseif r1~=1
           error('Game has to be a row vector of length 2^n-1.');
         end
         if strcmp(gformat,'mama')
            w=gameToMama(w);
         end
         mv=max(w);
         mnQ=mv>w(N);
         k=1:n;
         if N==1,
            Si=N;
         else
            Si=bitset(N,k,0);
         end
         obj.tuvalues = w;
         obj.tusize = N;
         obj.tuplayers = n;
         obj.tuvi=w(bitset(0,k));
         obj.tuessQ=sum(obj.tuvi)<=w(N);
         if islogical(w)
            obj.tutype = 'sv';
           else
            obj.tutype = gtype;
         end
         obj.tuformat = gformat;
         obj.tumnQ = mnQ;
         obj.tumv = mv;
         obj.tuSi = Si;
       end
     end


     function obj = set.tuvalues(obj,w)
         [r1, N]=size(w);
         [~, n]=log2(N);
         if N~=2^n-1
           error('Game has not the correct size! Must be a 2^n-1 vector.');
         elseif r1~=1
           error('Game has to be a row vector of length 2^n-1.');
         end
         obj.tuvalues = w;
     end

     function obj = set.tuformat(obj,gformat)
        gf = gformat;
        if strcmp(gf,'mattug')
           obj.tuformat = gformat;
        elseif strcmp(gf,'mama')
           obj.tuformat = gformat;
        else
          error('Game format not recognized!')
        end
     end

     function obj = set.tustpt(obj,pt)
          if all(size(pt)==[1, obj.tuplayers ])
              obj.tustpt = pt;
          elseif isempty(pt)
              obj.tustpt = pt; 
          else
              error('Payoff vector has not the correct size! It has to be (1,n).')
          end
     end


     function obj = set.tuessQ(obj,essQ,vi)
         if isempty(obj.tuessQ)
            n=obj.tuplayers;
            k=1:n;
            v=obj.tuvalues;
            obj.tuvi=v(bitset(0,k));
            obj.tuessQ=sum(obj.tuvi)<=v(end);
         elseif isempty(vi)
            n=obj.tuplayers;
            k=1:n;
            v=obj.tuvalues;
            obj.tuvi=v(bitset(0,k));
            obj.tuessQ=essQ;
         else
            obj.tuvi=vi;
            obj.tuessQ=essQ;
         end
     end

     function obj = set.tuvi(obj,vi)
         if isempty(vi)
            n=obj.tuplayers;
            k=1:n;
            v=obj.tuvalues;
            obj.tuvi=v(bitset(0,k));
         else
            obj.tuvi=vi;
         end
     end

    
     function obj = startpt(obj,pt)
     % STARTPT set a starting point to the TuGame class object.
     %
     % Usage: clv = startpt(clv,pt)
     %
     % Define variables:
     % output:
     %    clv          -- TuGame class object.
     %
     % input:
     %    clv          -- TuGame class object.
     %    pt           -- staring point for computation (optional).
       if nargin < 2 
            obj.tustpt = select_starting_pt(obj.tuvalues);
       elseif isempty(pt)
            obj.tustpt = select_starting_pt(obj.tuvalues);
       else
            if all(size(pt)==[1, obj.tuplayers ])
              obj.tustpt = pt;
            else
               error('Payoff vector has not the correct size! It has to be (1,n).')
            end 
       end 
     end
  end
end
