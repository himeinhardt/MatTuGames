function ABd=GetAll_Borel_dividends(clv,idxfm,method,tol)
% GETAll_BOREL_DIVIDENDS computes all possible arrays of Borel dividends from an index function matrix.
% For n>14 this function needs some time to complete.
%
% Usage: ABd=clv.GetAll_Borel_dividends(idxfm,method)
% Define structure variables:
%  output:
%  bd.dv1    -- Borel dividends w.r.t. similar basis v_sp1 that spans the game space.
%  bd.dv2    -- Borel dividends w.r.t. similar basis v_sp2 that spans the game space.
%  bd.BGbs   -- Borel bases obtained from two (T,k) intermediate linear game bases.
%  bd.rp1Q   -- Cross check if the first dividends are consistent with the default game v.
%               Return values are one (true), or zero (false).
%  bd.rp2Q   -- Cross check if the second dividends are consistent with the default game v.
%               Return values are one (true), or zero (false).
%  v         -- The default (input) TU-game of length 2^n-1.
%
%  input:
%  clv      -- TuGame class object.
%  idxfm    -- An index function matrix s.t. for each row we have with i(1)=1 and i(k)=i(k-1) or i(k-1) + 1 or i(k-1) -1.
%              For instance, for n=3 we can set for row one: idxfm(1,:)=[1 2 1].
%  method   -- A string to format the matrix. Permissible methods
%                to format the matrices are 'full','sparse' or the
%                empty string '', to invoke the default, which is 'full' for n < 8
%                otherwise 'sparse'.
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.
%
% Example:
% Define a bankruptcy game by
% d_vec=[27,65,88,43];
% E=148.6667;
% bv=bankruptcy_game(E,d_vec);
% 
% Then define the index function matrix by
% 
% >> idxfm
%
% idxfm =
%
%     1     1     2     3
%     1     1     1     2
%     1     1     2     2
%     1     2     2     2
%     1     2     2     3
%     1     2     3     3
%     1     2     3     2
%     1     1     1     1
%     1     2     3     4
%
% Then call
% >> ABd=GetAll_Borel_dividends(bv,idxfm)
% Select the pair (1,3), for instance, then we get
%
% >> ABd(1,3)
%
% ans = 
%
%  struct with fields:
%
%     dv1: [-44.7222 -20.1111 13.9815 -0.2222 10.8704 -2.5185 84.2315 -34.0556 19.8704 8.6481 75.2315 5.5370 78.3426 89.5648 -170.7500]
%     dv2: [-75.6944 -51.0833 24.3056 -31.1944 21.1944 7.8056 56.4815 -65.0278 30.1944 18.9722 47.4815 15.8611 50.5926 61.8148 75.3704]
%    BGbs: [1x1 struct]
%    rp1Q: 1
%    rp2Q: 1
%       v: [0 0 17.6667 13.6667 40.6667 78.6667 105.6667 0 0 33.6667 60.6667 56.6667 83.6667 121.6667 148.66
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/08/2020        1.9             hme
%                

N=clv.tusize;
n=clv.tuplayers;

if nargin<5
 tol=10^6*eps;
 if n>7
   method='sparse';
  else 
   method='full';
 end
end

[s1,~]=size(idxfm);
for ii=1:s1
    for jj=1:s1
        if ii~=jj
          ABd(ii,jj)=clv.get_Borel_dividends(idxfm(ii,:),idxfm(jj,:),method);
        end
    end
end

