function PS=getPSgame(n,cf,method,tol)
% GETPSGAME computes a PS game from the PS game basis.
%
% Define variables:
%  output:
% output:
%  PS       -- A structure element with the following contents:
%    psQ    -- Returns true (1) or false (0).
%    Q      -- Returns an array of ones and/or zeros of length n.
%    c      -- The sum of the marginal contribution and its
%              complement for each player. 
%    Mi     -- The matrix of marginal contribution and its
%              complement for each player. Each row must be constant. 
%    sh     -- The Shapley value of the PS game.
%    pnc    -- The pre-nucleolus which must be equal to the Shapley value.
%              It is (1/2)*c.
%    v      -- The constructed PS game.
%    gb     -- The PS game basis.
%    cf     -- The vector of coefficients. 
%
%  input:
%  n        -- The number of player involved.
%  cf       -- Coefficients of the PS game basis (optional). Can 
%              be an empty string '';
%  method   -- A string to format the matrix. Permissible methods
%              to format the matrices are 'full','sparse' or the
%              empty string '', to invoke the default, which is full.
%  tol      -- Tolerance value. Its default value is set to 10^6*eps
%              (optional).
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   04/17/2020        1.9             hme
%
%

if nargin<1
    error('At least the game must be given!');
elseif nargin<2
     cf='';
     tol=10^6*eps;
     method='full';
elseif nargin < 3
     tol=10^6*eps;	
     method='full';
elseif nargin < 4
     tol=10^6*eps;
     if isempty(method)
        method='full';
     end 	
end
gb=PS_GameBasis(n);
[s1,s2]=size(gb);
if strcmp(method,'sparse')
   gb=sparse(gb);
end   
seed=137;
rng(seed,'v5uniform');
if isempty(cf)
   bd=200;	
   cf=randi([1,bd],1,s2);
end
v=cf*gb';
PS=ps_gameQ(v,tol);
PS.v=v;
PS.gb=gb;
PS.cf=cf;
