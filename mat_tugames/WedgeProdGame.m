function v=WedgeProdGame(strc,n)
%WedgeProdGame computes from a MC-nets (concise) representation a cooperative production game (wedge production game).
%
% Source: Ieong, Samuel and Shoham, Yoav (2005), Marginal Contribution Nets: A Compact Representation Scheme for Coalitional Games, URL: https://doi.org/10.1145/1064009.1064030
%
%
%  Usage: v=WedgeProdGame(strc,n)
%
%    
% Define variables:
%  output:
%  v        -- A Tu-Game v of length 2^n-1 (wedge production game). 
%
%  input:    
%  strc     -- Rules set represented as a structure array.
%  n        -- Specifies the number of players involved in the game, must be an integer. 
%
%
%    
% Example:
% Let us consider a rule set given by 6 rules as follows
%    
% (1.) x1 \land x2 \land x3 -> 5
% (2.) x1 \land x4 -> 3
% (3.) x2 \land ~x3 -> -2
% (4.) x2 -> 1
% (5.) x3 -> 2
% (6.) x1 \land x3 \land ~x4 -> -3
%
% Then translating the rule set into a structure array through 
%    
% RulesSet(6)=struct('PositiveLiterals',[1 3],'NegativeLiterals',[4],'Value',-3);
% RulesSet(5)=struct('PositiveLiterals',[3],'NegativeLiterals',[],'Value',2);
% RulesSet(4)=struct('PositiveLiterals',[2],'NegativeLiterals',[],'Value',1);
% RulesSet(3)=struct('PositiveLiterals',[2],'NegativeLiterals',[3],'Value',-2);
% RulesSet(2)=struct('PositiveLiterals',[1 4],'NegativeLiterals',[],'Value',3);
% RulesSet(1)=struct('PositiveLiterals',[1 2 3],'NegativeLiterals',[],'Value',5);    
%
% Setting in the next step the number of players to     
% n=4;
%
% It allows a decoding of the MC-nets into their unique integer representation given by      
%
% v=WedgeProdGame(RulesSet,4)    
%    0   -1   -1    2   -1    3    5    0    3   -1    2    2    5    3   11    
%
%
% Example 2:
% Alternatively, the rules set can be represented as a cell array like
% clm={{[1 2 3],[],[5]},{[1 4],[],3},{[2],[3],-2},{[2],[],1},{[3],[],2},{[1 3],[4],-3}};
%
% Again set n=4    
%
% If the cell array respect the above format, the function can be called in connection with the integer by    
%     
% v2=WedgeProdGame(clm,4)    
%    0   -1   -1    2   -1    3    5    0    3   -1    2    2    5    3   11
%    

    
%    
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/26/2023        1.9.1           hme
%    
narginchk(2,2); % check for legal number of input arguments.
v=ReverseMCNetsRep(strc,n);
