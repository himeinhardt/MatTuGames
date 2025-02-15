function mcnrep=GetMCNetRules(clm)
%GETMCNETRULES transforms a cell array of the MC-nets representation of a TU game into a structure array. 
%    
%
% Source: Ieong, Samuel and Shoham, Yoav (2005), Marginal Contribution Nets: A Compact Representation Scheme for Coalitional Games, URL: https://doi.org/10.1145/1064009.1064030
%
%
%  Usage: mcnrep=GetMCNetRules(clm)
%
%    
% Define variables:
%  output:
%  mcnrep   -- Rules set represented as a structure array.
%
%  input:    
%  clm      -- Rules set represented as a cell array.
%  n        -- Specifies the number of players involved in the game, must be an integer. 
%
%
% Example:
% Let us consider a ruleset given by 6 rules as follows
%    
% (1.) x1 \land x2 \land x3 -> 5
% (2.) x1 \land x4 -> 3
% (3.) x2 \land ~x3 -> -2
% (4.) x2 -> 1
% (5.) x3 -> 2
% (6.) x1 \land x3 \land ~x4 -> -3
%
% Then translating the ruleset into a cell array through  
% clm={{[1 2 3],[],[5]},{[1 4],[],3},{[2],[3],-2},{[2],[],1},{[3],[],2},{[1 3],[4],-3}};
%    
% This cell array can then be transferred into a structure array by calling
% mcnrep=GetMCNetRules(clm) 
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
narginchk(1,1); % check for legal number of input arguments.    
s=numel(clm);
for kk=s:-1:1
    mcnrep(kk)=struct('PositiveLiterals',clm{kk}(1),'NegativeLiterals',clm{kk}(2),'Value',clm{kk}(3));
end    