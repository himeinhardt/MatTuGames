function [sv,wC]=tricameral_assembly(n,b,c)
% TRICAMERAL_ASSEMBLY computes from a set of parameters a simple game. 
%
%  Source: P. Dubey and L. S. Shapley (1954), Mathematical Properties of the Banzhaf Power Index.
%
% Usage: sv=tricameral_assembly(n,b,c)
%
% Define variables:
%  output:
%  sv       -- Simple game in length of 2^n-1.
%  wC       -- The list of winning coalitions. 
%
%  input:
%  n        -- number of players n=a+b+c, whereas a=1 (veto player)
%  b        -- Number of players in chamber B, default is b=3.
%  c        -- Number of players in chamber C, default is c=5;
%
%
% Example:
% Take the example from page 103 of the above references. There are three chambers A,B, and C
% with 1,3,5 players. The winning coalitions are those that include a majority of every  
% chamber. Hence, the quorum for chamber A is one, for B it is 2 and for C it is 3.
% Then call 
%
% [sv,wC]=tricameral_assembly(9,3,5)
%
% to get the simple game and the list of winning coalitions. We provide only the set of
% winning coalitions. 
%
% wC =
%
% Columns 1 through 35:
%
%   119   123   125   127   183   187   189   191   215   219   221   223   231   235   237   239   247   251   253   255   311   315   317   319   343   347   349   351   359   363   365   367   375   379   381
%
% Columns 36 through 64:
%
%   383   407   411   413   415   423   427   429   431   439   443   445   447   455   459   461   463   471   475   477   479   487   491   493   495   503   507   509   511 
%
% 
    
    
%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/10/2020        1.9             hme
%

if nargin<1
   n=9;
   a=1;
   b=3;
   c=5;   
else
   a=1;	
end
	
if sum([a,b,c])~=n
   error("The game is not tricameral");
%elseif a~=1
%   error("Chamber A has non-veto players");	
end
%% majority quorum for chambers B and C.
md2=mod([b,c],2);
if md2(1)==1
   thb=ceil(b/2);
elseif md2(1)==0
   thb=ceil(b/2)+1;
end

if md2(2)==1
   thc=ceil(c/2);
elseif md2(2)==0
   thc=ceil(c/2)+1;
end
%% Setting game parameters.
N=2^n-1;
S=1:N;
%A=1:a;
B=a+1:a+b;
C=b+2:b+c+1;
%% Selecting coalitions that contain veto player.
vC=S(bitget(S,a)==1);
lcv=length(vC);
bC=zeros(1,lcv);
sv=zeros(1,N);
%% Selecting coalitions that contain players of B and satisfy the quorum.
for k=1:lcv
   plb=B(bitget(vC(k),B)==1);
   if length(plb)>=thb
      bC(k)=vC(k);	   
   end	   
end
sbC=bC(bC>0);
%% Selecting coalitions that contain players of C and satisfy the quorum.
lbc=length(sbC);
cC=zeros(1,lbc);
for kk=1:lbc
   plc=C(bitget(sbC(kk),C)==1);
   if length(plc)>=thc
      cC(kk)=sbC(kk);
   end
end
%% Winning coalitions.
wC=cC(cC>0);
sv=ismember(S,wC);
