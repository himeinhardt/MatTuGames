function [nMCS,mf,x,v]=flow_probMinCut(E,c,ow)
% FLOW_GAMEMinCut computes from a flow problem (E,c,ow) the minimum cut using the optimization toolbox.
% Uses Dual-Simplex (Matlab R2015a). This function requires the
% grTheory - Graph Theory Toolbox that can be found here:
%
% http://www.mathworks.com/matlabcentral/fileexchange/4266-grtheory-graph-theory-toolbox?s_tid=srchtitle
%
% Usage: [nMCS,mf]=flow_probMinCut(E,c,ow)
%
% Define variables:
% output:
%  nMCS     -- The list of the numbers of arrows included 
%              in first minimal cut-set.
%  mf       -- The total flow through each minimal cut set.
%  x        -- The vector of flows of the minimal cut set.
%  v        -- A Tu-Game v of length 2^n-1. 
% input: 
%  E        -- An edge matrix of size (lx2) or a cell of numel l.
%              The source must be given by 0, and the sink by the
%              number of vertices plus one. However, the edge matrix
%              can also be of size (lx3) then c can be empty. This is
%              due to be compatible with the GrTheory toolbox. In this
%              case there is no need to adjust the edge matrix from GrTheory.
%  c        -- A capacity vector of max flow for each edge. Can be set
%              to [] if the edge matrix has size(lx3). The capacity
%              vector will then be obtained from E(:,3).
%  ow       -- An ownership vector of the edges.
%
%
% Example:
% Define a matrix of edges given by
% E =
%   0   0   1   1   1   4   2   3   3   4   5
%   1   2   2   3   4   2   5   4   6   5   6
%
% or equivalently
%
% E={[0 1];[0 2];[1 2];[1 3];[1 4];[4 2];[2 5];[3,4];[3,6];[4,5];[5,6]};
% Do not forget to set the separator by semicolon (;) not by comma.
%
% We have 11 edges here. Furthermore, in total we have 5 vertices, and one
% source given by number 0, and a sink by number 6.
% Then specify the capacity vector of max flow of each edge.
% c =
%   2   2   1   2   1   1   2   1   2   1   3
% Finally, determine to which player the edges are owned
% ow =
%   1   2   3   1   1   2   2   2   3   3   1
%
% Here, we have a player set of {1,2,3}, hence n=3.
%
% Now, invoke
%  [nMCS,mf,v]=flow_probMinCut(E,c,ow)
% to get
% v =
%     0     0     2     0     2     0     4
%


%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/15/2018        0.9             hme
%


if nargin <2
   error('The capacity constraints must be given');
elseif nargin <3
   warning('Simple ownership vector will be constructed! We get as many player as edges.')
   n=length(E);
   ow=1:n;
end
mE=E;
if iscell(E)
   n=length(E);
   E=cell2mat(E);
% if the source in E is given by zero.
   me=min(E(:,1));
   if me(1)==0
      s=1;
      E=E+s;
      t=max(E(:,2));
      E(:,3)=c';
   else
      s=min(E(:,1));
      t=max(E(:,2));
      E(:,3)=c';
   end
else
   n=length(E); 
   [m,cs]=size(E);
   if m==2
      E=E';
      me=min(E(:,1));
      if me(1)==0
         E0=E(:,1:2);
         s=1;
         t=max(E(:,2))+s;
         E0=E0+s;
         E=[E0,c'];
      else
        s=min(E(:,1));
        t=max(E(:,2));
        E(:,3)=c';
      end
   else
      me=min(E(:,1));
      if me(1)==0
         E0=E(:,1:2);
         s=1; 
         t=max(E(:,2))+s;
         E0=E0+s;
         E=[E0,c'];
      else   
        s=min(E(:,1));
        t=max(E(:,2));
        E(:,3)=c';
      end
   end
end


if nargout < 4
   [nMCS,mf]=grMinCutSet(E,s,t); % Here, first call of grMaxFlows.
   cfl=grMaxFlows(E,s,t);  % Not nice! We need to call grMaxFlows again.
   x=zeros(1,n);
   x(nMCS)=cfl(nMCS); %  Must be a core imputation of the game.!!!
else
   if nargin==3
      v=flow_game(mE,c,ow);
   else
      v=grMaxFlowGame(E,s,t,n);
   end
   [nMCS,mf]=grMinCutSet(E,s,t);
end

