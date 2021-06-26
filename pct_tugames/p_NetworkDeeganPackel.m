function [nDPPI,PV]=p_NetworkDeeganPackel(E,th,w)
% P_NETWORKDeeganPackel computes the network Deegan-Packel index from the set of winning coalitions of a network E
% while imposing a threshold of th using Matlab's PCT.
%
% Inspired by the Source: M. Holler and F. Rupp, Power in Networks: The Medici,2020.
%
% Usage: nDPPI=NetworkDeeganPackel(E,th)
% Define variables:
%  output:
%  nPGI     -- Network Public Good Index.
%  PV       -- Public value (PV) of players. 
%
%  input:
%  E        -- An edge matrix of size (lx2) or a cell of numel l.
%              The source must be given by 1, and the sink by the
%              number of the player set. However, the edge matrix
%              can also be of size (lx3) then c can be empty.
%  w        -- A vector of power weights.
%  th       -- Threshold to pass a bill (positive number) greater than 2 and 
%              not greater than or equal to n.
%
% Example:
% Define a matrix of edges given by
% E =
%   1   1   1   1   2   2   3   4
%   2   3   4   5   3   4   4   5
%
% or equivalently
%
% E={[1 2];[1 3];[1 4];[1 5];[2 3];[2,4];[3 4];[4 5]}
% Do not forget to set the separator by semicolon (;) not by comma.
%
% We have 8 edges here. Furthermore, in total we have 5 vertices, and one
% source given by number 1, and a sink by number 5. 
% Set th=3 for the threshold to pass a bill.
% Optinally, one can specify the voting weight of each player by a weighted vector, for instance,
% by w=[2 5 6 2 4]. Otherwise, the default of equal weight is assumed.
%
% Here, we have a player set of {1,2,3,4,5}, hence n=5.
%
% Now, invoke
%  nDPPI=p_NetworkDeeganPackel(E,th)
% to get
%  0.2222    0.1852    0.1852    0.2222    0.1852
%
% For the example with weights vector w and th=10, we get
% nDPPI=p_NetworkDeeganPackel(E,10,w)
%
%  0.1667    0.1944    0.2500    0.1667    0.2222
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/05/2020        1.9             hme
%

if nargin <2
   error('The threshold must be given');
elseif nargin < 3
   [~,mWC,~,n]=p_NetworkMajorityGame(E,th);
else
   [~,mWC,~,n]=p_NetworkMajorityGame(E,th,w);
end

m=length(mWC);
nDPPI=zeros(1,n);
parfor k=1:n;
    mWCk=mWC(bitget(mWC,k)==1);
    sk=length(mWCk);
    A=zeros(sk,1);
    for jj=1:n,
      A(:,jj) = bitget(mWCk,jj);
    end
    PV{k}=A*ones(n,1);
    ssk=sum(1./PV{k});
    nDPPI(k)=ssk/m;
end

