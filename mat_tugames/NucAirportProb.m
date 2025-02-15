function [nc,alp,tk]=NucAirportProb(C,nj)
% NUCAIRPORTPROB computes the nucleolus from a airport capital cost problem.
% 
% Source: S. C. Littlechild. A simple expression for the nucleolus in a special case. International Journal of Game Theory, 3:21–29, 1974b. URL https://doi.org/10.1007/BF01766216. 
%    
% Usage: [nc,alp,tk]=NucAirportProb(C,nj)
% Define variables:
%  output:
%  nc        -- The nucleolus of the associated airport capital cost game.
%  alp       -- Sequence of minmax excess value.
%  tk        -- Sequence of cumulative plane movements over subsets of types.
%
%  input:
%  C        -- Vector of capital costs per types of length m.
%  nj       -- Vector of plane movements per types of length m.
%  

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/29/2023        1.9.1             hme
% 

msg=nargchk(2,2,nargin); % check for legal number of input arguments.
error(msg);

% Pre-allocating output variables.        
m=length(C);
nc=zeros(1,m);
alp=zeros(1,m);
tk=zeros(1,m);
%
% Initializing iteration variables.
mta=0;
cnj=cumsum(nj); % Cumulative plane movements from 1 up to 1,2, till the m-type. 
cnjm=cnj(m); % Cumulative plane movements of all types.
cnj(m)=[]; % Cumulative plane movements from 1 up to 1,2, till the (m-1)-type.
tCrd=C; % Reassigning capital cost vector to a new variable.
Cm=C(m); % Selecting the capital cost of type m, and assigning to a variable.
tCrd(m)=[]; % Deleting capital cost of type m from the list.
%
% Initializing indices.
kk=1;
elm=0;
st=1;
% Solving iteratively the airport cost problem in form of stylized LPs (reduced form of LPs).
while  isempty(tCrd)==0
    % Solving cost allocation problem kk.
    cnj1=cnj+1;
    [C1,idx]=min((tCrd-sum(mta))./cnj1);
    C2=(Cm-sum(mta))/cnjm;
    alp(kk)=min([C1,C2]);     
    elm=idx+elm;
    % Assigning output values of problem kk.
    tk(kk)=cnj(idx); % Taking the cumulative plane movement from type st up to idx.
    nc(st:elm)=alp(kk); % Nucleolus cost distributions from type st up to elm.
    cnj=cumsum(nj(elm+1:end)); % Truncated cumulative plane movements from elm+1 up to elm+1, till the m-type.
    cnjm=cnj(end); % Truncated cumulative plane movements of all types.
    cnj(end)=[]; % Truncated cumulative plane movements from elm+1 up to elm+1 till the (m-1)-type.
    mta(kk)=tk(kk)*alp(kk);
    tCrd=tCrd(idx+1:end); % Singling out capital cost vector from idx+1 up to m-1.
    kk=kk+1;
    st=elm+1;
end    
% Assigning the final values w.r.t. type m.
tk(kk)=cnjm;
C2=(Cm-sum(mta))/cnjm;
alp(kk)=C2;
nc(end)=C2;
alp=-alp; % Changing to the correct sign of the minmax values. 
