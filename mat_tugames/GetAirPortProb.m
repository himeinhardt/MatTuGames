function [C,nj,rwuidx]=GetAirPortProb(m,cap,a,fm,seed,filename)
% GETAIRPORTPROP computes a pseudo-random airport cost allocation problem of size m+1.
% 
% Usage: [C,nj]=GetAirPortProb(types,cap)
% Define variables:
%  output:
%  C        -- Vector of annual capital costs per types of length m+1.
%  nj       -- Vector of annual plane movements per types of length m+1.
%  rwuidx   -- Vector of runway user index.
%    
%
%  input:
%  m           -- The number of aircraft types, must be an integer. 
%                 The largest type is added, hence m+1 types are considered.    
%  cap         -- Integer number of annual total capital cost.
%  a           -- Minimum of runway usage of aircraft types. Default is a=0.3.
%  fm          -- Maximum number of flight movements of a type. Default is fm=25000.
%                 The range of flight movments per type is in the range [20,fm].     
%  seed        -- A seed number. Default is seed=125.
%  filename    -- Name of the file, csv format is expected. Admissible filenames are
%                 'AirportProb01.csv' with suffix of the file extension or 
%                 'AirportProb01' without suffix of the file extension.   
%  
%  

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/04/2023        1.9.1             hme
% 
    
narginchk(2,6)
    
if nargin<3  
   a=0.3; % minimum of runway usage of types.
   fm=25000; % maximum number of flight movments of a type.
   seed=125;
   filename='AirportProb';
elseif nargin<4
   if isempty(a)
      a=0.3; % minimum of runway usage of types.
   end        
   fm=25000; % maximum number of flight movments of a type.
   seed=125;
   filename='AirportProb';   
elseif nargin<5
   if isempty(a) && isempty(fm) 
      a=0.3; % minimum of runway usage of types.
      fm=25000; % maximum number of flight movments of a type.
   elseif isempty(a) 
      a=0.3; % minimum of runway usage of types.
   end    
   seed=125;
   filename='AirportProb';   
elseif nargin<6
   if isempty(a) && isempty(fm) && isempty(seed)
      a=0.3; % minimum of runway usage of types.
      fm=25000; % maximum number of flight movments of a type.
      seed=125;
   elseif isempty(a) && isempty(fm) 
      a=0.3; % minimum of runway usage of types.
      fm=25000; % maximum number of flight movments of a type.
   elseif isempty(a) 
      a=0.3; % minimum of runway usage of types.
   end    
   filename='AirportProb';
else
   if isempty(a) && isempty(fm) && isempty(seed) && isempty(filename)
      a=0.3; % minimum of runway usage of types.
      fm=25000; % maximum number of flight movments of a type.
      seed=125;
      filename='AirportProb';
   elseif isempty(a) && isempty(fm) && isempty(seed)
      a=0.3; % minimum of runway usage of types.
      fm=25000; % maximum number of flight movments of a type.
      seed=125;
   elseif isempty(a) && isempty(fm) 
      a=0.3; % minimum of runway usage of types.
      fm=25000; % maximum number of flight movments of a type.
   elseif isempty(a) 
      a=0.3; % minimum of runway usage of types.
   end        
end    
rng(seed,'v5uniform');
b=1;   % maximum of runway usage of types. 
rd = a + (b-a).*rand(1,m);
rwuidx=sort(rd,'ascend'); % runway user index
C=cap*rwuidx; % Vector of capital cost per type of length m.
m=m+1;
C(m)=cap; % Vector of capital cost per type of length m+1.
nj=randi([20,fm],1,m); % vector of annual aircarft landings of types.
rwuidx=[rwuidx,1];
if nargin==6
    filename=[filename,'.csv']; 
    writematrix([C;nj;inf(1,m);inf(1,m);rwuidx],filename); % [C;nj;fj;aj;rwuidx]=(saving capital cost, plane movements, actual fee, maintance cost, runway user index)    
else    
    dt = char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss')); 
    filename=[filename,'-',date,'-',dt(end-7:end),'.csv']; 
    writematrix([C;nj;inf(1,m);inf(1,m);rwuidx],filename); 
end    
