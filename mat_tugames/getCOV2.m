function CMAT=getCOV2(obs,n,seed)
% GETCOV2 computes from a sample of observations obs and for n-assets the covariance matrix V of a portfolio with negative risk-return relationship.
%
% Source: Benjamin R. Auer and Tobias Hiller (2019), Cost gap, Shapley, or nucleolus allocation: Which is the best
%         game-theoretic remedy for the low-risk anomaly?
%
% Usage: CMAT=getCOV2(obs,n)
%
% Define field variables:
% output:
%  V        -- Covariance matrix.
%  NV       -- Normalized covariance/correlation matrix.
%  mu       -- Vector of mean values. 
%  std      -- Vector of standard deviations.
%  var      -- Vector of variances.   
%    
% input: 
%  obs      -- Integer to determine the number of random observation for each assets. 
%  n        -- Number of assets (integer value).
%  seed     -- Seed value, default is 137.
%
% Example:
% CMAT=getCOV2(100,4)
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
%   03/12/2021        1.9             hme
%
    
    
if nargin < 1
   obs=100;
   n=4;
   %   a=5; % std not larger than 25%.
   a=10; % std larger than 25%.
   b=25;
   seed=137;
elseif nargin < 2
   n=4; 
    % a=5;
   a=10; % std larger than 25%.    
   b=25; 
   seed=137;   
elseif nargin < 3
    % a=5;
   a=10; % std larger than 25%.    
   b=25; 
   seed=137;
end
rng(seed,'v5uniform');
%% pre-specifying mean values
av = sort(max(a.*randn(n,1) + b,2),'descend');    
%% Generating univariate random data
rdata=zeros(obs,n);
for kk=1:n %% lower return is related to higher risk (low risk puzzle)  
    %% a/2+kk-n increases std for asset kk w.r.t. kk-1
    rdata(:,kk)=(a/2+kk-n).*randn(obs,1) + av(kk);
end 
Sig=cov(rdata);
mu=mean(rdata);
R=mvnrnd(mu,Sig,obs);
smr=sum(R,2);
mvar=var(smr);
smu=mean(smr);
cts=smr-smu;
betfac=zeros(1,n);
rmu=zeros(1,n);
for kk=1:n
        %% computing multivariate covariance to determine beta.
        rmu(kk)=mean(R(:,kk));
        ctr=R(:,kk)-rmu(kk);
        covsr=cts'*ctr/(obs-1);
        %% cov is a matrix (univariate distributions).
        %% betfac(kk)=cov(R(:,kk),smr)/mvar;
        betfac(kk)=covsr/mvar;
end
CMAT.V=cov(R);
iV=inv(diag(sqrt(diag(CMAT.V))));
CMAT.NV=iV*CMAT.V*iV;
CMAT.mu=mean(R);
CMAT.std=std(R);
CMAT.var=var(R);
CMAT.beta=betfac;
