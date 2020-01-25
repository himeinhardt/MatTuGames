function CddKernelCatchers(varargin)
% CDDKERNELCATCHERS plots some kernel catchers of game v. 
%
% Requires the Multi-Parametric Toolbox 3
% http://people.ee.ethz.ch/~mpt/3/
%
% Usage: CddKernelCatchers(varargin)
% Define variables:
%  output:
%             -- A plot of the upper set, lower set, and reasonable set of game v.
%                All these set are plotted w.r.t. the core and imputation set.
%
%                Upper set contains the (pre-)kernel under certain conditions.
%                Lower set contains the (pre-)kernel under certain conditions.
%                Both sets could be empty.
%                Reasonable set contains the pre-kernel. 
%
%  input:     At most six input arguments are admissible. Any order of
%             the input arguments below is allowed. Nevertheless, at
%             least a game v must be specified. 
%
%
%  v          -- A Tu-Game v of length 2^n-1.
%  crit_val   -- A number as a string to specify the volume of the reasonable set
%                in relation to the upper set order to plot this set. If the
%                ratio is larger than 15 (default), then the reasonable set
%                will not be plotted. In case that you want to plot the set
%                that have a ratio, let say, of 16 invoke
%
%                CddKernelCatchers(v,'17')                
%
%  reas_sol   -- An integer to draw the reasonable set in connection with all other
%                kernel catchers. The reasonalbe set will not be plotted if it is
%                more than 15 times larger than the upper set.
%                The reasonable set will be drawn if the size is different from the 
%                upper set, otherwise it is ommitted or set a crit_val as above
%                to force a plot.
%
%                Permissible values are 1 (true) or 0 (false).
%                Default is 1 (true). 
%
%  add_sol    -- A string to invoke additional solutions into the plot.
%                Permissible solutions are:
%                'none', this is the default value.
%                'prk', a pre-kernel element will be incorporated.
%                'prn', the pre-nucleolus will be incorporated.
%                'shap', the Shapley value will be incorporated.
%                'all', all three solutions above will be incorporated.
%  vw_pt       -- A string command to determine the view point.
%                The default view point is
%                [120, 25]
%
%                To specify a different view point type, for instance
%                vw_pt='view(130,35)'
%                
%                Then invoke
%
%                CddKernelCatchers(v,'all',vw_pt,'17')
%
%  tol        -- A positive tolerance value. Its default value is set to 10^9*eps.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/13/2014        0.6             hme
%

narginchk(1,6);

v='';
add_sol='';
reas_sol='';
tol='';
vw_pt='';
crit_val='';

for i=1:nargin

  if ischar(varargin{i})
    if strcmp(varargin{i},'none')
       add_sol='none';
    elseif strcmp(varargin{i},'prk')
       add_sol='prk';
    elseif strcmp(varargin{i},'prn')
       add_sol='prn';
    elseif strcmp(varargin{i},'shap')
       add_sol='shap';
    elseif strcmp(varargin{i},'all')
       add_sol='all';
    elseif isempty(regexp(varargin{i},'view'))==0
       vw_pt=varargin{i};
    else
       crit_val=str2num(varargin{i});
       if isempty(crit_val)==0
          crit_val=str2num(varargin{i});
       end
    end
  else 
    if varargin{i}==1
      reas_sol=varargin{i};     
    elseif varargin{i}==0
      reas_sol=varargin{i};
    elseif varargin{i}>0 & varargin{i}<1
      tol=varargin{i};     
    elseif length(varargin{i})==7
      v=varargin{i};
    elseif length(varargin{i})==15
      v=varargin{i};
    else
     N=15;
     if varargin{i}>1
      if length(varargin{i})<N
        error('Game has not the correct size!'); 
      elseif length(varargin{i})>N
        error('Game has dimension larger than four!');
      else
        error(['Input argument at position ' int2str(i) ' not recognized']);
      end
     else
       error(['Input argument at position ' int2str(i) ' not recognized']);
     end
    end
  end
end


if isempty(v)
   error('At least the game must be given!');
end
if isempty(add_sol)
  add_sol='none';
end
crQ=CddCoreQ(v);
core_sol=crQ;
if isempty(tol)
  tol=10^9*eps; 
end
if isempty(vw_pt)
  vw_pt='view(120,25)';
end

N=length(v);
[~, n]=log2(N);

if isempty(crit_val)
   ctv=15;
else
   ctv=crit_val;
end
[impv,~,imp_vol,Pip]=CddImputationVertices(v);
y1=range(impv);
[~,idx]=min(y1);
if core_sol==1
  [~,~,cr_vol,Pv]=CddCoreVertices(v,idx);
end
[~,~,rs_vol,Prea]=CddReasonableSetVertices(v,idx);
[~,~,prs_vol,Pprea]=CddUpperSetVertices(v,idx);
sma=smallest_amount(v);
lwseQ=sum(sma)>v(N);
if lwseQ==1
else
  [~,~,lws_vol,Plws]=CddLowerSetVertices(v,idx);
end
rq_vol=rs_vol/prs_vol;
v_prk=PreKernel(v);
v_prn=CddPrenucl(v);
v_sh=ShapleyValue(v);

if reas_sol==1
  if rq_vol<ctv
    ms1=min(Pip.V);
    ml1=max(Pip.V);
    ms2=min(Prea.V);
    ml2=max(Prea.V);
  else
    ms1=min(Pip.V);
    ml1=max(Pip.V);
    ms2=min(Pprea.V);
    ml2=max(Pprea.V);
    warning('Vol:Exc1','Volume of the reasonable set is %s', num2str(ctv));
    warning('Vol:Exc2','times larger than the upper set!');
    warning('Vol:Exc3','If you want to plot the reasonable set change the crit_val to a value larger than %s',num2str(rq_vol));
  end
else
ms1=min(Pip.V);
ml1=max(Pip.V);
ms2=min(Pprea.V);
ml2=max(Pprea.V);
end
if n==4
   v_prk(:,idx)=[];
   v_prn(:,idx)=[];
   v_sh(:,idx)=[];
else
   [X1,X2]=ToSimplex(v_prk);
   v_prk=[X1,X2];
   [X1,X2]=ToSimplex(v_prn);
   v_prn=[X1,X2];
   [X1,X2]=ToSimplex(v_sh);
   v_sh=[X1,X2];
end

clf;
if n==4
% Plot Reasonable Set
 if reas_sol==1
% All Solutions
  if strcmp(add_sol,'all')
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot3(v_prk(1),v_prk(2),v_prk(3));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      hn=plot3(v_prn(1),v_prn(2),v_prn(3));
      hs=plot3(v_sh(1),v_sh(2),v_sh(3));
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      if rq_vol<ctv
         eqprQ=abs(prs_vol-rs_vol)<tol;
         if eqprQ==0
            hr=Prea.plot('linewidth', 0.9);
            set(hr,'FaceLighting','phong','FaceAlpha',0.2,'FaceColor',[0.7 0 0]);
         end
      end
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
%      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Pre-Kernel
  elseif strcmp(add_sol,'prk')
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot3(v_prk(1),v_prk(2),v_prk(3));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if rq_vol<ctv
         eqprQ=abs(prs_vol-rs_vol)<tol;
         if eqprQ==0
            hr=Prea.plot('linewidth', 0.9);
            set(hr,'FaceLighting','phong','FaceAlpha',0.2,'FaceColor',[0.7 0 0]);
         end
      end
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Pre-Nucleolus
  elseif strcmp(add_sol,'prn')
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hn=plot3(v_prn(1),v_prn(2),v_prn(3));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      if rq_vol<ctv
         eqprQ=abs(prs_vol-rs_vol)<tol;
         if eqprQ==0
            hr=Prea.plot('linewidth', 0.9);
            set(hr,'FaceLighting','phong','FaceAlpha',0.2,'FaceColor',[0.7 0 0]);
         end
      end
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Shapley Value
  elseif strcmp(add_sol,'shap')
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hs=plot3(v_sh(1),v_sh(2),v_sh(3));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      if rq_vol<ctv
         eqprQ=abs(prs_vol-rs_vol)<tol;
         if eqprQ==0
            hr=Prea.plot('linewidth', 0.9);
            set(hr,'FaceLighting','phong','FaceAlpha',0.2,'FaceColor',[0.7 0 0]);
         end
      end
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% No solution
  elseif strcmp(add_sol,'none')
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if rq_vol<ctv
         eqprQ=abs(prs_vol-rs_vol)<tol;
         if eqprQ==0
            hr=Prea.plot('linewidth', 0.9);
            set(hr,'FaceLighting','phong','FaceAlpha',0.2,'FaceColor',[0.7 0 0]);
         end
      end
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Default is none
  else
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if rq_vol<ctv
         eqprQ=abs(prs_vol-rs_vol)<tol;
         if eqprQ==0
            hr=Prea.plot('linewidth', 0.9);
            set(hr,'FaceLighting','phong','FaceAlpha',0.2,'FaceColor',[0.7 0 0]);
         end
      end
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
  end
%
% No Reasonable set
 else
  if strcmp(add_sol,'all')
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot3(v_prk(1),v_prk(2),v_prk(3));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      hn=plot3(v_prn(1),v_prn(2),v_prn(3));
      hs=plot3(v_sh(1),v_sh(2),v_sh(3));
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Pre-Kernel
  elseif strcmp(add_sol,'prk')
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot3(v_prk(1),v_prk(2),v_prk(3));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Pre-Nucleolus
  elseif strcmp(add_sol,'prn')
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hn=plot3(v_prn(1),v_prn(2),v_prn(3));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Shapley Value
  elseif strcmp(add_sol,'shap')
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hs=plot3(v_sh(1),v_sh(2),v_sh(3));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% No solution
  elseif strcmp(add_sol,'none')
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Default is none
  else
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
  end
 end
%
%
% 3-Perons 
else
% Plot Core
 if reas_sol==1
% All Solutions
  if strcmp(add_sol,'all')
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot(v_prk(1),v_prk(2));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      hn=plot(v_prn(1),v_prn(2));
      hs=plot(v_sh(1),v_sh(2));
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if rq_vol<ctv
         eqprQ=abs(prs_vol-rs_vol)<tol;
         if eqprQ==0
            hr=Prea.plot('linewidth', 0.9);
            set(hr,'FaceLighting','phong','FaceAlpha',0.2,'FaceColor',[0.7 0 0]);
         end
      end
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
%      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold off;
% Pre-Kernel
  elseif strcmp(add_sol,'prk')
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot(v_prk(1),v_prk(2));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if rq_vol<ctv
         eqprQ=abs(prs_vol-rs_vol)<tol;
         if eqprQ==0
            hr=Prea.plot('linewidth', 0.9);
            set(hr,'FaceLighting','phong','FaceAlpha',0.2,'FaceColor',[0.7 0 0]);
         end
      end
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      hold off;
% Pre-Nucleolus
  elseif strcmp(add_sol,'prn')
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hn=plot(v_prn(1),v_prn(2));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      hold on
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if rq_vol<ctv
         eqprQ=abs(prs_vol-rs_vol)<tol;
         if eqprQ==0
            hr=Prea.plot('linewidth', 0.9);
            set(hr,'FaceLighting','phong','FaceAlpha',0.2,'FaceColor',[0.7 0 0]);
         end
      end
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      hold off;
% Shapley Value
  elseif strcmp(add_sol,'shap')
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hs=plot(v_sh(1),v_sh(2));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if rq_vol<ctv
         eqprQ=abs(prs_vol-rs_vol)<tol;
         if eqprQ==0
            hr=Prea.plot('linewidth', 0.9);
            set(hr,'FaceLighting','phong','FaceAlpha',0.2,'FaceColor',[0.7 0 0]);
         end
      end
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      hold off;
% No solution
  elseif strcmp(add_sol,'none')
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if rq_vol<ctv
         eqprQ=abs(prs_vol-rs_vol)<tol;
         if eqprQ==0
            hr=Prea.plot('linewidth', 0.9);
            set(hr,'FaceLighting','phong','FaceAlpha',0.2,'FaceColor',[0.7 0 0]);
         end
      end
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      hold off;
% Default is none
  else
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if rq_vol<ctv
         eqprQ=abs(prs_vol-rs_vol)<tol;
         if eqprQ==0
            hr=Prea.plot('linewidth', 0.9);
            set(hr,'FaceLighting','phong','FaceAlpha',0.2,'FaceColor',[0.7 0 0]);
         end
      end
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      hold off;
  end
% Plot No Reasonable set 
 else
% All Solutions
  if strcmp(add_sol,'all')
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot(v_prk(1),v_prk(2));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      hn=plot(v_prn(1),v_prn(2));
      hs=plot(v_sh(1),v_sh(2));
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      hold off;
% Pre-Kernel
  elseif strcmp(add_sol,'prk')
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot(v_prk(1),v_prk(2));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      hold off;
% Pre-Nucleolus
  elseif strcmp(add_sol,'prn')
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hn=plot(v_prn(1),v_prn(2));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      hold off;
% Shapley Value
  elseif strcmp(add_sol,'shap')
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hs=plot(v_sh(1),v_sh(2));
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      hold off;
% No solution
  elseif strcmp(add_sol,'none')
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      set(hpr,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0 .25 .75]);
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      hold off;
% Default is none
  else
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold on
      if core_sol==1
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.4,'FaceColor',[0 .5 1]);
      end
      hpr=Pprea.plot('linewidth', 0.9);
      if lwseQ==0
         eqilQ=abs(imp_vol-lws_vol)<tol;
         if eqilQ==0
            hlw=Plws.plot('linewidth', 1.3);
            set(hlw,'FaceLighting','phong','FaceAlpha',0.35,'FaceColor',[0.7 0.5 1]);
         end
      end      
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Kernel Catchers']);
      hold off;
  end
 end 
end
