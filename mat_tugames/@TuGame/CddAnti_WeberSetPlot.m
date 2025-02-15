function CddAnti_WeberSetPlot(varargin)
% CDDANTI_WEBERSETPLOT plots the Weber set of game v, 
%
% Requires the Multi-Parametric Toolbox 3
% http://people.ee.ethz.ch/~mpt/3/
%
% Usage: CddAnti_WeberSetPlot(varargin)
% Define variables:
%  output:
%
%            - A plot of the anti Weber set of game v, it is identical to 
%              Weber set.
%
%  input:     At most five input arguments are admissible. Any order of
%             the input arguments below is allowed. Nevertheless, at
%             least a game v must be specified. 
%
%
%  clv        -- TuGame class object.
%  core_sol   -- An integer to draw the weber set in connection with the anti-core.
%                The anti-core will be drawn if the size is different from the 
%                Weber set, otherwise it is ommitted.
%                Permissible values are 1 (true) or 0 (false).
%                Default is 1 (true). 
%  add_sol    -- A string to invoke additional solutions into the plot.
%                Permissible solutions are:
%                'none', this is the default value.
%                'aprk', a anti-pre-kernel element will be incorporated.
%                'aprn', the anti-pre-nucleolus will be incorporated.
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
%                CddAnti_WeberSetPlot(v,'all',vw_pt)
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
%   02/14/2022        1.9.1           hme
%

narginchk(1,5);

add_sol='';
core_sol='';
tol='';
vw_pt='';

for i=1:nargin

  if ischar(varargin{i})
    if strcmp(varargin{i},'none')
       add_sol='none';
    elseif strcmp(varargin{i},'aprk')
       add_sol='aprk';
    elseif strcmp(varargin{i},'aprn')
       add_sol='aprn';
    elseif strcmp(varargin{i},'shap')
       add_sol='shap';
    elseif strcmp(varargin{i},'all')
       add_sol='all';
    elseif isempty(regexp(varargin{i},'view'))==0
       vw_pt=varargin{i};
    else
      error(['Input argument at position ' int2str(i) ' not recognized']);
    end
  elseif isa(varargin{i},'TuGame')
    clv=varargin{i};
    v=clv.tuvalues;
    N=clv.tusize;
    n=clv.tuplayers;
  else 
    if varargin{i}==1
      core_sol=varargin{i};     
    elseif varargin{i}==0
      core_sol=varargin{i};
    elseif varargin{i}>1
       error(['Input argument at position ' int2str(i) ' not recognized']);
    else
       error(['Input argument at position ' int2str(i) ' not recognized']);
    end
  end
end


if isempty(add_sol)
  add_sol='none';
end
if isempty(core_sol)
  acrQ=CddAntiCoreQ(clv);
  core_sol=acrQ;
else
  acrQ=CddAntiCoreQ(clv);
  if acrQ==0
    core_sol=acrQ;
  end
end
if isempty(tol)
  tol=10^9*eps; 
end
if isempty(vw_pt)
  vw_pt='view(120,25)';
end



Mgc=clv.AllMarginalContributions();
y1=range(Mgc);
[~,idx]=min(y1);
web_v=Mgc;
if core_sol==1
  [crv_vert,~,cr_vol,Pv]=clv.CddAntiCoreVertices(idx);
end
[impv,~,~,Pip]=clv.CddAntiImputationVertices(idx);
v_prk=clv.Anti_PreKernel();
v_prn=clv.CddAntiPrenucl();
v_sh=clv.ShapleyValue();
ms1=min(Pip.V);
ml1=max(Pip.V);
if n==4
   v_prk(:,idx)=[];
   v_prn(:,idx)=[];
   v_sh(:,idx)=[];
   web_v(:,idx)=[];
else
   [X1,X2]=ToSimplex(v_prk);
   v_prk=[X1,X2];
   [X1,X2]=ToSimplex(v_prn);
   v_prn=[X1,X2];
   [X1,X2]=ToSimplex(v_sh);
   v_sh=[X1,X2];
   [X1,X2]=ToSimplex(web_v);
   web_v=[X1,X2];
end
Pweb=Polyhedron(web_v);
web_vol=volume(Pweb);

clf;
if n==4
% Plot Core
 if core_sol==1
% All Solutions
  if strcmp(add_sol,'all')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot3(v_prk(1),v_prk(2),v_prk(3));
      grid on;
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      hn=plot3(v_prn(1),v_prn(2),v_prn(3));
      hs=plot3(v_sh(1),v_sh(2),v_sh(3));
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceLighting','phong','FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      eqwcQ=abs(cr_vol-web_vol)<tol;
      if eqwcQ==0 
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0.1 0.4 0]);
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Pre-Kernel
  elseif strcmp(add_sol,'prk')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot3(v_prk(1),v_prk(2),v_prk(3));
      grid on;
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceLighting','phong','FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      eqwcQ=abs(cr_vol-web_vol)<tol;
      if eqwcQ==0
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0.1 0.4 0]);
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Pre-Nucleolus
  elseif strcmp(add_sol,'prn')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hn=plot3(v_prn(1),v_prn(2),v_prn(3));
      grid on;
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceLighting','phong','FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      eqwcQ=abs(cr_vol-web_vol)<tol;
      if eqwcQ==0 
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0.1 0.4 0]);
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Shapley Value
  elseif strcmp(add_sol,'shap')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hs=plot3(v_sh(1),v_sh(2),v_sh(3));
      grid on;
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceLighting','phong','FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      eqwcQ=abs(cr_vol-web_vol)<tol;
      if eqwcQ==0 
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0.1 0.4 0]);
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% No solution
  elseif strcmp(add_sol,'none')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceLighting','phong','FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      eqwcQ=abs(cr_vol-web_vol)<tol;
      if eqwcQ==0 
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0.1 0.4 0]);
      end
      title(['Weber Set of the Game']);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Default is none
  else
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceLighting','phong','FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      eqwcQ=abs(cr_vol-web_vol)<tol;
      if eqwcQ==0
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0.1 0.4 0]);
      end
      title(['Weber Set of the Game']);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
  end
%
% No Core
 else
  if strcmp(add_sol,'all')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot3(v_prk(1),v_prk(2),v_prk(3));
      grid on;
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      hn=plot3(v_prn(1),v_prn(2),v_prn(3));
      hs=plot3(v_sh(1),v_sh(2),v_sh(3));
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceLighting','phong','FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Pre-Kernel
  elseif strcmp(add_sol,'prk')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot3(v_prk(1),v_prk(2),v_prk(3));
      grid on;
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceLighting','phong','FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Pre-Nucleolus
  elseif strcmp(add_sol,'prn')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hn=plot3(v_prn(1),v_prn(2),v_prn(3));
      grid on;
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceLighting','phong','FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Shapley Value
  elseif strcmp(add_sol,'shap')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hs=plot3(v_sh(1),v_sh(2),v_sh(3));
      grid on;
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceLighting','phong','FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% No solution
  elseif strcmp(add_sol,'none')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceLighting','phong','FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      title(['Weber Set of the Game']);
      eval(vw_pt);
      hold off;
      camlight('headlight');
      lighting phong;
      material shiny;
% Default is none
  else
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceLighting','phong','FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      title(['Weber Set of the Game']);
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
 if core_sol==1
% All Solutions
  if strcmp(add_sol,'all')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot(v_prk(1),v_prk(2));
      grid on;
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      hn=plot(v_prn(1),v_prn(2));
      hs=plot(v_sh(1),v_sh(2));
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      eqwcQ=abs(cr_vol-web_vol)<tol;
      if eqwcQ==0 
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceAlpha',0.3,'FaceColor',[0.1 0.4 0]);
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold off;
% Pre-Kernel
  elseif strcmp(add_sol,'prk')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot(v_prk(1),v_prk(2));
      grid on;
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      eqwcQ=abs(cr_vol-web_vol)<tol;
      if eqwcQ==0 
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceAlpha',0.3,'FaceColor',[0.1 0.4 0]);
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold off;
% Pre-Nucleolus
  elseif strcmp(add_sol,'prn')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hn=plot(v_prn(1),v_prn(2));
      grid on;
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      eqwcQ=abs(cr_vol-web_vol)<tol;
      if eqwcQ==0 
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceAlpha',0.3,'FaceColor',[0.1 0.4 0]);
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold off;
% Shapley Value
  elseif strcmp(add_sol,'shap')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hs=plot(v_sh(1),v_sh(2));
      grid on;
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      eqwcQ=abs(cr_vol-web_vol)<tol;
      if eqwcQ==0 
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceAlpha',0.3,'FaceColor',[0.1 0.4 0]);
      end
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold off;
% No solution
  elseif strcmp(add_sol,'none')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      eqwcQ=abs(cr_vol-web_vol)<tol;
      if eqwcQ==0 
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceAlpha',0.3,'FaceColor',[0.1 0.4 0]);
      end
      title(['Weber Set of the Game']);
      hold off;
% Default is none
  else
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      eqwcQ=abs(cr_vol-web_vol)<tol;
      if eqwcQ==0
         hc=Pv.plot('linewidth', 0.9);
         set(hc,'FaceAlpha',0.3,'FaceColor',[0.1 0.4 0]);
      end
      title(['Weber Set of the Game']);
      hold off;
  end
% Plot No Core
 else
% All Solutions
  if strcmp(add_sol,'all')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot(v_prk(1),v_prk(2));
      grid on;
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      hn=plot(v_prn(1),v_prn(2));
      hs=plot(v_sh(1),v_sh(2));
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold off;
% Pre-Kernel
  elseif strcmp(add_sol,'prk')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hp=plot(v_prk(1),v_prk(2));
      grid on;
      set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold off;
% Pre-Nucleolus
  elseif strcmp(add_sol,'prn')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hn=plot(v_prn(1),v_prn(2));
      grid on;
      set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold off;
% Shapley Value
  elseif strcmp(add_sol,'shap')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      hs=plot(v_sh(1),v_sh(2));
      grid on;
      set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      Pip.plot('alpha',0,'linewidth', 0.7);
      title(['Weber Set of the Game']);
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold off;
% No solution
  elseif strcmp(add_sol,'none')
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      title(['Weber Set of the Game']);
      hold off;
% Default is none
  else
      ms2=min(Pweb.V);
      ml2=max(Pweb.V);
      mrg=min(max(ml2)/4,1/2);
      sm=floor(min(ms1,ms2))-mrg;
      lr=ceil(max(ml1,ml2))+mrg;
      Pip.plot('alpha',0,'linewidth', 0.7);
      grid on;
      set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
      hold on
      h=Pweb.plot('linewidth', 1.3);
      set(h,'FaceAlpha',0.3,'FaceColor',[0.3 0.7 0]);
      title(['Weber Set of the Game']);
      hold off;
  end
 end 
end
