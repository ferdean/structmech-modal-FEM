function h=plotFEM(varargin)
%plotFEM Mesh and solution plot for 2-D PDE
%   H=PLOTFEM(PDEM) plots the mesh for a 2-D PDE defined by the PDEModel
%   object PDEM. The handles to the plotted graphics objects are returned 
%   in the optional output argument H.
%
%   H=PLOTFEM(NODES, ELEMENTS) plots the mesh for a 2-D PDE defined by the 
%   NODES and ELEMENTS matrices. NODES is a 2-by-NumNodes matrix representing 
%   nodal coordinates. ELEMENTS is a 3-by-NumElements or 6-by-NumElements 
%   matrix representing element connectivity in terms of node IDs. 
%
%   H=PLOTFEM(P,E,T) plots the mesh for a 2-D PDE defined by PET mesh
%   format. P is a 2-by-NumNodes matrix representing nodal coordinates.
%   E is a matrix representing the association between the geometry and
%   the mesh. T is a matrix representing the element connectivity in terms 
%   of node IDs. The end row of T represents the geometry face ID to which
%   the element belongs.
%
%   H=PLOTFEM(...,'Name1',Value1, 'Name2',Value2,...) plots the solution
%   to a 2-D PDE over the specified mesh or elements. 
%
%   The Name-Value pairs are used to specify the solution values and
%   plotting styles are as follows:
%
%   Name1                   Value/{Default}        Description
%   ----------------------------------------------------------------
%   XYData                  data                   - Solution data to plot 
%                                                   for example, u, abs(c*grad u)
%   XYStyle                 off | flat | {interp}  - Shaded color style
%   Contour                 {off} | on             - Plot contour lines
%   ZData                   data                   - Data for Z-height plot
%   ZStyle                  off | {continuous} | discontinuous - Plot style
%   FlowData                data                   - Data for quiver plot option
%   FlowStyle               off | {arrow}          - Quiver plot display option
%   ColorMap                name of valid colormap {'cool'} or color matrix
%   XYGrid                  {off} | on             - Option for grid plot format
%   GridParam               [tn; a2; a3] triangle index and interpolation params
%                                        from earlier call to tri2grid
%   Mesh                    {off} | on             - Display mesh edges
%   ColorBar                off | {on}             - Display the color bar
%   Title                   string {''}            - Title for the plot
%   Levels                  no of contour levels or a vector specifying levels {10}
%   NodeLabels              {off} | on    - Display node labels on the mesh.                                 
%   ElementLabels           {off} | on    - Display element labels on the mesh.                                  
%   FaceAlpha               data          - Transparency of shaded plot.
%                                           Scalar in the range [0, 1], default = 1
%   Deformation             Displacement  - For a structural analysis model, 
%                                           use this option to plot the
%                                           deformed shape. Value must be
%                                           the displacement structure,
%                                           available in
%                                           StaticStructuralResults object.
%  DeformationScaleFactor  data           - Defines the scaling factor for
%                                           plotting deformed shape. Use
%                                           this option in conjunction with
%                                           Deformation to override the
%                                           default. Default value is
%                                           computed internally based on
%                                           the geometry dimension and
%                                           magnitude of deformation.
%
%
%   Note: Plot options may only be applicable in relevant context. For example, 
%         ZStyle is ignored without ZData. Also, NodeLabels and ElementLabels 
%         are  ignored when used with ZData.
%
%       Copyright 1994-2019 The MathWorks, Inc.
%       
% =========================================================================
% Modified by Ferran de Andrés (4.2021)
% =========================================================================

% Error checks
if nargin > 0
    [varargin{:}] = convertStringsToChars(varargin{:});
end

nargs = nargin;
if nargs<1
  error(message('pde:pdeplot:nargin'));
end

if isa(varargin{1}, 'pde.EquationModel')    
    if isempty(varargin{1}.Mesh)
      error(message('pde:pdeplot:UnmeshedPdeModel'));
    end
    if rem(nargs,2)~=1
        error(message('pde:pdeplot:NoParamPairs'));
    end
elseif isa(varargin{1}, 'pde.FEMesh')    
    if rem(nargs,2)~=1
        error(message('pde:pdeplot:NoParamPairs'));
    end         
else
    narginchk(2, 33);  
    validateattributes(varargin{1},{'numeric'},{'real', 'finite', 'nonsparse', 'nonnan'}, inputname(1));
    validateattributes(varargin{2},{'numeric'},{'real', 'finite', 'nonsparse', 'nonnan'},inputname(2));
    if(nargin > 2  && isnumeric(varargin{3}))
       validateattributes(varargin{3},{'numeric'},{'real', 'finite', 'nonsparse', 'nonnan'},inputname(3)); 
       if rem(nargs,2)~=1
            error(message('pde:pdeplot:NoParamPairs'));
       end  
    else
       if rem(nargs,2)~=0
            error(message('pde:pdeplot:NoParamPairs'));
       end  
    end
end

if isa(varargin{1}, 'pde.EquationModel')   
    themsh = varargin{1}.Mesh;             
    [p,e,t] = themsh.meshToPet();
    nargs = nargs+2;
    varargin = {varargin{2:end}};
elseif isa(varargin{1}, 'pde.FEMesh')    
    themsh = varargin{1};          
    [p,e,t] = themsh.meshToPet();
    nargs = nargs+2;
    varargin = {varargin{2:end}};    
else      
    if(nargin > 2  && isnumeric(varargin{3}))
        p = varargin{1};
        e = varargin{2};
        t = varargin{3};
        validateattributes(t,{'numeric'},{'integer'},inputname(3));
        varargin = {varargin{4:end}}; 
    else       
        p = varargin{1};
        e = [];
        t = varargin{2};
        t(end+1,:) = 1;
        validateattributes(t,{'numeric'},{'integer'},inputname(2));
        varargin = {varargin{3:end}};        
    end
end

if size(p,1) ~= 2
   error(message('pde:pdeplot:NotTwoD'))
end   
     
% Default values
xydata=[];
xystyle='interp';
zdata=[];
zstyle='continuous';
flowdata=[];
flowstyle='arrow';
cmap='cool';
intool='off';
xygrid='off';
mh='off';
falpha=1;
ndlabels='off';
ellabels='off';
cbar='on';
title='';
levels=10;
tn=[]; a2=[]; a3=[];
znodedata=0;
cont='off';
deformation = [];
scaleFactor =[];
edgecolor = [0.8 0.8 0.8];
linewidth = 1;

% Recover Param/Value pairs from argument list
numVarArgs = nargs - 3;
for ii=1:2:numVarArgs
    Param = varargin{ii}; Value = varargin{ii+1};
    if ~ischar(Param)        
        error(message('pde:pdeplot:ParamNotString'))
    elseif size(Param,1)~=1
        error(message('pde:pdeplot:ParamNumRowsOrEmpty'))
    end
    switch lower(Param)
        case 'xydata'
            xydata = Value;
            if ischar(xydata)
                error(message('pde:pdeplot:xydataNotVector'))
            end
        case 'xystyle'
            xystyle=lower(Value);
            if ~ischar(xystyle)
                error(message('pde:pdeplot:GenericNotString', 'xystyle'))
            elseif ~(strcmp(xystyle,'off') || strcmp(xystyle,'flat') || ...
                    strcmp(xystyle,'interp'))
                error(message('pde:pdeplot:xystyleInvalidOption'))
            end
        case 'zdata'
            zdata = Value;
            if ischar(zdata)
                error(message('pde:pdeplot:ZdataNotVector'))
            end
        case 'zstyle'
            zstyle=lower(Value);
            if ~ischar(zstyle)
                error(message('pde:pdeplot:GenericNotString','zstyle'))
            elseif ~(strcmp(zstyle,'off') || strcmp(zstyle,'continuous') || ...
                    strcmp(zstyle,'discontinuous'))
                error(message('pde:pdeplot:ZstyleInvalidOption'))
            end
        case 'flowdata'
            flowdata = Value;
            if ischar(flowdata)
                error(message('pde:pdeplot:FlowdataNotMatrix'))
            end
        case 'flowstyle'
            flowstyle=lower(Value);
            if ~ischar(flowstyle)
                error(message('pde:pdeplot:GenericNotString','flowstyle'))
            elseif ~(strcmp(flowstyle,'off') || strcmp(flowstyle,'arrow'))
                error(message('pde:pdeplot:FlowstyleInvalidOption'))
            end
        case 'colormap'
            if ischar(Value)
                cmap = lower(Value);
            elseif size(Value,2)~=3
                error(message('pde:pdeplot:ColormapSizeOrNotString'))
            else
                cmap=Value;
            end
        case 'intool'
            intool = lower(Value);
            if ~ischar(intool)
                error(message('pde:pdeplot:GenericNotString','intool'))
            elseif ~(strcmp(intool,'off') || strcmp(intool,'on'))
                error(message('pde:pdeplot:GenericMustBeOffOn','intool'))
            end
        case 'xygrid'
            xygrid = lower(Value);
            if ~ischar(xygrid)
                error(message('pde:pdeplot:GenericNotString','xygrid'))
            elseif ~(strcmp(xygrid,'off') || strcmp(xygrid,'on'))
                error(message('pde:pdeplot:GenericMustBeOffOn','xygrid'))
            end
        case 'mesh'
            mh = lower(Value);
            if ~ischar(mh)
                error(message('pde:pdeplot:GenericNotString','mesh'))
            elseif ~(strcmp(mh,'off') || strcmp(mh,'on'))
                error(message('pde:pdeplot:GenericMustBeOffOn','mesh'))
            end
        case 'nodelabels'
            ndlabels = lower(Value);
            if ~ischar(ndlabels) || ~(strcmp(ndlabels,'off') || strcmp(ndlabels,'on'))
                error(message('pde:pdeplot:GenericMustBeOffOn','NodeLabels'))
            end
        case 'elementlabels'
            ellabels = lower(Value);
            if ~ischar(ellabels) || ~(strcmp(ellabels,'off') || strcmp(ellabels,'on'))
                error(message('pde:pdeplot:GenericMustBeOffOn','ElementLabels'))
            end
        case 'facealpha'
            falpha = Value;
            if ~isnumeric(falpha) || ~isscalar(falpha) || ~isfinite(falpha)
                error(message('pde:pdeplot:InvalidAlphaOption'))
            end
        case 'colorbar'
            cbar = lower(Value);
            if ~ischar(cbar)
                error(message('pde:pdeplot:GenericNotString','colorbar'))
            elseif ~(strcmp(cbar,'off') || strcmp(cbar,'on'))
                error(message('pde:pdeplot:GenericMustBeOffOn','colorbar'))
            end
        case 'title'
            title = Value;
            if ~ischar(title)
                error(message('pde:pdeplot:GenericNotString','title'))
            end
        case 'levels'
            levels = Value;
            if isempty(levels)
                levels=10;
            end
        case 'gridparam'
            gridparam = Value;
            if ischar(gridparam)
                error(message('pde:pdeplot:GridparamChar'))
            elseif rem(size(gridparam,1),3)
                error(message('pde:pdeplot:InvalidGridparam'))
            end
            n=size(gridparam,1)/3;
            tn=gridparam(1:n,:);
            a2=gridparam(n+1:2*n,:);
            a3=gridparam(2*n+1:3*n,:);
        case 'contour'
            cont = lower(Value);
            if ~ischar(cont)
                error(message('pde:pdeplot:GenericNotString','contour'))
            elseif ~(strcmp(cont,'off') || strcmp(cont,'on'))
                error(message('pde:pdeplot:GenericMustBeOffOn','contour'))
            end
        case 'deformation'
            deformation = Value;
        case 'deformationscalefactor'
            scaleFactor = Value;
        case 'edgecolor'
            edgecolor = Value;
        case 'linewidth'
            linewidth = Value;
        case 'facecolor'
            error(message('pde:pdeplot:FaceColorNotSupported2D'))
        otherwise
            error(message('pde:pdeplot:InvalidParam', Param))
    end
end

bbox = [min(p(1,:)) max(p(1,:));
       min(p(2,:)) max(p(2,:))];


if ~isempty(deformation) && isempty(scaleFactor)
     MaxDeformationMag = max(sqrt(deformation.ux.^2+deformation.uy.^2));
     if(MaxDeformationMag > eps(10))
        scaleFactor = min(bbox(:,2) - bbox(:,1)) /MaxDeformationMag; % Based on lowest bounding box dimension
     else
        scaleFactor = 0;
     end
end

if ~isempty(deformation)
    p(1,:) = p(1,:) +  scaleFactor*deformation.ux';
    p(2,:) = p(2,:) +  scaleFactor*deformation.uy';
end


% A few more checks
if isempty(xydata) || (strcmp(xystyle,'off') && strcmp(cont,'off'))
    plotxy=0;
else
  plotxy=1;
end
if isempty(zdata) || strcmp(zstyle,'off')
  plotz=0;
else
  plotz=1;
end
if isempty(flowdata) || strcmp(flowstyle,'off')
  plotflow=0;
else
  plotflow=1;
end

if strcmp(intool,'on')
  intool=1;
elseif strcmp(intool,'off')
  intool=0;
end

if intool
  showhidd = get(0,'ShowHiddenHandles');
  set(0,'ShowHiddenHandles','on')
  solutionpos =[0.13 0.1 0.8 0.75];
  pde_fig = findobj(get(0,'Children'),'flat','Tag','PDETool');
  ax = findobj(get(pde_fig,'Children'),'flat','Tag','PDEAxes');
  if isempty(ax)
    error(message('pde:pdeplot:IntoolNoAxes'))
  end
  hfile=findobj(get(pde_fig,'Children'),'flat','Tag','PDEFileMenu');
  flags=get(hfile,'UserData');
  if flags(2)==3
    solvemode=1;
  else
    solvemode=0;
  end
end

if ~plotxy && ~plotz && ~plotflow

  if ~isempty(t)   
    np=size(p,2);
    if size(t, 1) == 4
        T=sparse(t([1 2 3],:),t([2 3 1],:),1,np,np);  
    else
        T=sparse(t([1 4 2 5 3 6],:),t([4 2 5 3 6 1],:),1,np,np);          
    end
    if ~isempty(e)
      E=sparse(e(1,:),e(2,:),1,np,np);
      T=T>(E|E');
    end
    [I,J]=find(T|T');
    K=find(I>=J);
    I=I(K);
    J=J(K);
    X=[p(1,I); p(1,J); NaN*ones(length(I),1)'];
    Y=[p(2,I); p(2,J); NaN*ones(length(I),1)'];

    X=X(:);
    Y=Y(:);
  else
    X=[];
    Y=[];
  end

  if ~isempty(e)
    ik1=e(1,:);
    ik2=e(2,:);

    XX=[p(1,ik1)' p(1,ik2)' NaN*ones(size(ik1'))]';
    YY=[p(2,ik1)' p(2,ik2)' NaN*ones(size(ik1'))]';
    XX=XX(:);
    YY=YY(:);
  else
    XX=[];
    YY=[];
  end

  if intool
    if solvemode && ~isempty(t)
      figure(pde_fig);
      % Clean up axes:
      hndls=get(ax,'UserData');
      if ~isempty(hndls)
        delete(hndls)
        set(ax,'UserData',[]);
      end
      set(get(ax,'Children'),'Visible','off')
      pos=get(ax,'Pos');
      if any(abs(pos(1:2)-solutionpos(1:2))>100*eps)
        set(ax,'Pos',solutionpos);
      end
      set(get(ax,'Title'),...
          'Color','k',...
          'String',title,...
          'Visible','on')
    else
      set(0,'CurrentFigure',pde_fig)
    end
    set(pde_fig,'CurrentAxes',ax)
    hh1 = line(X,Y,'Color',edgecolor,'Parent',ax);
    hh2 = line(XX,YY,'Color','r','Parent',ax);
    if ~isempty(hh1)
        hh = [hh1;hh2];
    else
        hh = hh2;
    end
    if ~nargout
      % Save line handles as this is a solution plot
      set(ax,'UserData',hh)
    end
  else
    % Plot the mesh edges in blue and the free boundary in red.
    hh=plot(X,Y,'Color',edgecolor,'LineWidth',linewidth);   
    shownodelabels(p, ndlabels);
    showelemlabels(p,t, ellabels); 
    axis equal;
  end
  if nargout==1
    h=hh;
  end
  if intool
    set(0,'ShowHiddenHandles',showhidd)
  end
  return
end

ntri=size(t,2); nnode=size(p,2);
if plotxy
  xys=size(xydata);
  if xys(1)==nnode
    xynodedata=1;
    xydata=xydata(:,1);
  elseif xys(2)==ntri
    xynodedata=0;
    xydata=xydata(1,:);
  elseif xys(2)==nnode
    xydata=xydata';
    xynodedata=1;
    xydata=xydata(:,1);
  elseif xys(1)==ntri
    xydata=xydata';
    xynodedata=0;
    xydata=xydata(1,:);
  else
    error(message('pde:pdeplot:xydataLength'))
  end
end
if plotz
  zs=size(zdata);
  if zs(1)==nnode
    znodedata=1;
    zdata=zdata(:,1);
  elseif zs(2)==ntri
    znodedata=0;
    zdata=zdata(1,:);
  elseif zs(2)==nnode
    zdata=zdata';
    znodedata=1;
    zdata=zdata(:,1);
  elseif zs(1)==ntri
    zdata=zdata';
    znodedata=0;
    zdata=zdata(1,:);
  else
    error(message('pde:pdeplot:ZdataLength'))
  end
end
if plotflow
  flows=size(flowdata);
  if flows(2)==ntri
    if flows(1)<2
      error(message('pde:pdeplot:FlowdataSize'))
    else
      flowdata=flowdata(1:2,:);
    end
  elseif flows(1)==ntri
    flowdata=flowdata';
    if flows(2)<2
      error(message('pde:pdeplot:FlowdataSize'))
    else
      flowdata=flowdata(1:2,:);
    end
  elseif flows(1)==nnode
    if flows(2)<2
      error(message('pde:pdeplot:FlowdataNumCols'))
    else
      flowdata=flowdata(:,1:2);
    end
  else
    error(message('pde:pdeplot:FlowdataSizeTri'))
  end
end

if intool
  if ~plotz
    set(0,'CurrentFigure', pde_fig); 
    %
    %Clean up axes
    hndls=get(ax,'UserData');
    if ~isempty(hndls)
      delete(hndls)
      set(ax,'UserData',[]);
    end
    set(get(ax,'Children'),'Visible','off')
    pos=get(ax,'Pos');
    if any(abs(pos(1:2)-solutionpos(1:2))>100*eps)
      set(ax,'Pos',solutionpos);
    end
  else
    % Find a figure to play the movie in:
    figs=findobj(get(0,'Children'),'flat','HandleVisibility','on');
    pfig=[];
    for i=1:length(figs)
      npl=get(figs(i),'Nextplot');
      if npl(1)=='a'
        if isempty(findobj(get(figs(i),'Children')))
          pfig=figure(figs(i));
          break;
        end
      elseif npl(1)=='r'
        pfig=figure(figs(i));
        clf reset
        break;
      end
    end
    if isempty(pfig)
      figure
    end
    ax = newplot;
    view(3);
  end
else
  % not called from pdeModeler
  ax = newplot;
  hold on
  if plotz
    view(3)
  else
    view(2)
  end
end

if strcmp(xygrid,'on')
  % Use x-y grid:
  if plotxy
    if ~xynodedata
      % convert triangle data to node data
      xydata=pdeprtni(p,t,xydata);
      xynodedata=1;
    end
  end
  if plotz
    % must be nodedata if tri2grid is to be used:
    if ~znodedata
      % convert triangle data to node data
      zdata=pdeprtni(p,t,zdata);
      znodedata=1;
    end
  end

  if isempty(tn)
    % Determine xy-grid from geometry:    
    bbox = meshbbox(p, t);    
    nt=size(t,2);
    nxy=ceil(sqrt(nt/2))+1;
    x=linspace(bbox.xmin,bbox.xmax,nxy);
    y=linspace(bbox.ymin,bbox.ymax,nxy);

    if plotxy
      xydata=tri2grid(p,t,xydata,x,y);
    end
    if plotz
      zdata=tri2grid(p,t,zdata,x,y);
    end
  else
    % We have interpolation parameters
    if plotxy
      xydata=tri2grid(p,t,xydata,tn,a2,a3);
    end
    if plotz
      zdata=tri2grid(p,t,zdata,tn,a2,a3);
    end

    % Determine xy-grid from triangle geometry:
    bbox = meshbbox(p, t);       
    x=linspace(bbox.xmin,bbox.xmax,size(tn,2));
    y=linspace(bbox.ymin,bbox.ymax,size(tn,1));

  end
end

colormap(cmap)
hh=[];

% OK, now sort out all the plot cases:

% case: mesh plot (3-d only)
if ~plotxy && plotz

  if strcmp(xygrid,'on')
    % use x-y grid
    [xx,yy]=meshgrid(x,y);
    colormap([1 1 0]);
    hh=mesh(xx,yy,zdata);
  else
    % use triangular grid
    if ~znodedata
      % convert triangle data to node data
      zdata=pdeprtni(p,t,zdata);
      znodedata=1;
    end

    if size(t, 1) == 4
        it1=t(1,:);
        it2=t(2,:);
        it3=t(3,:);
        X=[p(1,it1)' p(1,it2)' p(1,it3)' p(1,it1)' NaN*ones(size(it1'))]';
        Y=[p(2,it1)' p(2,it2)' p(2,it3)' p(2,it1)' NaN*ones(size(it1'))]';
        Z=[zdata(it1) zdata(it2) zdata(it3) zdata(it1) NaN*ones(size(it1'))]';
    else    
        it1=t(1,:);
        it2=t(4,:);
        it3=t(2,:);
        it4=t(5,:);
        it5=t(3,:);
        it6=t(6,:);
        X=[p(1,it1)' p(1,it2)' p(1,it3)' p(1,it4)' p(1,it5)' p(1,it6)' p(1,it1)' NaN*ones(size(it1'))]';
        Y=[p(2,it1)' p(2,it2)' p(2,it3)' p(2,it4)' p(2,it5)' p(2,it6)' p(2,it1)' NaN*ones(size(it1'))]';
        Z=[zdata(it1) zdata(it2) zdata(it3) zdata(it4) zdata(it5) zdata(it6) zdata(it1) NaN*ones(size(it1'))]';
    end    
    X=X(:);
    Y=Y(:);
    Z=Z(:);
    hh=plot3(X,Y,Z);
  end

  % case: flat or interpolated plots
elseif (strcmp(xystyle,'flat') || strcmp(xystyle,'interp')) && plotxy

  if strcmp(xygrid,'on')

    if ~plotz
      if intool
        set(pde_fig,'CurrentAxes',ax)
      end
      zdata = zeros(size(xydata));
    end
    if strcmp(xystyle,'interp')
      if intool && ~plotz
        hold on
      end
      if strcmp(mh,'on')
        hh=surf(x,y,zdata,xydata, 'FaceAlpha',falpha);
        set(hh,'Facecolor','interp','Edgecolor','k')       
      elseif strcmp(mh,'off')
        hh=surf(x,y,zdata,xydata, 'FaceAlpha',falpha);
        set(hh,'Facecolor','interp','Edgecolor','none')
      end
      if sum(zdata(:)) == 0
            shownodelabels(p, ndlabels);
            showelemlabels(p,t, ellabels);    
      end
    elseif strcmp(xystyle,'flat')
      if intool && ~plotz
        hold on
      end
      if strcmp(mh,'on')
        hh=surf(x,y,zdata,xydata, 'FaceAlpha',falpha);        
      elseif strcmp(mh,'off')
        hh=surf(x,y,zdata,xydata, 'FaceAlpha',falpha);
        set(hh,'Facecolor','flat','Edgecolor','none')
      end
      if sum(zdata(:)) == 0
            shownodelabels(p, ndlabels);
            showelemlabels(p,t, ellabels);    
      end
    end

  else

    if plotxy
      if ~xynodedata && strcmp(xystyle,'interp')
        % convert triangle data to node data
        xydata=pdeprtni(p,t,xydata);
        xynodedata=1;
      end
    end
    if plotz
      if znodedata && strcmp(zstyle,'discontinuous')
        % convert node data to triangle data
        zdata=pdeintrp(p,t,zdata);
        znodedata=0;
      elseif ~znodedata && strcmp(zstyle,'continuous')
        % convert triangle data to node data
        zdata=pdeprtni(p,t,zdata);
        znodedata=1;
      end
    else
      zdata=zeros(size(xydata));
      znodedata=xynodedata;
    end

     if size(t, 1) == 4
        it1=t(1,:);
        it2=t(2,:);
        it3=t(3,:);
        X=[p(1,it1); p(1,it2); p(1,it3)];
        Y=[p(2,it1); p(2,it2); p(2,it3)];
        if ~znodedata
            Z=[zdata; zdata; zdata];
        else
            Z=[zdata(it1)'; zdata(it2)'; zdata(it3)'];
        end
        if ~xynodedata
            C=[xydata; xydata; xydata];
        else
            C=[xydata(it1)'; xydata(it2)'; xydata(it3)'];
        end
     else
        it1=t(1,:);
        it2=t(4,:);
        it3=t(2,:);
        it4=t(5,:);
        it5=t(3,:);
        it6=t(6,:);
        X=[p(1,it1); p(1,it2); p(1,it3); p(1,it4); p(1,it5); p(1,it6)];
        Y=[p(2,it1); p(2,it2); p(2,it3); p(2,it4); p(2,it5); p(2,it6)];
        if ~znodedata
            Z=[zdata; zdata; zdata; zdata; zdata; zdata];
        else
            Z=[zdata(it1)'; zdata(it2)'; zdata(it3)'; zdata(it4)'; zdata(it5)'; zdata(it6)'];
        end
        if ~xynodedata
            C=[xydata; xydata; xydata; xydata; xydata; xydata;];
        else
            C=[xydata(it1)'; xydata(it2)'; xydata(it3)'; xydata(it4)'; xydata(it5)'; xydata(it6)'];
        end
     end
     
    if intool && ~plotz
      set(pde_fig,'CurrentAxes',ax)
    end
    if strcmp(xystyle,'interp')
      if strcmp(mh,'on')
        hh=patch(X,Y,Z,C,'Parent',ax,'FaceAlpha',falpha);
       
      elseif strcmp(mh,'off')
        hh=patch(X,Y,Z,C,'Parent',ax, ...
            'Edgecolor','none','FaceAlpha',falpha);
      end
      if sum(Z(:)) == 0
            shownodelabels(p, ndlabels);
            showelemlabels(p,t, ellabels);    
       end
    elseif strcmp(xystyle,'flat')
      if strcmp(mh,'on')
        hh=patch(X,Y,Z,mean(C),'Parent',ax,'FaceAlpha',falpha);       
      elseif strcmp(mh,'off')
        hh=patch(X,Y,Z,mean(C),'Parent',ax,...
            'Edgecolor','none','FaceAlpha',falpha);
      end
      if sum(Z(:)) == 0
           shownodelabels(p, ndlabels);
           showelemlabels(p,t, ellabels);   
      end
    end
  end
end

% case: contour plot
if strcmp(cont,'on') && plotxy

  if ~xynodedata
    % convert triangle data to node data
    xydata=pdeprtni(p,t,xydata);
    xynodedata=1;
  end

  if ~plotz
    zdata=xydata;
  elseif ~znodedata
    % convert triangle data to node data
    zdata=pdeprtni(p,t,zdata);
  end

  xymin=min(min(xydata));
  xymax=max(max(xydata));
  zmin=min(min(zdata));
  zmax=max(max(zdata));
  if xymax==xymin, xymax=xymin+1; end
  if zmax==zmin, zmax=zmin+1; end
  if numel(levels)==1
    n=levels;
    if plotz
      lmin=(n*zmin+zmax)/(n+1);
      lmax=(zmin+n*zmax)/(n+1);
      levmin=zmin; levmax=zmax;
    else
      lmin=(n*xymin+xymax)/(n+1);
      lmax=(xymin+n*xymax)/(n+1);
      levmin=xymin; levmax=xymax;
    end
    zlmin=(n*zmin+zmax)/(n+1);
    zlmax=(zmin+n*zmax)/(n+1);
    lev=linspace(lmin,lmax,n);
    zlev=linspace(zlmin,zlmax,n);
  else
    levels=sort(levels);
    n=length(levels);
    lmin=levels(1);
    lmax=levels(n);
    zlmin=lmin;
    zlman=lmax;
    lev=levels;
    zlev=levels;
    if plotz  
      levmin=zmin; levmax=zmax;
    else
      levmin=xymin; levmax=xymax;
    end
  end

  cm=colormap;
  ncm=size(cm,1);
  icm=floor(((lev-levmin)/(levmax-levmin))*(ncm-1)+0.5)+1;
  if max(icm)>ncm || min(icm)<1
    icmindx=find(icm<=ncm & icm>=1);
    icm=icm(icmindx);
  end
  ccm=cm(icm,:);

  % Ensure that overlayed contour is drawn on top of surface plot
  if ~strcmp(xystyle,'off') && ~plotz
    set(gca,'SortMethod','childorder')
  end

  if strcmp(xygrid,'on')
    [xx,yy]=meshgrid(x,y);
    %if ~plotz
      %zdata = zeros(size(xx));
    %end
    if ~intool || plotz
      hold on      
      if strcmp(xystyle,'off')
        [~,hhc]=contour3(xx,yy,zdata,levels);
      else
        [~,hhc]=contour(xx,yy,zdata,levels);
      end
      % plot geometry boundaries:
      h1=pdeplot(p,e,[]);
    else
      set(pde_fig,'CurrentAxes',ax)
      hold on
      [~,hhc] = contour(xx,yy,zdata,levels);
      % plot geometry boundaries:
      h1=pdeplot(p,e,[],'intool','on');
    end
    set(h1,'color','k');
    d=[];
 
    if strcmp(xystyle,'off')
      hhc=[hhc(:);h1(:)];
    else
      % Black overlayed contours:
      set(hhc,'EdgeColor','k')
      hhc=[hhc(:); h1];
    end
  else
   
     if size(t, 1) == 4
       tlq = t;
    else
       tlq = [t([1, 4, 6], :), t([4, 2, 5], :), t([5, 3, 6], :), t([4, 5, 6], :)];
    end
      
    nt=size(tlq,2);
    zt=reshape(zdata(tlq(1:3,:)),3,nt);
    xyt=reshape(xydata(tlq(1:3,:)),3,nt);
    ztmax=max(zt); ztmin=min(zt);
    XX=[]; YY=[]; ZZ=[];
    for j=1:length(lev)
      jlev=zlev(j);
      it=find(ztmin<=jlev & ztmax>=jlev);
      if size(it)

        z1=zt(1,it);
        z2=zt(2,it);
        z3=zt(3,it);

        a21=zeros(1,length(it));
        itt=find(z2~=z3);       % This kludge is to avoid the warning message
        a21(itt)=(jlev-z3(itt))./(z2(itt)-z3(itt));
        itt=find(z2==z3);
        a21(itt)=NaN*ones(size(itt));
        a32=zeros(1,length(it));
        itt=find(z3~=z1);
        a32(itt)=(jlev-z1(itt))./(z3(itt)-z1(itt));
        itt=find(z3==z1);
        a32(itt)=NaN*ones(size(itt));
        a13=zeros(1,length(it));
        itt=find(z1~=z2);
        a13(itt)=(jlev-z2(itt))./(z1(itt)-z2(itt));
        itt=find(z1==z2);
        a13(itt)=NaN*ones(size(itt));

        a2=NaN*ones(2,length(it));
        a3=NaN*ones(2,length(it));
        ii=ones(1,length(it));          % 1+the number of points found so far

        itt=find(a21>=0 & a21<=1);      % On side 1
        a2(ii(itt)+2*(itt-1))=a21(itt);
        a3(ii(itt)+2*(itt-1))=1-a21(itt);
        ii(itt)=ii(itt)+ones(size(itt));
        itt=find(a32>=0 & a32<=1);      % On side 2
        a2(ii(itt)+2*(itt-1))=zeros(size(itt));
        a3(ii(itt)+2*(itt-1))=a32(itt);
        %  ii(itt)=ii(itt)+ones(size(itt));
        itt=find(a13>=0 & a13<=1);      % On side 3
        % This must be the second endpoint
        a2(2,itt)=1-a13(itt);
        a3(2,itt)=zeros(size(itt));

        X=[(1-a2(1,:)-a3(1,:)).*p(1,tlq(1,it))+ ...
                a2(1,:).*p(1,tlq(2,it))+a3(1,:).*p(1,tlq(3,it)); ...
            (1-a2(2,:)-a3(2,:)).*p(1,tlq(1,it))+ ...
                a2(2,:).*p(1,tlq(2,it))+a3(2,:).*p(1,tlq(3,it)); ...
            NaN*ones(size(it))];
        Y=[(1-a2(1,:)-a3(1,:)).*p(2,tlq(1,it))+ ...
                a2(1,:).*p(2,tlq(2,it))+a3(1,:).*p(2,tlq(3,it)); ...
            (1-a2(2,:)-a3(2,:)).*p(2,tlq(1,it))+ ...
                a2(2,:).*p(2,tlq(2,it))+a3(2,:).*p(2,tlq(3,it)); ...
            NaN*ones(size(it))];
        Z=[jlev*ones(size(it)); jlev*ones(size(it)); NaN*ones(size(it))];
        X=X(:);
        Y=Y(:);
        Z=Z(:);

        nxx=size(XX,1);
        nx=size(X,1);
        if nxx>nx
          nn=NaN*ones(nxx-nx,1);
          X=[X; nn];
          Y=[Y; nn];
          Z=[Z; nn];
        elseif nxx<nx
          nn=NaN*ones(nx-nxx,size(XX,2));
          XX=[XX; nn];
          YY=[YY; nn];
          ZZ=[ZZ; nn];
        end
        XX=[XX, X];
        YY=[YY, Y];
        ZZ=[ZZ, Z];
      end                               % size(it)
    end

    % plot geometry boundaries:
    
    if intool && ~plotz
      set(pde_fig,'CurrentAxes',ax)
    end
    hold on

    if ~plotz
      ZZ = zeros(size(XX));
    end
    for i=1:size(XX,2)
      if strcmp(xystyle,'off')
      % Colored contours:
        hndl(i)=line(XX(:,i),YY(:,i),ZZ(:,i),...
            'Parent',ax,...
            'color',ccm(i,:));
      else
      % Overlayed black contours:
        contc='k';
        hndl(i)=line(XX(:,i),YY(:,i),zeros(size(XX(:,i))),...
            'Parent',ax,...
            'color',contc);
      end
    end

    if intool && ~plotz
      h1=pdeplot(p,e,[],'intool','on');
    else
      h1=pdeplot(p,e,[]);
    end
    if strcmp(mh,'on')
        set(h1,'color','k')
    else
        set(h1,'color','none')
    end
    
    hold off

    if(exist('hndl', 'var'))
      hhc=[hndl(:);h1(:)];
    else
      hhc=h1(:);
    end
  end

  hh=[hh; hhc];  
  
end % if strcmp(cont,'on') && plotxy,


% case: add vector arrows to plot
if plotflow
  % convert triangle data to node data
  if size(flowdata,2)==ntri
    flowdata=pdeprtni(p,t,flowdata);
  end

  % Determine xy-grid from geometry:
  bbox = meshbbox(p, t);  
  % We hope 21 arrows per row looks good
  na=21;
  x=linspace(bbox.xmin,bbox.xmax,na);
  y=linspace(bbox.ymin,bbox.ymax,na);

  u=tri2grid(p,t,flowdata(:,1),x,y);
  v=tri2grid(p,t,flowdata(:,2),x,y);
  [msg,x,y]=xyzchk(x,y,u,v);
  if ~isempty(msg)
 %   if ischaruct(msg)
 %       msg = msg.message;
 %   end      
    pdeModeler('error',msg)
    if intool
      set(0,'ShowHiddenHandles',showhidd)
    end
    return
  end

  if plotxy && ~plotz
    % Setting the SortMethod of the current axes to 'childorder' will
    % ensure that the arrows are drawn on top. Not performed for
    % 3-D, though, since the disabled back to front ordering
    % destroys the appearance.
    set(gca,'SortMethod','childorder')
  end
  hold on
  oks=find(~isnan(u));
  hq=quiver(x(oks),y(oks),u(oks),v(oks),'r-');
  hh=[hh; hq];
  hold off
end

if intool && ~plotz
  set(pde_fig,'nextplot','add')
  set(ax,'UserData',hh)
end

% Finally, if there are no patches, plot an invisible patch to avoid
% problems with colorbar scaling:
if isempty(findobj(hh,'flat','Type','patch')) && plotz
  patch(0,0,0,'Parent',ax,'Visible','off')
end

% Turn on colorbar
if plotxy && strcmp(cbar,'on')
  if strcmp(xystyle,'off') && strcmp(cont,'off')
    cmax=max(max(zdata));
    cmin=min(min(zdata));
  else
    cmax=max(max(xydata));
    cmin=min(min(xydata));
  end
  if cmin~=cmax
    caxis([cmin cmax]);
  end
  if intool
    hc = findobj(get(pde_fig,'Children'),'flat','Tag','Colorbar');
    if ~isempty(hc)
      location = get(hc,'location');
    else
      location = 'EastOutside';
    end
    hc = colorbar(location,'UIContextMenu',uicontextmenu); % Disable context menu
  else
    hc=colorbar('UIContextMenu',uicontextmenu); % Disable context menu
  end
  hh=[hh; hc];
end

% turn on mouse-based 3-D rotation:
if plotz
  rotate3d on
  set(get(ax,'parent'),'currentaxes',ax)
end

% Finally, set the axes title
col = 'k';

set(get(ax,'Title'),...
    'Color',col,...
    'String',title,...
    'Visible','on')

if intool
  set(ax,'Nextplot','replace')
  set(pde_fig,'Nextplot','replace')
  if gcf~=pde_fig, set(gcf,'NextPlot','replace'), end
  set(0,'ShowHiddenHandles',showhidd)
else
  hold off
end

if nargout==1
  h=hh;
end

end
 
%
% meshbbox - Bounding box around the given elements
%
function bbox = meshbbox(p, t)
    nr = size(t,1);
    if (nr == 4 || nr == 7) 
        bbox.xmin=min(p(1,t(1:end-1,:))); 
        bbox.xmax=max(p(1,t(1:end-1,:)));
        bbox.ymin=min(p(2,t(1:end-1,:))); 
        bbox.ymax=max(p(2,t(1:end-1,:)));
    else
        bbox.xmin=min(p(1,t)); 
        bbox.xmax=max(p(1,t));
        bbox.ymin=min(p(2,t)); 
        bbox.ymax=max(p(2,t));
    end
end

function shownodelabels(p, ndloption)
   if ~strcmp(ndloption,'on')
       return;
   end
   if isempty(p)
      return;
   end
   numnds = size(p,2);
   nlabs = 'n' + string(1:numnds); nlabs = cellstr(nlabs');
   p(3,:) = 0;
   nt = matlab.graphics.primitive.world.Text('VertexData', single(p), 'String', nlabs, 'HorizontalAlignment','center', 'Clipping','on','Interpreter','none', 'Parent', gca);
   nt.Font.Name = 'Helvetica';
   nt.Font.Size = 10;
end


function showelemlabels(p,t, eloption)
   if ~strcmp(eloption,'on')
       return;
   end
   if isempty(p) || isempty(t)
      return;
   end
   numt = size(t,2);
   warnState = warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');
   tr = triangulation(t(1:3,:)', p');
   warning(warnState);
   ic = tr.incenter();
   elabs = 'e' + string(1:numt); elabs = cellstr(elabs');
   ic(:,3) = 0;
   et = matlab.graphics.primitive.world.Text('VertexData', single(ic)', 'String', elabs, 'HorizontalAlignment','center', 'Clipping','on','Interpreter','none', 'Parent', gca);
   et.Font.Name = 'Helvetica';
   et.Font.Size = 10;
end

