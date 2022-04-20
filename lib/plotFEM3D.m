function h = plotFEM3D( varargin )
%PDEPLOT3D Plot mesh and solution for 3-D PDE
%   H=PDEPLOT3D(PDEM) plots the mesh for a 3-D PDE defined by the PDEModel
%   object PDEM. The handles to the plotted graphics objects are returned 
%   in the optional output argument H.
%
%   H=PDEPLOT3D(NODES, ELEMENTS) plots the mesh for a 2-D PDE defined by the 
%   NODES and ELEMENTS matrices. NODES is a 3-by-NumNodes matrix representing 
%   nodal coordinates. ELEMENTS is a 4-by-NumElements or 10-by-NumElements 
%   matrix representing element connectivity in terms of node IDs. 
%
%   H=PDEPLOT3D(...,'Name1',Value1, 'Name2',Value2,...) plots the solution
%   to a 3-D PDE over the specified mesh or elements. 
%
%   The Name-Value pairs are used to specify the solution values and
%   plotting styles are as follows:
%
%   Name        Value/{Default}   Description
%   ----------------------------------------------------------------
%   ColormapData            data          - Colormap plot.
%                                           ColormapData is a vector of nodal 
%                                           solution values.
%   Mesh                    {off} | on    - Display mesh edges. This option
%                                           is used in conjunction with
%                                           ColormapData. It is ignored if
%                                           ColormapData it not specified.
%   NodeLabels              {off} | on    - Display node labels on the 
%                                           mesh boundary. Use in
%                                           conjunction with FaceAlpha.
%   ElementLabels           {off} | on    - Display element labels on the 
%                                           mesh boundary. Use in
%                                           conjunction with FaceAlpha.
%   FaceAlpha               data          - Transparency of faces.
%                                           Scalar in the range [0, 1],
%                                           default = 1
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
%
%   See also: CREATEPDE, pde.PDEMODEL, PDECONT, PDEGPLOT, pde.FEMesh, PDESURF
%
%   Copyright 2014-2019 The MathWorks, Inc.
%
% =========================================================================
% Modified by Ferran de Andrés (7.2021)
% =========================================================================

if nargin > 0
    [varargin{:}] = convertStringsToChars(varargin{:});
end

nargs = nargin;
if nargs < 1
  error(message('pde:pdeplot:nargin'))
end

thepde = [];
nodeAndElementDef = false;
if isa(varargin{1}, 'pde.EquationModel') 
  if rem(nargs,2)~=1
        error(message('pde:pdeplot:NoParamPairs'))
  end
  thepde = varargin{1};
  if ~isa(thepde.Mesh, 'pde.FEMesh')
      error(message('pde:pdeplot:UnmeshedPdeModel'))
  elseif thepde.IsTwoD
      error(message('pde:pdeplot:NotThreeD'))
  else
      themsh = thepde.Mesh;
      [p,~,t] = themsh.meshToPet();    
      varargin = {varargin{2:end}};      
  end
elseif isa(varargin{1}, 'pde.FEMesh') 
  themsh = varargin{1};
  if rem(nargs,2)~=1
      error(message('pde:pdeplot:NoParamPairs'))
  end
  if size(themsh.Nodes,1) ~= 3
      error(message('pde:pdeplot:NotThreeD'))
  else      
      [p,~,t] = themsh.meshToPet();    
      varargin = {varargin{2:end}};      
  end  
else
    if nargs < 2
        error(message('pde:pdeplot:InvalidArgs'))
    elseif rem(nargs,2)~=0
        error(message('pde:pdeplot:NoParamPairs'))
    end
    p = varargin{1};
    t = varargin{2};
    nodeAndElementDef = true;
    validateattributes(varargin{1},{'numeric'},{'real', 'finite', 'nonsparse', 'nonnan'},'pdeplot3D',inputname(1));
    validateattributes(varargin{2},{'numeric'},{'real', 'integer'},'pdeplot3D',inputname(2));
    varargin = {varargin{3:end}};    
end

if size(p,1) ~= 3
   error(message('pde:pdeplot:NotThreeD'))
end   

parser = inputParser;
addParameter(parser,'colormapdata', [], @isValidU);
addParameter(parser,'flowdata', [], @isnumeric);
addParameter(parser,'Deformation', []);
addParameter(parser,'DeformationScaleFactor', [], @isnumeric);
addParameter(parser,'NodeLabels', 'off', @isValidNdLabelOption);
addParameter(parser,'ElementLabels', 'off', @isValidElLabelOption);
addParameter(parser,'Mesh', 'on', @isValidMeshDispOption);
addParameter(parser,'FaceAlpha', 0.8, @isnumeric);
addParameter(parser,'EdgeColor', [0.6 0.6 .6]);
addParameter(parser,'FaceColor', [0.9 0.9 0.9]);
parse(parser,varargin{:});
colormapdata = parser.Results.colormapdata;
flowdata = parser.Results.flowdata;
faceAlpha = parser.Results.FaceAlpha;
showMesh = parser.Results.Mesh;
showNodeLabels = parser.Results.NodeLabels;
showElemLabels = parser.Results.ElementLabels;
deformation = parser.Results.Deformation;
scaleFactor = parser.Results.DeformationScaleFactor;
faceColor = parser.Results.FaceColor;
edgeColor = parser.Results.EdgeColor;

if(strcmpi(showElemLabels, 'on') && nodeAndElementDef && size(p,1) == 3) 
    error(message('pde:pdeplot:ElemLabelsNeedsCompleteMesh'));
end

if(size(t,1) == 5 || size(t,1) == 11)
    t(end,:)=[];
end   

numElemNodes = size(t,1);
if(numElemNodes ~= 4 && numElemNodes ~= 10)
  error(message('pde:pdeModel:invalidT', numElemNodes));
end
numNodes = size(p,2);
if ~isempty(colormapdata)
  if(length(colormapdata) ~= numNodes)
    error(message('pde:pdeplot:colormapdataLength'));
  end
end

ha = newplot;
hf = get(ha,'Parent');
set(hf,'Color','white');
set(ha,'ClippingStyle','rectangle');
hold on;

bbox = [min(p(1,:)) max(p(1,:));
        min(p(2,:)) max(p(2,:));
        min(p(3,:)) max(p(3,:))];
    
if ~isempty(deformation) && isa(thepde,'pde.StructuralModel')
    if isempty(scaleFactor)
        MaxDeformationMag = max(sqrt(deformation.ux.^2+deformation.uy.^2+deformation.uz.^2));
        if(MaxDeformationMag > eps())
            scaleFactor = min(bbox(:,2) - bbox(:,1)) /MaxDeformationMag; % Based on lowest bounding box dimension
        else
            scaleFactor = 0;
        end
    end
    p(1,:) = p(1,:) +  scaleFactor*deformation.ux';
    p(2,:) = p(2,:) +  scaleFactor*deformation.uy';
    p(3,:) = p(3,:) +  scaleFactor*deformation.uz';
end

bbox = [min(p(1,:)) max(p(1,:));
        min(p(2,:)) max(p(2,:));
        min(p(3,:)) max(p(3,:))];
 
plotAxis3D(bbox)


if ~isempty(colormapdata) 
  [ltri] = tetBoundaryFacets(p,t);  
  if(numElemNodes == 10)
      ltri = splitQuadraticTri(ltri);    
  end
  h1=colorbar(ha); ht = hgtransform(ha); colormap(ha,'parula');
  if min(colormapdata) ~= max(colormapdata)
      caxis(ha, [min(colormapdata) max(colormapdata)]);
  end
  
  h2=patch('Faces',ltri, 'Vertices', p', 'FaceVertexCData', colormapdata(:), ...
    'AmbientStrength', .75,  ...
    'EdgeColor', 'none', 'FaceColor', 'interp', 'parent',ht,'FaceAlpha',faceAlpha, 'Clipping','off');
  if strcmpi(showMesh, 'on')   
    [ltri, lp] = tetBoundaryFacets(p,t(1:4,:));  
    tr = triangulation(ltri, lp);
    tre = (tr.edges())';
    x = lp(:,1);
    y = lp(:,2);
    z = lp(:,3);
    xt=x(tre); xt(3,:) = NaN; 
    yt=y(tre); yt(3,:) = NaN;
    zt=z(tre); zt(3,:) = NaN;
    plot3(xt(:),yt(:),zt(:),'-k');    
  end
  if(strcmpi(showNodeLabels, 'on'))
    plotNodeLabels(p, t);
  end
  if(strcmpi(showElemLabels, 'on') && ~isempty(thepde))      
      plotElementLabels(thepde.Mesh);
  end     

  if nargout==1
    h = [h1 h2];
  end
elseif ~isempty(flowdata)
    if(length(flowdata(:)) ~= 3*size(p,2))
       error('Length of flowdata must be 3*number of points');
    end
    if(isvector(flowdata)==1)
      flowdata = reshape(flowdata, numNodes, 3);
    end
    quiver3(p(1,:)', p(2,:)', p(3,:)', flowdata(:,1), flowdata(:,2), flowdata(:,3));
else
  %Plot the boundary triangulation:
  [ltri, lp] = tetBoundaryFacets(p,t(1:4,:));  
  set(gcf, 'renderer', 'opengl');
  hh=trisurf(ltri, lp(:,1),lp(:,2),lp(:,3), ...
    'FaceColor',faceColor, 'EdgeColor',edgeColor,'FaceAlpha', faceAlpha, 'Clipping','off');
  if nargout==1
    h = hh;
  end    
  if(strcmpi(showNodeLabels, 'on'))
    plotNodeLabels(p,t);
  end
  if(strcmpi(showElemLabels, 'on'))
     plotElementLabels(themsh);
  end     
end
axis equal; axis tight; axis off; view(30,30);
hold off;
end


function plotNodeLabels(p,t)
    tri = tetBoundaryFacets(p,t);      
    vispts = unique(tri(:));   
    warnState = warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');    
    tr = triangulation(tri(:,1:3),p');
    if size(tri,2) >= 6
        tr2 = triangulation(tri(:,4:6),p');
    end
    warning(warnState);
    [~, ric] = tr.incenter();  
    ric = mean(ric);
    vn = tr.vertexNormal(vispts);
    if size(tri,2) >= 6
        vn = vn + tr2.vertexNormal(vispts);
    end
    vn = (vn.*(0.15*ric))';
    px = p(1,vispts)+ vn(1,:);
    py = p(2,vispts)+ vn(2,:);
    pz = p(3,vispts)+ vn(3,:);
    labels = 'n' + string(1:max(vispts)); labels = cellstr(labels');
    labels = labels(vispts);
    pos = single([px ; py; pz]);
    nt = matlab.graphics.primitive.world.Text('VertexData', pos, 'String', labels, 'HorizontalAlignment','center', 'Clipping','on','Interpreter','none', 'Parent', gca);
    nt.Font.Name = 'Helvetica';
    nt.Font.Size=8;
end

function plotElementLabels(msh)
    ma = msh.MeshAssociation;
    fa = ma.FaceAssociativity;
    facets = cell2mat(fa');
    n = msh.Nodes';
    e = msh.Elements(1:4,:)';
    lfc = ma.FaceCodes(:,1:3); % Linear Face codes
    facets = facets';
    i = sub2ind(size(e), facets(:,1), lfc(facets(:,2),1));
    j = sub2ind(size(e), facets(:,1), lfc(facets(:,2),2));
    k = sub2ind(size(e), facets(:,1), lfc(facets(:,2),3));
    ei = e(i);
    ej = e(j);
    ek = e(k);
    nf = [ei(:), ej(:), ek(:)];
    warnState = warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');    
    tr = triangulation(nf,n);
    warning(warnState);
    [ic, ric] = tr.incenter();
    fn = faceNormal(tr);
    ic = ic + (fn.*ric)*0.15;
    labels = 'e' + string(1:max(facets(:,1))); labels = cellstr(labels');
    labels = labels(facets(:,1));
    pos = single(ic');
    et = matlab.graphics.primitive.world.Text('VertexData', pos, 'String', labels, 'HorizontalAlignment','center', 'Clipping','on','Interpreter','none', 'Parent', gca);
    et.Font.Name = 'Helvetica';
    et.Font.Size=8;
end


function ok = isValidU(udata)
    if ~isreal(udata) || issparse(udata) || any(~isfinite(udata))
        error(message('pde:pdeplot:InvalidColormapdata'));   
    end     
    ok = true;
end

function ok = isValidMeshDispOption(dispopt)
    if ~ischar(dispopt)
       error(message('pde:pdeplot:GenericMustBeOffOn','mesh'));   
    end
    validatestring(dispopt,{'on', 'off'}) ;
    ok = true;
end

function ok = isValidNdLabelOption(labelopt)
    if ~ischar(labelopt)
       error(message('pde:pdeplot:GenericMustBeOffOn', 'NodeLabels'));   
    end
    validatestring(labelopt,{'on', 'off'}) ;
    ok = true;
end

function ok = isValidElLabelOption(labelopt)
    if ~ischar(labelopt)
       error(message('pde:pdeplot:GenericMustBeOffOn', 'ElementLabels'));   
    end
    validatestring(labelopt,{'on', 'off'}) ;
    ok = true;
end

function t4=splitQuadraticTri(t)
t4Nodes = [1 4 6; 4 5 6; 4 2 5; 6 5 3];
t4 = [t(:,t4Nodes(1,:)); t(:,t4Nodes(2,:)); t(:,t4Nodes(3,:)); t(:,t4Nodes(4,:))];
end

