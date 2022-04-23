% -------------------------------------------------------------------------
% File: main.m
%
% Author: Ferran de Andres. deandresvertferran(at)gmail(dot)com
%
% Discription: Linear FEM solver for structural mechanics applications.
% Specifically, the current version solves a modal problem on a given
% externally generated mesh (giving as an output the natural frequencies
% and characeristic modes of the solid). Allows for multimaterial single
% body domains. Structural damping is implemented.
% 
% Environment: MATLAB
%
% MIT License
% 
% Copyright (c) 2022 Ferran de Andrés Vert
% -------------------------------------------------------------------------

clear all; close all; clc;

addpath lib
addpath tests
addpath plots

run materialdata.m

load('meshPlate.mat');
% load('meshWheel.mat');

fixednodes = [];     % Homogeneous dirichlet BC are set in these nodes

plotMesh   = true;

if plotMesh   
    
    figure(1)
    plotFEM3D(mesh.nodes.', mesh.top(1:4,:))
    hold on

    set(gca,'fontname','times')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    daspect([1 1 1]);
    
    xlabel('Width [m]')
    ylabel('Height [m]')
    
    title({['Mesh plot'],...
        ['\fontsize{12}\color{gray} Number of nodes: ',num2str(size(mesh.nodes,1)), '. Number of elements: ',num2str(size(mesh.top,2))]},'fontsize',14)
    
    hold off
end

%% Mass and stiffness matrices - Fast assembly algorithm

verbose = false;         % Note: increases computation time
[M, K, DOF, nDOF, DOFfree] = fastMatrices(mesh, E, v, rho, eta, fixednodes, verbose);

% spy(M)

%% Natural frequencies and modes 

tr             = 10;    % Number of natural frequencies extracted
[wnHz, phi, ~] = modal(K, M, tr);

%% Natural mode representation
%%% Representation parameters:
mode        = 8;           % Mode number
escaler     = 1e-2;        % [1e-4, 1]. If representation does not make sense,
                           % try reducing the scaler.

original    = true;        % Plot undeformed mesh 

modefig = figure(2);
modedraw3D(mode, mesh, DOF, DOFfree, abs(wnHz), real(phi), escaler, original)

% saveas(modefig,['mode_',num2str(mode)],'png')

%% Natural mode animation

numberLoops = 1;
modeAnimate3D(mode, mesh, DOF, DOFfree, wnHz, real(phi), escaler, numberLoops, original)

%% Mode clasification

modeClasification(mesh.nodes, fixednodes, phi, wnHz)

%% Nodal transformation

Mq  = phi.'*M*phi;
Kq  = phi.'*K*phi;

% Chop
Mq(Mq<1e-5) = 0;      
Mq(abs((Mq-1))<1e-5) = 1;
Kq(Kq<1) = 0; 

%% FRF

verbose  = true;
[~, ~]   = FRF(Mq, Kq, DOF, nDOF, phi, [10, 10],'axial', 'radial', verbose);

% save('HraW.mat','HraW')







