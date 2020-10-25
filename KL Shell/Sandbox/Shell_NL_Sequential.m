                                                                
% Shell script JOEL HILMERSSON

%feature accel off

% figure(4)
%     
%     grid minor
%     PlotGeometry(a_tot, B2, Xi, Eta, res_xi, res_eta, n, m, p, q)
%     
%     view(70,20)


% Get Geometry

% Read NURBS surface from Rhino and refine as desired
clear all; close all; clc
% Add path to Rhino data (NOTE: UPDATE!)
path        = '\\file00.chalmers.se\home\joelhi\.win\My Documents\Thesis\MultiElementShellFlat.csv';
%path =        'MultiElementShellTradeMark.xls';  

% Set plot resolution properties (does not affect analysis)
res_xi      = 35;
res_eta     = 10;
willPlot    = 1;    % Set to 1 if you want to plot surface, 
                    % otherwise set to 0.

% For k-refinement use both p-refinement and h-refinement
% (i.e. set elevDeg = 1 and refineKnot = 1)
% Set p-refinement
elevDeg     = 0;
d_xi        = 1;    % Set to as many degree elevations as you want to 
                    % perform in each direction. I.e. 1 if you want to 
                    % raise p to p+1 etc.
d_eta       = 1;    % Same goes for the eta direction.
% Set h-refinement
refineKnot  = 0;
ref_xi      = 0;    % Setting refinement to n generates 2^n knot spans.
ref_eta     = 1;    % As above for eta-direction.

% Get NURBS surface
[p, q, Xi, Eta, n, m, cPts, ~, B2] = GetRhinoGeometry(path, res_xi, res_eta, willPlot, elevDeg, d_xi, d_eta, refineKnot, ref_xi, ref_eta);

% Re-arrange pt structure
pts = zeros(size(cPts,1),3);
weights = zeros(size(cPts,1),1);
ind = 1;
for j = 1 : m
    for i = 1 : n
        a = B2(:,i,j);
        pts(ind,:) = a(1:3)';
        weights(ind) = a(4);
        ind = ind + 1;
    end
end

% Store Data
% ---------------------------------------------------------------------------------------

% Store data in desired structure
nPar = cell(1,3);
nPar{1,1} = weights;
nPar{1,2} = Xi;
nPar{1,3} = Eta;

% Material properties
E           = 12e9; %[N/m2]
nu          = 0;
t           = 0.1;
ep          = [E nu t];
ngp         = 3;

% Generate connectivities

%MATLAB EXPORT FUNCTION POSSIBLY
% -------------------------------------------------------------------------------------------
[index, element, sctr, numOfEl, elRangeXi, elRangeEta] = Connectivities(p, q, Xi, Eta, m, n);
%---------------------------------------------------------------------------------------------
% Total number of control points and degrees of freedom
numCPts = n*m;
dofPerCPt = 3;
numDofs = numCPts*dofPerCPt;

% Allocate matrices
f = zeros(numDofs, 1);

% Find nodes to constrain
tol = 0.001;
              
%Get free dofs

constrDofs = [];

freeDofs = setdiff(1:numDofs,constrDofs);

bc = [constrDofs', zeros(length(constrDofs),1)];

numConstrainedDofs = length(constrDofs);

%  Lagrange constraints

c_params = [0 0;
            0 0;
            0 0;
            0 1;
            0 1;
            0 1];
        
        
c_dof = [ 1;
          2;
          3;
          1;
          2;
          3];
      
c_vals = [0;0;0;0;0;0;0;-1;0];

numLG = length(c_dof);

%Loop control

numSteps =1;
a_tot = zeros(numDofs,1);

refPts = pts;
maxIter =1;
fullNR = 1;
all_iter = zeros(numSteps,1);
loadDisp = zeros(numSteps,5);
initLoadDisp = zeros(numSteps,5);
convergence = zeros(numSteps*maxIter,1);

h = figure('units','normalized','outerposition',[0 0 1 1]);

axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testSurface.gif';

s = 1;

%Pre evaluate ShapeFunctions
% MATLAB EXPORT FUNCTION
% -------------------------------------------------------------------------------------------------------
[BasisFun,GpW,GpJ,u,v] = BasisFunctionsEvaluate(p,q,ngp,numOfEl,index,elRangeXi,elRangeEta,Xi,Eta,weights);
%--------------------------------------------------------------------------------------------------------

%Pre evaluate ShapeFunctions
% MATLAB EXPORT FUNCTION
% -------------------------------------------------------------------------------------------------------
[Bu,u_c] = applyLaGrangeBC2D(c_params,c_dof,c_vals,n,m,p,q,Xi,Eta,index,element,weights,numDofs);
% -------------------------------------------------------------------------------------------------------

% a_c = [nonzero_bc(:,1) nonzero_bc(:,2)/numSteps];
fi_current = zeros(numDofs,1);
f_e = zeros(numDofs,1);

%f_e(121:3:132) = 100;

f_e(3:3:end) =  10;

couplDofs = [1;
             2;
             3];

Bcoupl = zeros(numDofs,length(couplDofs));
         
for i = 1:length(couplDofs)
   % -------------------------------------------------------------------------------------------------------
   Bcoupl(:,i) = computeCouplingConstraints([1,1],[1 0],couplDofs(i),[n,n],[m,m],[p,p],[q,q],Xi,Xi,Eta,Eta,index,index,element,element,weights,weights,numDofs); 
   % -------------------------------------------------------------------------------------------------------
end

K = zeros(numDofs+numLG);
bc_init = bc;
uc_init = u_c;
% Loop
for step = 1:numSteps
    step
    tic
    a_c = bc_init/numSteps;
    u_c = uc_init / numSteps;
    
    if step >= 2
       f_e = zeros(numDofs,1);
    end
    
    %Get stiffness (Now linear only)
    % MATLAB EXPORT FUNCTION
    tic
% -------------------------------------------------------------------------------------------------------
    [K] = getPatchKLinear(pts,refPts,sctr,element,ngp,ep,BasisFun,GpW,GpJ);
% -------------------------------------------------------------------------------------------------------    
    toc
    f_i = fi_current;

    % MATLAB EXPORT FUNCTION
% -------------------------------------------------------------------------------------------------------
    [a_tot,delta_a,delta_r,r_LG] = solveEquations(a_tot,numDofs,freeDofs,K,f_e,f_i,Bu,-1,a_c,u_c);
% -------------------------------------------------------------------------------------------------------

  
    [pts] = updateGeometry(refPts,a_tot, 1);
    
    % Compute Internal Forces
    % MATLAB EXPORT FUNCTION
    
% -------------------------------------------------------------------------------------------------------
    [fi_current] = computePatchInternalForces(pts,refPts,sctr,element,ngp,ep,BasisFun,GpW,GpJ);
% -------------------------------------------------------------------------------------------------------
    f_i = fi_current;

    [maxDiff,fDiff] = checkEq(freeDofs,f_i,f_e,r_LG);
    
    convergence(s) = maxDiff;
    s = s + 1;
 
    % Equilibrium iteration
    % ----------------------------------------------------------------------------------------------------------    initLoadDisp(step,:) = [fi_current(33) fi_current(66) fi_current(99) fi_current(132) a_tot(132)];
    
    tol = 0.0000001*maxDiff;
    
    if tol < 0.1
        tol = 0.1;
    end
    iter = 1;
    toc
    while maxDiff > tol
 
        f_i = fi_current;
        if fullNR == 1
        
        %Stiffness matrix
    % MATLAB EXPORT FUNCTION
% -------------------------------------------------------------------------------------------------------
        [K] = getPatchKLinear(pts,refPts,sctr,element,ngp,ep,BasisFun,GpW,GpJ);
        end
    
            % MATLAB EXPORT FUNCTION
% -------------------------------------------------------------------------------------------------------    
        [a_tot,delta_a,delta_r,r_LG] = solveEquations(a_tot,numDofs,freeDofs,K,f_e,f_i,Bu,-1,[a_c(:,1) zeros(size(a_c,1),1)],zeros(size(u_c,1),1));
% -------------------------------------------------------------------------------------------------------    
    
        pts = updateGeometry(refPts,a_tot, 1);

    %Internal Forces
        % MATLAB EXPORT FUNCTION
% -------------------------------------------------------------------------------------------------------
        [fi_current] = computePatchInternalForces(pts,refPts,sctr,element,ngp,ep,BasisFun,GpW,GpJ);
% -------------------------------------------------------------------------------------------------------
        
        [maxDiff,fDiff] = checkEq(freeDofs,f_i,f_e,r_LG);
    
        convergence(s) = maxDiff;
        s = s + 1;
    
    if iter == maxIter
        print = 'Max iterations yo!'
        iter
        maxDiff
        step
        break;       
    end
    iter = iter+1;
    end
    
    % END EQ ITERATIONS
    % -------------------------------------------------------------------------------------------------------

    all_iter(step) = iter;
    
    % PLOT AND GIF
    % --------------------------------------------------------------------------------------------------------
    
    figure(2)
    clf

    hold on
    plot3(refPts(:,1),refPts(:,2),refPts(:,3),'kx')
    
    PlotGeometry(a_tot, B2, Xi, Eta, res_xi, res_eta, n, m, p, q,ep,nPar,elRangeXi,elRangeEta,element,index,1,1)
    tic
    
    includeKnots = false;
    
    t_arr = ones(11,36);
    
    [Xi_pos,Eta_pos] = getEvalPts(res_xi,res_eta);
    
    %[stressVal, Sx, Sy ,Sz] = getPatchStress(a_tot, refPts, Xi, Eta, Xi_pos, Eta_pos, n, m, p, q,ep,weights,element,index,1,-1,t_arr);
    toc
    view(70,20)
    drawnow
    hold off
    
    % GIF machen
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if step == 1 
       imwrite(imind,cm,filename,'gif','DelayTime',0.2,'Loopcount',inf); 
    elseif step == numSteps
        imwrite(imind,cm,filename,'gif','DelayTime',2,'WriteMode','append');
    else 
       imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
    end
    
    loadDisp(step,:) = [fi_current(2) fi_current(66) fi_current(99) fi_current(132) a_tot(2)];
    
    
end

%     figure(2)
%     clf
% 
%     hold on
%     plot3(refPts(:,1),refPts(:,2),refPts(:,3),'kx')
%     
%     PlotGeometry(a_tot, B2, Xi, Eta, res_xi, res_eta, n, m, p, q,ep,nPar,elRangeXi,elRangeEta,element,index)
%     view(70,20)
%     
%     hold off

    
figure(3)

plot(1:numSteps,all_iter,'k--x')

figure(4)
hold on
% plot([0; abs(loadDisp(:,5))],[0; loadDisp(:,1)],'k-x')
% plot([0; abs(loadDisp(:,5))],[0; loadDisp(:,2)],'r-x')
% plot([0; abs(loadDisp(:,5))],[0; loadDisp(:,3)],'b-x')
% plot([0; abs(loadDisp(:,5))],[0; loadDisp(:,4)],'g-x')
plot([0; abs(loadDisp(:,5))],[0; loadDisp(:,1)],'k-x')


%plot(abs(initLoadDisp(:,2)),initLoadDisp(:,1),'ro')
hold off
figure(5)

plot(1:length(convergence),convergence,'k--x')

