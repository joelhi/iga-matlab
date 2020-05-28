
% Shell script JOEL HILMERSSON

%feature accel off

% figure(4)
%     
%     grid minor
%     PlotGeometry(a_tot, B2, Xi, Eta, res_xi, res_eta, n, m, p, q)
%     
%     view(70,20)


% Get Geometry

clear all; close all; clc
% Add path to Rhino data (NOTE: UPDATE!)
%path        = '\\file00.chalmers.se\home\joelhi\.win\My Documents\Thesis\MultiElementShellTradeMark.xls';
path =        'MultiElementShellTradeMark.xls';  

% Set plot resolution properties (does not affect analysis)
res_xi      = 20;
res_eta     = 20;
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
E           = 1e7;
nu          = 0;
t           = 0.05;
ep          = [E nu t];
ngp         = 3;

% Generate connectivities
[index, element, sctr, numOfEl, elRangeXi, elRangeEta] = Connectivities(p, q, Xi, Eta, m, n);

% Total number of control points and degrees of freedom
numCPts = n*m;
dofPerCPt = 3;
numDofs = numCPts*dofPerCPt;

% Allocate matrices
f = zeros(numDofs, 1);

% Find nodes to constrain
tol = 0.001;

j = 1;
k = 1;

fixDofs = zeros(3*1,1);
moveDofs = zeros(3*1*1,2);

moveDist = -0.1;

for i = 1:length(pts)

    %X values
    
%     if abs(pts(i,1)) < tol
%         fixDofs(j:j) = [(3*i)-2];
%         j = j+1;
%     end
%     
    if abs(pts(i,2)) < tol
        fixDofs(j:j+1) = [(3*i)-2 (3*i)];
        j = j+2;
    end
    
    if abs(pts(i,2)) < tol && abs(pts(i,1)) < tol
        fixDofs(j:j+1) = [(3*i)-2 (3*i)];
        j = j+2;
    end
    
    if abs(pts(i,2)) < tol
        moveDofs(k,1) = (3*i)-1;
        moveDofs(k,2) = -moveDist + (moveDist*pts(i,1))/8;
        k = k+1;
    end
    
    %X values
    if abs(pts(i,2)-max(pts(:,2))) < tol
        fixDofs(j:j+1) = [(3*i)-2 (3*i)];
        j = j+2;
    end
    
    if abs(pts(i,2)-max(pts(:,2))) < tol && abs(pts(i,1)) < tol
        fixDofs(j:j+1) = [(3*i)-2 (3*i)];
        j = j+2;
    end
%     
%     if abs(pts(i,2)-0.5*max(pts(:,2))) < tol 
%         fixDofs(j:j+1) = [(3*i)-2 (3*i)-1];
%         j = j+2;
%     end
    
    if abs(pts(i,2)-max(pts(:,2))) < tol
        moveDofs(k,1) = (3*i)-1;
        moveDofs(k,2) = moveDist - (moveDist*pts(i,1))/8;
        k = k+1;
    end
     
end

% Boundary Conditions


[uniqueMoveDofs,index1,index2] = unique(moveDofs(:,1));

nonzero_bc = moveDofs(index1,:); 
              
%Get free dofs

freeDofs = setdiff(1:numDofs,nonzero_bc(:,1))'; 

nonZeroDofs = [freeDofs;
               nonzero_bc(:,1)];
              
sortedNonZeroDofs = sort(nonZeroDofs);

numConstrainedDofs = length(nonzero_bc(:,1));
%Loop control

numSteps =20;
a_tot = zeros(numDofs,1);

refPts = pts;
maxIter = 100;
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
[BasisFun,GpW,GpJ] = BasisFunctionsEvaluate(p,q,ngp,numOfEl,index,elRangeXi,elRangeEta,Xi,Eta,weights);
%--------------------------------------------------------------------------------------------------------

a_c = [nonzero_bc(:,1) nonzero_bc(:,2)/numSteps];
fi_current = zeros(numDofs,1);
f_e = zeros(numDofs,1);

f_e(3:3:numDofs) = 0.5;
% Loop
for step = 1:numSteps
    step
    a_c = [nonzero_bc(:,1) (nonzero_bc(:,2)/numSteps)];
    
    if step == 2
       f_e = zeros(numDofs,1);
        
    end
    
    tic
    
    %Get stiffness (Now linear only)
    % MATLAB EXPORT FUNCTION
% -------------------------------------------------------------------------------------------------------
    [K] = shellKL_IGA(pts,refPts,sctr,element,ngp,ep,BasisFun,GpW,GpJ);
% -------------------------------------------------------------------------------------------------------    
    % APPLY LAGRANGE CONSTRAINTS
    
    % fun Apply2DLagrange(x,y,z)
   
    % Compute first step

    f_i = fi_current;
    
    f = f_e - f_i;
    % MATLAB EXPORT FUNCTION
% -------------------------------------------------------------------------------------------------------
    [delta_a,r] = solveq(K,f,a_c);
    
    a_tot = a_tot + delta_a;
  
    [pts] = updateGeometry(refPts,a_tot, 1);
    
    % Compute Internal Forces
    % MATLAB EXPORT FUNCTION
% -------------------------------------------------------------------------------------------------------
    [fi_current] = shellKL_Fint(pts,refPts,sctr,element,ngp,ep,BasisFun,GpW,GpJ);
    
    fe_f = f_e(freeDofs);
    fi_f = fi_current(freeDofs);
    max_fi_f = max(abs(fe_f-fi_f));
    convergence(s) = max_fi_f;
    s = s + 1;
 
    % Equilibrium iteration
    % ----------------------------------------------------------------------------------------------------------
    
    initLoadDisp(step,:) = [fi_current(33) fi_current(66) fi_current(99) fi_current(132) a_tot(132)];
    
    tol = 0.01;
    iter = 1;

    while max_fi_f > tol
 
        if fullNR == 1
        
        %Stiffness matrix
    % MATLAB EXPORT FUNCTION
% -------------------------------------------------------------------------------------------------------
        [K] = shellKL_IGA(pts,refPts,sctr,element,ngp,ep,BasisFun,GpW,GpJ);
    % MATLAB EXPORT FUNCTION
% -------------------------------------------------------------------------------------------------------    
        % Rewrite this 
        [Kff,Kfc,Kcc] = partK(nonZeroDofs,sortedNonZeroDofs,K,numConstrainedDofs);
    
        end
        
    % fun Apply2DLagrange(x,y,z)
    fe_f = f_e(freeDofs);
    f_f = -(fe_f-fi_f);
    
     % MATLAB EXPORT FUNCTION
% -------------------------------------------------------------------------------------------------------
    delta_a_f = -solveq(Kff,f_f);
% -------------------------------------------------------------------------------------------------------    

    a_tot(freeDofs) = a_tot(freeDofs) + delta_a_f;
    
    pts = updateGeometry(refPts,a_tot, 1);

    %Internal Forces
        % MATLAB EXPORT FUNCTION
% -------------------------------------------------------------------------------------------------------
    [fi_current] = shellKL_Fint(pts,refPts,sctr,element,ngp,ep,BasisFun,GpW,GpJ);
% -------------------------------------------------------------------------------------------------------
    fi_f = fi_current(freeDofs);
    
    max_fi_f = max(abs(fi_f-fe_f));
    convergence(s) = max_fi_f;
    s = s + 1;
    
    if iter == maxIter
        print = 'Max iterations yo!'
        iter
        max_fi_f
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
    toc
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


