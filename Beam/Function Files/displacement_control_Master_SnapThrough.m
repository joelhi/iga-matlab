%% Beam Joel
clear all; clc;close all

feature accel off

% Geometry
p = 4;
tic
% if y0 is set to one all Y values will be set to zero
y0 = 0.25;
x0 = 2;
% Get geometry of a curve
[pts, Xi, Weights] = getGeometry(1,y0,x0);

%pts = [0 0;
%       1 0;
%       2 0];

%Xi = [0 0 0 1 1 1];

%Weights = [1 1 1];
   
%Material/Properties
E           = 10*10^6;
h           = 1; 
b           = 1;
A           = b*h;
I           = h^3*b/12;

Pcrit = E*I*pi^2/(10*x0)^2;

Pcrit = 0.5*Pcrit

%Number of gauss points
ngp         = 4;

% Build thickness cell

t = cell(3,1);

t_vals = [1 1];
t_knot = [0 0 1 1];
t_p = 1;

t{1}= t_vals;
t{2} = t_knot;
t{3} = t_p;

reductionFactor = 8;

% Connectivity
%------ Get Connectivity

numElems = length(unique(Xi))-1;

[elRange,elConn,index,numElems] = getConnectivity(Xi,p);

[numCpts,numDofsPerNode,numDofs,nodesPerElem,sctr] = getSctr(elConn,pts,numElems); 

% Control

zero_bc = [1;
           21;
           22];

nonzero_bc = [1 0;
              2 -1.2;
              21 0;
              22 0];

%Get free dofs

[freeDofs] = getFreeDofs(numDofs,zero_bc,nonzero_bc(:,1));

nonZeroDofs = [freeDofs;
               nonzero_bc(:,1)];
                
           
cDofs = sort([nonzero_bc(:,1)]);
sortedNonZeroDofs = sort(nonZeroDofs);
               
numConstrainedDofs = length(nonzero_bc(:,1));
%Loop control
numSteps =20;
figure(1)
hold on
a_tot = zeros(numDofs,1);
f_ext = zeros(numDofs,1);
refPts = pts;
step_force = zeros(numSteps,1);
step_force_init = zeros(numSteps,1);
a_step = zeros(numSteps,1);
a_step_step = zeros(numSteps,1);
all_feq = zeros(length(freeDofs),numSteps);
all_fif = zeros(length(freeDofs),numSteps);
plotInterval = 1;
fullNR = 0;
redSteps =0;
plot(refPts(:,1),refPts(:,2),'o--k')
tic
for step = 1:numSteps
    
    if step == numSteps
        test = 1;     
    end

    % Get tangent stiffness
    
    [K,Kmm,Kbb,~,~] = initK(numDofs);
    KNL= K;

    %Get stiffness at Tangent

    [K,~] = computeTanK(ngp,pts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,elRange,elConn,numElems,sctr,K,KNL,t);
    
    % Partition K into two parts (prescribed + nonpresctribed.
    
    %Remove fixed dofs
    
    %[K1,f1]=statcon(K,f_ext,zero_bc);
    
    
    [Kff,Kfc,Kcc] = partK(freeDofs,cDofs,K);
    
    Kcf = Kfc';
    
    %Constrained displacement vector
    a_c = nonzero_bc(:,2)/numSteps;
    
    if step < redSteps
    a_c = a_c / reductionFactor;
    end
    
    
    % Compute first step
    [fi_current] = getFinternal(ngp,pts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,numElems,elConn,elRange,sctr,t)
    % Equivalen force vector from disp
    fi_f = fi_current(freeDofs);
    f_eq = -Kfc*a_c;
    all_feq(:,step)=f_eq;
    all_fif(:,step)=fi_f;
    
    a_f = solveq(Kff,f_eq-fi_f);
    
    a_tot(freeDofs) = a_tot(freeDofs) + a_f;
    a_tot(nonzero_bc(:,1)) = a_tot(nonzero_bc(:,1))+ a_c;
    
    [defPts] = updateGeometry(refPts,a_tot, 1);
    
    [fi_current] = getFinternal(ngp,defPts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,numElems,elConn,elRange,sctr,t);
    a_step_init(step) = a_tot(2);    
    step_force_init(step) = fi_current(2);
    
    fi_f = fi_current(freeDofs);
    max_fi_f = max(abs(fi_f));
    
    % iterate for equilibrium
    
    maxIter = 100;
    pts = defPts;
    
    tol = 0.0001;
    iter = 1;
    
    while max_fi_f > tol
    
        if fullNR == 1
            
        [K,Kmm,Kbb,~,~] = initK(numDofs);
        KNL= K;
        
        [K,~] = computeTanK(ngp,pts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,elRange,elConn,numElems,sctr,K,KNL,t);
        
        [Kff,Kfc,Kcc] = partK(freeDofs,cDofs,K);
    
        Kcf = Kfc';
        end
 
    delta_a_f = -solveq(Kff,fi_f);
    a_tot(freeDofs) = a_tot(freeDofs) + delta_a_f;
    
    [defPts] = updateGeometry(refPts,a_tot, 1);
    
    [fi_current] = getFinternal(ngp,defPts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,numElems,elConn,elRange,sctr,t);
    fi_f = fi_current(freeDofs);
    
    max_fi_f = max(abs(fi_f));
    
    if iter == maxIter
        print = 'Max iterations yo!'
        iter
        max_fi_f
        step
        break;
        
    end
    iter = iter+1;
    if step == numSteps
        test = 1;     
    end
    
    pts = defPts;
    end
    iter
    if plotInterval == 2
    defPtsPlot = updateGeometry(refPts,a_tot,1);
    defCurvePlot = nrbmak(defPtsPlot',Xi);
    nrbplot(defCurvePlot,100)
    
    plotInterval = 1;
    else
    plotInterval = plotInterval + 1;    
    end
    a_step(step) =a_tot(2);
    step_force(step) = fi_current(2);
    
    
end
toc
defPtsPlot = updateGeometry(refPts,a_tot,1);
defCurvePlot = nrbmak(defPtsPlot',Xi);
nrbplot(defCurvePlot,100)

CurvePlot = nrbmak(refPts',Xi);

plot(defPtsPlot(:,1),defPtsPlot(:,2),'o--k')
nrbplot(CurvePlot,100)
nrbplot(defCurvePlot,100)
grid on
grid minor
hold off

figure(2)
hold on
%plot(abs(a_step_init),-step_force_init,'r--x') 
plot([0; abs(a_step)],[0; -step_force],'k-x')
grid on
grid minor
a_step = [0;
          a_step];
step_force = [0;
              step_force];

for i = 1:length(a_step)-2
plot([abs(a_step(i+1)) abs(a_step_init(i))], [-step_force(i+1) -step_force_init(i)],'b-')
end
hold off
% 
% % Control
% 
% zero_bc = [1;
%            2;
%            21;
%            22];
% 
% dist1=1.6-pts(2,1);
% dist2=10-1.6-pts(10,1);
%        
% nonzero_bc = [3 dist1;
%               19 dist2];
% 
% %Get free dofs
% 
% [freeDofs] = getFreeDofs(numDofs,zero_bc,nonzero_bc(:,1));
% 
% nonZeroDofs = [freeDofs;
%                nonzero_bc(:,1)];
%                 
% sortedNonZeroDofs = sort(nonZeroDofs);
%                
% numConstrainedDofs = length(nonzero_bc(:,1));
% %Loop control
% numSteps = 10;
% hold on
% f_ext = zeros(numDofs,1);
% %refPts = pts;
% step_force = zeros(numSteps,1);
% step_force_init = zeros(numSteps,1);
% plotInterval = 1;
% fullNR = 1;
% redSteps = 0;
% figure(1)
% hold on
% plot(pts(:,1),pts(:,2),'o--k')
% tic
% for step = 1:numSteps
%     
%     if step == numSteps
%         test = 1;     
%     end
% 
%     % Get tangent stiffness
%     
%     [K,Kmm,Kbb,~,~] = initK(numDofs);
%     KNL= K;
% 
%     %Get stiffness at Tangent
% 
%     [K,~] = computeTanK(ngp,pts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,elRange,elConn,numElems,sctr,K,KNL);
%     
%     % Partition K into two parts (prescribed + nonpresctribed.
%     
%     %Remove fixed dofs
%     
%     %[K1,f1]=statcon(K,f_ext,zero_bc);
%     
%     K(zero_bc,:) = [];
%     K(:,zero_bc) = [];
%     
%     %Partition
%     
%     K1 = K;
%     
%     [Kff,Kfc,Kcc] = partK(nonZeroDofs,sortedNonZeroDofs,K1,numConstrainedDofs);
%     
%     Kcf = Kfc';
%     
%     %Constrained displacement vector
%     a_c = nonzero_bc(:,2)/numSteps;
%     
%     % Compute first step
%     
%     % Equivalen force vector from disp
%     
%     f_eq = -Kfc*a_c;
%     
%     a_f = solveq(Kff,f_eq);
%     
%     a_tot(freeDofs) = a_tot(freeDofs) + a_f;
%     a_tot(nonzero_bc(:,1)) = a_tot(nonzero_bc(:,1))+ a_c;
%     
%     [defPts] = updateGeometry(refPts,a_tot, 1);
%     
%     [fi_current] = getFinternal(ngp,defPts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,numElems,elConn,elRange,sctr);
%         
%     step_force_init(step) = abs(fi_current(21));
%     
%     fi_f = fi_current(freeDofs);
%     max_fi_f = max(abs(fi_f));
%     
%     % iterate for equilibrium
%     
%     maxIter = 100;
%     pts = defPts;
%     
%     tol = 0.01;
%     iter = 1;
%     
%     while max_fi_f > tol
%  
%     if fullNR == 1
%             
%         [K,Kmm,Kbb,~,~] = initK(numDofs);
%         KNL= K;
%         
%         [K,~] = computeTanK(ngp,pts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,elRange,elConn,numElems,sctr,K,KNL);
%         
%         K(zero_bc,:) = [];
%         K(:,zero_bc) = [];
%     
%         %Partition
%     
%         K1 = K;
%     
%         [Kff,Kfc,Kcc] = partK(nonZeroDofs,sortedNonZeroDofs,K1,numConstrainedDofs);
%     
%         Kcf = Kfc';
%     end
%         
%     delta_a_f = -solveq(Kff,fi_f);
%     a_tot(freeDofs) = a_tot(freeDofs) + delta_a_f;
%     
%     [defPts] = updateGeometry(refPts,a_tot, 1);
%     
%     [fi_current] = getFinternal(ngp,defPts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,numElems,elConn,elRange,sctr);
%     fi_f = fi_current(freeDofs);
%     
%     max_fi_f = max(abs(fi_f));
%     
%     if iter == maxIter
%         print = 'Max iterations yo!'
%         iter
%         max_fi_f
%         step
%         break;
%         
%     end
%     iter = iter+1;
%     if step == numSteps
%         test = 1;     
%     end
%     end
%     
%     if plotInterval == 100
%     defPtsPlot = updateGeometry(refPts,a_tot,1);
%     defCurvePlot = nrbmak(defPtsPlot',Xi);
%     nrbplot(defCurvePlot,100)
%     
%     plotInterval = 1;
%     else
%     plotInterval = plotInterval + 1;    
%     end
%     step_force(step) = abs(fi_current(21));
%     
% end
% 
% toc
% defPtsPlot = updateGeometry(refPts,a_tot,1);
% defCurvePlot = nrbmak(defPtsPlot',Xi);
% nrbplot(defCurvePlot,100)
% 
% CurvePlot = nrbmak(refPts',Xi);
% 
% 
% nrbplot(CurvePlot,100)
% nrbplot(defCurvePlot,100)
% 
% plot(defPtsPlot(:,1),defPtsPlot(:,2),'o--r')
% 
% hold off
% 
% figure(4)
% hold on
% plot(1:length(step_force_init),step_force_init,'k--x') 
% plot(1:length(step_force),step_force,'k--o') 
% hold off

  