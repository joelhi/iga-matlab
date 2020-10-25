%% Beam Joel
clear all; clc;close all

feature accel off

% Geometry
p = 4;
tic
% if y0 is set to one all Y values will be set to zero
y0 = 0.02;
x0 = 1;
% Get geometry of a curve
[pts, Xi, Weights] = getGeometry(p,y0,x0);

%pts = [0 0;
%       1 0;
%       2 0];

%Xi = [0 0 0 1 1 1];

%Weights = [1 1 1];
   
%Material/Properties
E           = 210e9;
h           = 1; 
b           = 0.1;
t           = 1;
A           = b*h;
I           = t^3*b/12;

Pcrit = E*I*pi^2/(10*x0)^2;
%Number of gauss points
ngp         = 4;

% Build thickness cell

t = cell(3,1);

t_vals = [0.08 0.08];
t_knot = [0 0 1 1];
t_p = 1;

t{1}= t_vals;
t{2} = t_knot;
t{3} = t_p;

reductionFactor = 1;

% Connectivity
%------ Get Connectivity

numElems = length(unique(Xi))-1;

[elRange,elConn,index,numElems] = getConnectivity(Xi,p);

[numCpts,numDofsPerNode,numDofs,nodesPerElem,sctr] = getSctr(elConn,pts,numElems); 

% Control

zero_bc = [2];

nonzero_bc = [1 3.4;
              2 0;
              11 0;
              21 -3.4;
              22 0];

cDofs = [1,2,11,21,22];


freeDofs = setdiff(1:numDofs,cDofs)
%Get free dofs

[freeDofs] = getFreeDofs(numDofs,zero_bc,nonzero_bc(:,1));

nonZeroDofs = [freeDofs;
               nonzero_bc(:,1)];
                
sortedNonZeroDofs = sort(nonZeroDofs);

numConstrainedDofs = length(nonzero_bc(:,1));
%Loop control
numSteps = 300;
h = figure('units','normalized','outerposition',[0 0 1 1]);

%hold on

a_tot = zeros(numDofs,1);
f_ext = zeros(numDofs,1);
refPts = pts;
step_force = zeros(numSteps,1);
step_force_init = zeros(numSteps,1);
a_step = zeros(numSteps,1);
a_step_step = zeros(numSteps,1);
all_feq = zeros(length(freeDofs),numSteps);
all_fif = zeros(length(freeDofs),numSteps);
all_iter_2 = zeros(numSteps,1);
plotInterval = 1;
fullNR = 1;
redSteps = 0;

f_ext = zeros(numDofs,1);

f_ext(3:3:end) =0;

tic

equi = zeros(numSteps,1);
counter = 1;

axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
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
    [fi_current] = getFinternal(ngp,pts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,numElems,elConn,elRange,sctr,t);
    % Equivalen force vector from disp
    fi_f = fi_current(freeDofs);
    f_eq = -Kfc*a_c;
    all_feq(:,step)=f_eq;
    all_fif(:,step)=fi_f;
    
    a_f = solveq(Kff,f_eq-fi_f+f_ext(freeDofs));
    
    a_tot(freeDofs) = a_tot(freeDofs) + a_f;
    a_tot(nonzero_bc(:,1)) = a_tot(nonzero_bc(:,1))+ a_c;
    
    [defPts] = updateGeometry(refPts,a_tot, 1);
    
    [fi_current] = getFinternal(ngp,defPts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,numElems,elConn,elRange,sctr,t);
    a_step_init(step) = a_tot(1);    
    step_force_init(step) = fi_current(1);
    
    fi_f = fi_current(freeDofs);
    max_fi_f = max(abs(fi_f-f_ext(freeDofs)));
    
    equi(counter) = max_fi_f;
    counter = counter + 1;
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
 
    delta_a_f = -solveq(Kff,fi_f + f_ext(freeDofs));
    a_tot(freeDofs) = a_tot(freeDofs) + delta_a_f;
    
    [defPts] = updateGeometry(refPts,a_tot, 1);
    
    [fi_current] = getFinternal(ngp,defPts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,numElems,elConn,elRange,sctr,t);
    fi_f = fi_current(freeDofs);
 
    
    max_fi_f = max(abs(fi_f-f_ext(freeDofs)));
    
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
   
    equi(counter) = max_fi_f;
    counter = counter+1;
    end
    
    if plotInterval == 1
        
    
    plot(refPts(:,1),refPts(:,2),'k--o')
    hold on
    defPtsPlot = updateGeometry(refPts,a_tot,1);
    defCurvePlot = nrbmak(defPtsPlot',Xi);
    nrbplot(defCurvePlot,100)
  
    plot(defPtsPlot(:,1),defPtsPlot(:,2),'k--o')
       
    xlim([-2 12])
     ylim([-1 6])
     drawnow
    hold off
    
          frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if step == 1 
          imwrite(imind,cm,filename,'gif','DelayTime',0.1,'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
      end
    
    
    plotInterval = 1;
    else
    plotInterval = plotInterval + 1;    
    end
    
    
    a_step(step) =a_tot(1);
    step_force(step) = fi_current(1);
    iter
    all_iter_2(step) = iter;
end
toc

%defPtsPlot = updateGeometry(refPts,a_tot,1);
%defCurvePlot = nrbmak(defPtsPlot',Xi);
%nrbplot(defCurvePlot,100)

%CurvePlot = nrbmak(refPts',Xi);

%plot(defPtsPlot(:,1),defPtsPlot(:,2),'o--k')
%nrbplot(CurvePlot,100)
%nrbplot(defCurvePlot,100)

%figure(2)
%hold on
% plot(abs(a_step_init),step_force_init,'k--x') 
%plot([0; abs(a_step)],[0; step_force],'k--o')

a_step = [0;
          a_step];
step_force = [0;
              step_force];

plot(1:numSteps,all_iter_2,'x-r')
          
% for i = 1:length(a_step)-2
% plot([abs(a_step(i+1)) abs(a_step_init(i))], [step_force(i+1) step_force_init(i)],'b-')
% end

diff = step_force_init(1:length(step_force_init)-1) - step_force(2:length(step_force_init));

% Control

dist1=2.7-pts(2,1);
dist2=10-2.7-pts(10,1);
dist3=3.1-pts(3,1);
dist4=10-3.1-pts(9,1);
dist5=3.5-pts(4,1);
dist6=10-3.5-pts(8,1);
       
zero_bc = [2;
           11;
           22];

nonzero_bc = [1 0;
              21 0];

%Get free dofs

[freeDofs] = getFreeDofs(numDofs,zero_bc,nonzero_bc(:,1));

nonZeroDofs = [freeDofs;
               nonzero_bc(:,1)];
                
sortedNonZeroDofs = sort(nonZeroDofs);
               
numConstrainedDofs = length(nonzero_bc(:,1));
%Loop control
numSteps = 10;
f_ext = zeros(numDofs,1);
%refPts = pts;
step_force = zeros(numSteps,1);
step_force_init = zeros(numSteps,1);
plotInterval = 1;
fullNR = 1;
redSteps = 0;
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
    
    K(zero_bc,:) = [];
    K(:,zero_bc) = [];
    
    %Partition
    
    K1 = K;
    
    [Kff,Kfc,Kcc] = partK(nonZeroDofs,sortedNonZeroDofs,K1,numConstrainedDofs);
    
    Kcf = Kfc';
    
    %Constrained displacement vector
    a_c = nonzero_bc(:,2)/numSteps;
    
    % Compute first step
    
    % Equivalen force vector from disp
    
    f_eq = -Kfc*a_c;
    
    a_f = solveq(Kff,f_eq);
    
    a_tot(freeDofs) = a_tot(freeDofs) + a_f;
    a_tot(nonzero_bc(:,1)) = a_tot(nonzero_bc(:,1))+ a_c;
    
    [defPts] = updateGeometry(refPts,a_tot, 1);
    
    [fi_current] = getFinternal(ngp,defPts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,numElems,elConn,elRange,sctr,t);
        
    step_force_init(step) = abs(fi_current(21));
    
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
        
        K(zero_bc,:) = [];
        K(:,zero_bc) = [];
    
        %Partition
    
        K1 = K;
    
        [Kff,Kfc,Kcc] = partK(nonZeroDofs,sortedNonZeroDofs,K1,numConstrainedDofs);
    
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
    
    if plotInterval == 1
    
    plot(0,0,'xr')
        axis equal
    hold on
    plot(3,0,'xr')
    plot(7,0,'xr')
    plot(10,0,'xr')
    
    
%     plot(refPts(:,1),refPts(:,2),'k--o')
%     axis equal
%     xlim([-1 11])
%     ylim([-1 5])
%     grid on 
%     grid minor
%     drawnow
    
    
    defPtsPlot = updateGeometry(refPts,a_tot,1);
    defCurvePlot = nrbmak(defPtsPlot',Xi);
    nrbplot(defCurvePlot,20)


    
    plot(defPtsPlot(:,1),defPtsPlot(:,2),'k--o')
    
    [kappa,epsilon,A3,a3,A_pos,a_pos] = computeStrain(0:0.01:1,numCpts,p,Xi,Weights,refPts,pts);
    scaleFact = 1;
    col = kappa / max(abs(kappa));
    tempXi = 0;
    for i = 1:length(kappa)
        
        a_upd = [0 0];
        a_upd = a_pos(i,:)+scaleFact*kappa(i)*a3(i,:);
        tempColor = col(i);
        
        if tempColor <0.01
            tempColor = 0;
        end
        
        plot([a_pos(i,1) a_upd(1,1)],[a_pos(i,2) a_upd(1,2)],'color',[tempColor 0.5 0.5])
        tempXi = tempXi + 0.05;

    
    end
        xlim([-1 11])
        ylim([-1 5])
        grid on 
        grid minor
        drawnow

    hold off


    hold off
    
              frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 


    
    
    plotInterval = 1;
    else
    plotInterval = plotInterval + 1;    
    end
    step_force(step) = abs(fi_current(21));
    
end

toc
%defPtsPlot = updateGeometry(refPts,a_tot,1);
%defCurvePlot = nrbmak(defPtsPlot',Xi);
%nrbplot(defCurvePlot,100)

%CurvePlot = nrbmak(refPts',Xi);


%nrbplot(CurvePlot,100)
%nrbplot(defCurvePlot,100)

%plot(defPtsPlot(:,1),defPtsPlot(:,2),'o--r')

hold off

%figure(4)

%plot(1:length(step_force_init),-step_force_init,'k--x') 
%plot(1:length(step_force),-step_force,'k--o') 


% Plot bending Moment/Normal Force on curved geometry


Xi_coordinates = linspace(0,1,0.1);

[kappa_xi,epsilon_x] = computeStrain(Xi_coordinates,numCpts,p,Xi,Weights,refPts,defPts);


  