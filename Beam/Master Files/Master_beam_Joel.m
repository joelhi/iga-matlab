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

%Material/Properties
E           = 12*10^9;
h           = 1; 
b           = 1;
t           = 0.08;
A           = b*h;
I           = t^3*b/12;

Pcrit = E*I*pi^2/(10*x0)^2;

%Number of gauss points
ngp         = 4;

% Connectivity
%------ Get Connectivity

[elRange,elConn,index,numElems] = getConnectivity(Xi,p);

[numCpts,numDofsPerNode,numDofs,nodesPerElem,sctr] = getSctr(elConn,pts,numElems); 

% Control

% Build thickness cell

t = cell(3,1);

t_vals = [0.08 0.08];
t_knot = [0 0 1 1];
t_p = 1;

t{1}= t_vals;
t{2} = t_knot;
t{3} = t_p;

bc = [2 0;
      11 0;
      22 0];

% External load
  
% Apply lagrange constraints on parameter
Xi_b = [0.8] ;
Xi_dof = [1];
Xi_b_val = [-0.001];

%find nodes affected of lagrange constraint

knotSpanIndex = FindSpan(numCpts-1,p,1,Xi) + 1
element = knotSpanIndex - p;
activeCpts = element:element+p;
c_dofs = zeros(p+1,1);

for i = 1:length(activeCpts)
    
    if Xi_dof(1) == 1
    c_dofs(i) = (2*activeCpts(i))-1;
    end
    if Xi_dof(1) == 2
    c_dofs(i) = (2*activeCpts(i));
    end
    
end
%If 1 lagrange constraints are included
applyLaGrange =0;

%[Xi,pts] = RefineKnotVectCurve(length(pts)-1,p,Xi,pts,X,length(X)-1);

%plot(Qw(:,1),Qw(:,2),'k--x')
%axis equal


%------------------------------------------------------------------------------------------------
tic
            
if applyLaGrange == 1
f_step = zeros(numDofs+length(Xi_b),1);
else
f_step = zeros(numDofs,1);
end


% Loop over Elements
numSteps = 2000;
fMax = -116000;
fStepVAL = fMax/numSteps;
fdof = [21];
    
f_step(fdof) = fStepVAL;
f_step(1) = -fStepVAL;


all_a = zeros(numDofs,numSteps);

%Initial Parameters
if applyLaGrange == 1
fi = zeros(numDofs+length(Xi_b),1);
else
fi = zeros(numDofs,1);
end
f_ext = zeros(numDofs,1);
a_tot = zeros(numDofs,1);

initBC = bc;

all_iter = zeros(numSteps,1);
all_detK = zeros(1,numSteps);
value = 1;
loadVal = zeros(numSteps,1);
loadValInt =zeros(numSteps,1);
refPts = pts;
iter_one = zeros(numSteps,1);

equi = zeros(100*numSteps,1);
counter = 1;
for loadFact = 1:numSteps

if loadFact == numSteps
    
    random = 1;
    
end
[K,Kmm,Kbb,~,~] = initK(numDofs);
KNL= K;

%Get stiffness at Tangent

[K,~] = computeTanK(ngp,pts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,elRange,elConn,numElems,sctr,K,KNL,t);

all_detK(loadFact) = det(K);
if applyLaGrange == 1
[K,f_ext] = applyLaGrangeBC(Xi_b,Xi_dof,Xi_b_val,p,Xi,Weights,K,f_ext(1:numDofs));
end

all_detK(loadFact) = det(K);
%Force increment
f_ext = f_ext + f_step;

%Disp increment
bc(:,2) = bc(:,2) + initBC(:,2);

tol = 0.001;
%end
maxIter = 100;
max_delta_F_steps = zeros(maxIter,1);
max_delta_f = 0;
iter = 1;
iterBC = bc;
old_max_delta_f = 1;
diff = tol+1;
while diff > tol
f = f_ext - fi;

[delta_a_temp,r_temp] = solveq(K,f,iterBC);
if applyLaGrange == 1
rest = sum(delta_a_temp(numDofs+1:numDofs+length(Xi_b)));
else
rest = 0;
end
[delta_a,delta_ax,delta_ay,r,f_ext] = getDispForce(delta_a_temp,r_temp,f_ext,numDofs);

a_tot = a_tot + delta_a;

% Update Geometry
[defPts] = updateGeometry(pts,delta_a, 1);

% Compute internal force vector
[fi_current] = getFinternal(ngp,defPts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,numElems,elConn,elRange,sctr,t);

%Update internal forces
fi = fi_current;
loadValInt(loadFact) = abs(fi(1));
%Remove BC forces
for tempBC = 1:length(iterBC)
    if iterBC(tempBC,2) == 0
    fi(iterBC(:,1)) = 0;
    else
    fi(iterBC(:,1)) = 0;
    end
    iterBC(tempBC,2) = 0;
end
%Calculate difference
delta_f = f_ext - fi;
equi(counter) = max(abs(delta_f));
counter = counter+1;
%neglect constrained nodes
delta_f(c_dofs) = 0;

old_max_delta_f = max_delta_f;
max_delta_f = abs(max(delta_f));

max_delta_F_steps(iter) = max_delta_f;

if applyLaGrange == 1
[~,f_ext] = applyLaGrangeBC(Xi_b,Xi_dof,zeros(1,length(Xi_b)),p,Xi,Weights,zeros(numDofs),f_ext);
[~,fi] = applyLaGrangeBC(Xi_b,Xi_dof,zeros(1,length(Xi_b),1),p,Xi,Weights,zeros(numDofs),fi);
end

%Update for next step
pts = defPts;
if iter == 1

    iter_one(value) = abs(a_tot(22));
    loadVal(value) = max(abs(f_ext));
    value = value+1;
    
end
if iter == maxIter
    break
end
iter = iter+1;
diff = max_delta_f;

end

all_a(:,loadFact) = a_tot;
diff_f = max_delta_F_steps(1:iter);
max_delta_f;
iter;

% if iter == maxIter
%     print = 'Max iterations reached'
%     %iter
%     %max_delta_f
%     %loadFact
%     break
% end
loadFact;
all_iter(loadFact) = iter;


%max(abs(f_ext))/Pcrit
%toc
end
max(abs(f_ext))/Pcrit;
toc
all_iter;
% PLOT_--------------------------------------------------------------------------------



%Get final scaled deformation

figure(1)
hold on
plot(refPts(:,1),refPts(:,2),'k--o');
curve = nrbmak(refPts',Xi);
nrbplot(curve,100)
grid on
grid minor


for i = 5:50:numSteps
    
defPtsPlot = updateGeometry(refPts,all_a(:,i),1);


defCurve = nrbmak(defPts',Xi);
defCurvePlot = nrbmak(defPtsPlot',Xi);


nrbplot(defCurvePlot,100)

end


defPtsPlot = updateGeometry(refPts,all_a(:,numSteps),1);
defCurvePlot = nrbmak(defPtsPlot',Xi);


nrbplot(defCurvePlot,100)

plot(defPts(:,1),defPts(:,2),'k--o');

hold off

xlim([-2 12])
ylim([-1 6])

figure(6)
plot(1:counter,equi(1:counter),'r--x')
grid on
grid minor
figure(7)
hold on

figure(9)

plot(1:numSteps,all_iter,'rx--')

figure(8)
plot([0 abs(all_a(21,:))],[0 loadValInt'],'x--k')
hold on
plot([0 abs(all_a(21,:))],ones(length(all_a(21,:))+1,1)*Pcrit,'--k')
grid on
grid minor
hold off
% figure(8)
% Ux_pts = [refPts(:,1) defPts(:,1) - refPts(:,1)];
% Uy_pts = [refPts(:,1) defPts(:,2) - refPts(:,2)];
% 
% Ux_crv = nrbmak(Ux_pts',Xi);
% Uy_crv = nrbmak(Uy_pts',Xi);
% 
% [kappa_xi,epsilon_x] = computeStrain(0:0.1:1,numCpts,p,Xi,Weights,refPts,defPts);
% 
% sigma_xx = E*epsilon_x;
% 
% legend('Undeformed Control Polygon','Deformed Control Polygon','Undeformed Geo','Deformed Geo')
% 
% figure(2)
% nrbplot2(Ux_crv,100)
% legend('U_x')
% 
% figure(3)
% 
% nrbplot2(Uy_crv,100)
% legend('U_y')
% 
% figure(4)
% 
% plot(0:0.01:1,sigma_xx)
% %ylim([-20 1]);
% legend('sigma_x_x')
% 
% figure(5)
% 
% plot(0:0.01:1,kappa_xi,'-x')
% legend('kappa_x_i')

