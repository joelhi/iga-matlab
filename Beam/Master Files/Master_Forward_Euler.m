%% Beam Joel
clear all; clc;close all
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
rho         = 3;
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

%Connectivity

[elRange,elConn,index,numElems] = getConnectivity(Xi,p);

[numCpts,numDofsPerNode,numDofs,nodesPerElem,sctr] = getSctr(elConn,pts,numElems); 

% Control

delta_t = 0.000097;

bc = [2;
      22];

f_ext = zeros(numDofs,1); 
  
f_dof = 28;
f_val = -116000;

spring_end = [5 0];
springK = 600;

numSteps =6000;

%refPts are undeformed, pts are current, defPts are next step
refPts = pts;

acc = zeros(numDofs,1); % Acceleration
v   = zeros(numDofs,1); % Speed
k   = 0; % Kinetic energy
u   = zeros(numDofs,1);

max_v_steps = zeros(numSteps,1);
min_v_steps = zeros(numSteps,1);
i = 1;

%Gif
h = figure('units','normalized','outerposition',[0 0 1 1]);
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';

all_u = zeros(numDofs,numSteps);
f_ext(numDofs-1) = f_val;
f_ext(1) = -f_val;
f_ext(2:2:end) = 10000;
tic
for step = 1:numSteps
    
    if step == 1400
        f_ext(2:2:end) = 0;
    end
    
    %Assign Spring Force
    %f_val = (spring_end(1,1) - pts(11,1))*springK;
    
    step
    % Get internal Forces from previous step
    
    f_int = getFinternal(ngp,pts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,numElems,elConn,elRange,sctr,t);
    
    delta_f = f_ext - f_int;
    max(abs(delta_f));
    % set fixed nodes forces to zero
    
    delta_f(bc') = 0;
    
    % Mass matrix
    
    M = eye(numDofs) * rho;
    
    % Get accelerations
    
    acc = (M\delta_f);
    
    % Update Speed
    
    v = 0.99*v + acc*delta_t;
    
    %Update Positions
    delta_u = (v*delta_t)/2;
    
    delta_u(bc') = 0;
    
    %Total displacements
    u = u+delta_u;
    all_u(:,step) = u;
    %Update geometry
    defPts = updateGeometry(pts,delta_u,1);

    defPtsPlot = updateGeometry(pts,delta_u,1);
    
    %Plot new position
    tempCrv = nrbmak(defPts',Xi);
    nrbplot(tempCrv,100);
    hold on
    %plot(spring_end(1,1),spring_end(1,2),'xr')
    %plot([spring_end(1,1) pts(11,1)],[spring_end(1,2) pts(11,2)],'--r')
    plot(refPts(:,1),refPts(:,2),'k--o')
    hold on
    defPtsPlot = updateGeometry(refPts,u,1);
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
          imwrite(imind,cm,filename,'gif','DelayTime',1,'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','DelayTime',0.01,'WriteMode','append'); 
      end
    
    
    min_v_steps(step) = min(v);
    max_v_steps(step) = max(v);
    
    pts = defPts;
    i = i + 1;
    if max(abs(v)) < 0.1 && max(abs(acc)) < 0.1
        break;
    end
    
end
toc


figure(2)
plot(1:numSteps,max_v_steps)
hold on
plot(1:numSteps,min_v_steps)
hold off
grid on
grid minor


figure(3)
hold on
plot(refPts(:,1),refPts(:,2),'k--o');
curve = nrbmak(refPts',Xi);
nrbplot(curve,100)

grid on
grid minor
axis equal
xlim([0 12])
ylim([-3 2])
for i = 200:200:numSteps
    
defPtsPlot = updateGeometry(refPts,all_u(:,i),1);


defCurve = nrbmak(defPts',Xi);
defCurvePlot = nrbmak(defPtsPlot',Xi);


nrbplot(defCurvePlot,100)

end
axis equal
xlim([0 12])
ylim([-3 2])

defPtsPlot = updateGeometry(refPts,all_u(:,numSteps),1);
plot(defPtsPlot(:,1),defPtsPlot(:,2),'k--o');
defCurvePlot = nrbmak(defPtsPlot',Xi);

axis equal
xlim([0 12])
ylim([-3 2])
hold off


figure(4)
plot(all_u(12,:),1:numSteps)
grid on
grid minor

figure(5)
plot(-all_u(21,:),1:numSteps)
grid on
grid minor

