%
% Modified to codes Matlab by :
% Hung Nguyen Xuan
%
% Faculty of Mathematics & Informatics, University of Natural Sciences
% Vietnam   National University–HCM

function [shb,phyc,normal,detjb]=nurbedge(coord,ni,nj,u_knot,v_knot,b_net,p,q,mcp,ncp,corient)
%global nsd nshl p q mcp ncp;

%global nsd nshl p q mcp ncp;

nsd = 2; 
nshl = (p+1)*(q+1);


% coord: Gauss point in local coordinates [-1 1]
% ni: element number in u-direction
% nj: element number in v-direction
% u_knot, v_knot: knot vectors
% b_net: control points and weights
% p,q: polynomial degrees
% mcp, ncp: # of basis functions/control points
% corient: orientation of the edge (1-bottom, 2-right, 3-top, 4-left)
%ifac,p,q,r,mcp,ncp,ocp,gp(igaussb),gp(igaussb),ni,nj,nk,u_knot,v_knot,w_knot,b_net,face_or
%nsd=3;
%nshl=(p+1)*(q+1)*(r+1);
normal=zeros(1,nsd);
% the details are dependent of face orientation.;

if(corient==1)
   u_hat = coord;
   v_hat = -1;  % -1.0 ensure w will be at beginning knot vector
   normal(2) = -1;

elseif(corient==2) 
   u_hat = 1;
   v_hat = coord;
   normal(1) = 1;

elseif(corient==3) 
   u_hat = coord;     % 1.0 ensures u will be at end of knot
   v_hat = 1; %  vector
   normal(2) = 1;

else%if(edge_or(f)==4) 
   u_hat = -1;
   v_hat = coord;
   normal(1) = -1;
end

%  find integration point in parameter space;

u =((u_knot(ni+1)-u_knot(ni))*u_hat +u_knot(ni+1) + u_knot(ni))/2;
v =((v_knot(nj+1)-v_knot(nj))*v_hat +v_knot(nj+1) + v_knot(nj))/2;

%   evaluate 1d size functions;

M=dersbasisfuns(ni,p,mcp,u,u_knot);% calculate in u direction
N=dersbasisfuns(nj,q,ncp,v,v_knot);% calculate in v direction

%     form basis functions;
icount = 0;
denom_sum = 0;
shb=zeros(1,nshl);
derv_sum_u = 0;
derv_sum_v = 0;
shbl=zeros(nshl,nsd);
phyc = zeros(nsd,1);

for j = 0:q
    for i = 0:p
        icount = icount+1;
        %  basis functions;
        shb(icount) = M(1,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j,nsd+1);
        denom_sum = denom_sum + shb(icount);
        shbl(icount,1) = M(2,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j,nsd+1);%u
        derv_sum_u = derv_sum_u + shbl(icount,1);
        shbl(icount,2) = M(1,p+1-i)*N(2,q+1-j)*b_net(ni-i,nj-j,nsd+1);%v
        derv_sum_v = derv_sum_v + shbl(icount,2);
    end
end
clear icount                   
% denom_sum1=denom_sum
% divide through by denominator;
 shbl(:,1) = shbl(:,1)/denom_sum -(shb(:)*derv_sum_u)/(denom_sum^2);
 shbl(:,2) = shbl(:,2)/denom_sum -(shb(:)*derv_sum_v)/(denom_sum^2);
 shb = shb/denom_sum;

%  now we calculate the face jacobian;
dxdxi=zeros(nsd,nsd);
icount = 0;

for j = 0: q
    for i = 0: p
       icount = icount + 1;
       phyc(1) = phyc(1) + b_net(ni-i,nj-j,1)*shb(icount);
       phyc(2) = phyc(2) + b_net(ni-i,nj-j,2)*shb(icount);
       dxdxi(1,1) = dxdxi(1,1) + b_net(ni-i,nj-j,1)*shbl(icount,1);
       dxdxi(1,2) = dxdxi(1,2) + b_net(ni-i,nj-j,1)*shbl(icount,2);
       dxdxi(2,1) = dxdxi(2,1) + b_net(ni-i,nj-j,2)*shbl(icount,1);
       dxdxi(2,2) = dxdxi(2,2) + b_net(ni-i,nj-j,2)*shbl(icount,2);
   end
end

%  jacobian of edge mapping;
if((corient==1)||(corient==3))
   ee = dxdxi(1,1)^2+dxdxi(2,1)^2;
else
   ee = dxdxi(1,2)^2+dxdxi(2,2)^2;
end
% Jacobian of face mapping
detjb = sqrt(ee);%Length of edge
