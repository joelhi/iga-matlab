function [N, dN, coord, J, normal]=nurbshaped(e,u_hat,v_hat,u_knot,v_knot,b_net,p,q,mcp,ncp,element_nod,coord_ij, corient)

%calculate the shape function and second derivatives


nsd = 2; 
% ------------------------------------------------------------------
%     get nurbs coordinates forml local node 1;
ni = coord_ij(element_nod(e,end),1);
nj = coord_ij(element_nod(e,end),2);
%     get u and v coordinates of integration point;
u =((u_knot(ni+1)-u_knot(ni))*u_hat +u_knot(ni+1) + u_knot(ni))/2;
v =((v_knot(nj+1)-v_knot(nj))*v_hat +v_knot(nj+1) + v_knot(nj))/2;

%     evaluate 1d size functions and derivatives each direction;
M=dersbasisfuns(ni,p,mcp,u,u_knot); % calculate in u direction
P=dersbasisfuns(nj,q,ncp,v,v_knot); % calculate in v direction

coord = zeros(nsd,1);
w= zeros((p+1)*(q+1), 1);
cpts = zeros((p+1)*(q+1), 2);
B = zeros(1, (p+1)*(q+1));
dBxi = zeros(1, (p+1)*(q+1));
dBeta = zeros(1, (p+1)*(q+1));

N = zeros(1, (p+1)*(q+1));
dN = zeros(2, (p+1)*(q+1));

% Form tensor products
icount = 0;
 for j = 1:q+1
    for i = 1:p+1
        icount = icount+1;
     
        w(icount) = b_net(ni+i-p-1,nj+j-q-1,nsd+1);
        cpts(icount, 1) = b_net(ni+i-p-1,nj+j-q-1,1);
        cpts(icount, 2) = b_net(ni+i-p-1,nj+j-q-1,2);
    %   basis functions;    
        B(icount) = M(1,i)*P(1,j);
        dBxi(icount) = M(2,i)*P(1,j);
        dBeta(icount) = M(1,i)*P(2,j);               
    end
 end 

% Multiply each B-spline function with corresponding weight
N(1,:) = B .* w';
dN(1,:) = dBxi .* w';
dN(2,:) = dBeta .* w';

% Compute the sums of B-spline functions
w_sum = sum(N(1,:));
dw_xi = sum(dN(1,:));
dw_eta = sum(dN(2,:));

% Compute NURBS basis functions and its first and second derivatives in
% local coordinates
dN(1,:) = dN(1,:)/w_sum - N*dw_xi/w_sum^2;
dN(2,:) = dN(2,:)/w_sum - N*dw_eta/w_sum^2;
N = N/w_sum;

% calculate coordinates in physical space
icount = 0;
for j=1:q+1
    for i=1:p+1
        icount = icount+1;
        coord(1) = coord(1) + N(icount)*cpts(icount, 1);
        coord(2) = coord(2) + N(icount)*cpts(icount, 2);
    end
end

% Compute Jacobian matrix
dxdxi = [dN(1,:)*cpts(:,1), dN(2,:)*cpts(:,1); dN(1,:)*cpts(:,2), dN(2,:)*cpts(:,2)];
J = det(dxdxi);

% Solve for first derivatives in global coordinates 
dN = dxdxi'\dN;

%  computation of normal, if function argument corient is given
normal = [0,0];
if nargin==13
    if(corient==1)
      nor(1) = dxdxi(2,1);% dy/dxi 
      nor(2) = -dxdxi(1,1);%dx/dxi
    elseif(corient==2) 
      nor(1) = dxdxi(2,2);
      nor(2) = -dxdxi(1,2);
    elseif(corient==3) 
      nor(1) = -dxdxi(2,1);
      nor(2) = dxdxi(1,1);
    else  
      nor(1) = -dxdxi(2,2);
      nor(2) = dxdxi(1,2);
    end

    tmp = sqrt(nor(1)^2 + nor(2)^2);
    normal = nor/tmp; % normal vector in two dimensions
end