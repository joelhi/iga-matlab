function [kappa_xi,epsilon_x,nVecRef,nVecDef,A_pos,a_pos] = computeStrain(xi_coord,numCpts,p,Xi,Weights,refPts,currentPts)
% Normal Strain

%Evaluate Derivatives

dNdxi = zeros(length(xi_coord),numCpts);
dN2dxi = zeros(length(xi_coord),numCpts);
N = zeros(length(xi_coord),numCpts);

ind = 1;

for i = 1:length(xi_coord)
    
    xi_p = xi_coord(i);
    
    span_ind = FindSpan(numCpts-1,p,xi_p,Xi);    
    
    [R, dRdxi, dR2dxi] = NURBS1DBasis2ndDers(xi_p,p,Xi,Weights);
    dN2dxi(ind,span_ind-p+1:span_ind+1) = dR2dxi;
    dNdxi(ind,span_ind-p+1:span_ind+1) = dRdxi;
    N(ind,span_ind-p+1:span_ind+1) = R;
    
    ind = ind + 1;
end

A_pos = N*refPts;
a_pos = N*currentPts;

A_xi = dNdxi*refPts;
a_xi = dNdxi*currentPts;

dot_Axi = zeros(length(A_xi),1);
dot_axi = zeros(length(A_xi),1);

A_xixi = dN2dxi*refPts;
a_xixi = dN2dxi*currentPts;

A3_ = [-A_xi(:,2),A_xi(:,1)];
a3_ = [-a_xi(:,2),a_xi(:,1)];

A3_norm = normr(A3_);
a3_norm = normr(a3_);

nVecDef = a3_norm;
nVecRef = A3_norm;

b_xi = zeros(length(A_xixi),1);
B_xi = zeros(length(A_xixi),1);

for i = 1:length(A_xi)
    
    dot_Axi_temp = (A_xi(i,1)*A_xi(i,1)+A_xi(i,2)*A_xi(i,2));
    dot_axi_temp = (a_xi(i,1)*a_xi(i,1)+a_xi(i,2)*a_xi(i,2));
    
    dot_Axi(i) = dot_Axi_temp/dot_Axi_temp;
    dot_axi(i) = dot_axi_temp/dot_Axi_temp;    
    
    b_xi_temp = a_xixi(i,1)*a3_norm(i,1) + a_xixi(i,2)*a3_norm(i,2);
    B_xi_temp = A_xixi(i,1)*A3_norm(i,1) + A_xixi(i,2)*A3_norm(i,2);
    
    b_xi(i) = b_xi_temp/dot_Axi_temp;
    B_xi(i) = B_xi_temp/dot_Axi_temp;
end

epsilon_x = 1/2*(dot_axi - dot_Axi);

kappa_xi = B_xi - b_xi;