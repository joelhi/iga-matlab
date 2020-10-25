function [Km,Kb] = computeK_Lin(ngp,Ex,Ey,numDofsPerNode,xiE,Xi,p,Weights,E,A,I,t_pts,t_knot,t_p)

varyingT = 1;

 [Q,W] = GetIntegrationConstants(ngp);
    
    if varyingT == 0
    memStiff = E*A;
    benStiff = E*I;
    end
    
    elemPts = [Ex Ey];
    Km = zeros(numDofsPerNode*length(Ex),numDofsPerNode*length(Ex));
    Kb = zeros(numDofsPerNode*length(Ex),numDofsPerNode*length(Ex));
    fb = zeros(numDofsPerNode*length(Ex),1);
    
    for i = 1:ngp
       
        Q1 = Q(i);
        W1 = W(i);
        XiCoord      = parent2ParametricSpace(xiE, Q1);
        
        if varyingT == 1
            
            N_t = zeros(1,length(t_pts));
        
            t_weights = ones(1,length(t_pts));
            
           for j = 1:length(t_pts)
     
                     N_t(1,j) = NURBSbasis(j, t_p, XiCoord, t_knot, t_weights);
                 

           end
        
           t = N_t*t_pts(:,2);
           
           At = A*t; % A doubles as width in this case
           It = A*t^3/12;
           
           memStiff = E*At;
           benStiff = E*It;
           
           
        end
        
        J1 = jacobianPaMapping(xiE);
        
        [R, dRdxi, dR2dxi] = NURBS1DBasis2ndDers(XiCoord,p,Xi,Weights);
        
        jacob1 = dRdxi * elemPts;
        jacob2 = dR2dxi * elemPts;
        
        % Base Vector
        
        a1 = jacob1;
        
        a3 = [-a1(2),a1(1)];
        
        norma1 = norm(a1);
        norma3 = norm(a3);
        
        a3 = a3/norma3;
        J2 = norma1;
        
        % C "MATRIX"
        
        a11 = dot(a1,a1);
        
        C = 1/a11^2;
        
        noBasis = length(R);
        Bm = zeros(1,noBasis*numDofsPerNode);
        Bb = zeros(1,noBasis*numDofsPerNode);
        
        for k = 1:noBasis
            
            dRIdx = dRdxi(k);
            dRI2dx = dR2dxi(k);
            
            id = [(2*k)-1 2*k];
            
            Bm(1,id) = [dRIdx*a1(1) dRIdx*a1(2)];
            Bb(1,id) = [-dRI2dx*a3(1)+1/norma3*jacob2(2)*dRIdx -dRI2dx*a3(2)+1/norma3*jacob2(1)*dRIdx];         
            
        end      
        
            Km = Km + memStiff * Bm' * C * Bm * J1 * J2 * W1;
            Kb = Kb + benStiff * Bb' * C * Bb * J1 * J2 * W1;
    end

end