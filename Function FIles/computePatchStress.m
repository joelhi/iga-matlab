function [n,m,dE1,dE2] = computePatchStress(BasisFun,pts,defPts,ep,E1,E2)
    
            E = ep(1);
            nu = ep(2);
            t = ep(3);
            
            
    
            memStiff = E*t/(1-nu^2);
            benStiff = E*t^3/12/(1-nu^2);

            dRdxi = BasisFun(2,:);
            dRdeta = BasisFun(3,:);
            dR2dxi = BasisFun(4,:);
            dR2det = BasisFun(5,:);
            dR2dxe = BasisFun(6,:);

            
            jacob1  = [dRdxi; dRdeta] * pts;                % 2x3(?) matrix
            jacob2 = [dR2dxi; dR2det; dR2dxe] * pts;        % 3x3(?) matrix

            djacob1  = [dRdxi; dRdeta] * defPts;                % 2x3(?) matrix
            djacob2 = [dR2dxi; dR2det; dR2dxe] * defPts;        % 3x3(?) matrix
            
            
            
%             dxdxi = jacob1(1,1);         dydxi = jacob1(1,2);
%             dxdet = jacob1(2,1);         dydet = jacob1(2,2);
            
            % a1, a2 and a3 vectors (surface basis vectors)
            % a1 & a2 defines the tangent plane at a specific point
        
            a1    = jacob1(1,:);
            a2    = jacob1(2,:);
            a3    = cross(a1,a2); 
            norma = norm(a3);
            a3    = a3/norma; 
            
            a1_norm = a1/norm(a1);
            a1_perp = cross(a1_norm,a3);
            
            if nargin == 4
            
            E1 = a1_norm;
            E2 = a1_perp;
            
            end
            
            J1 = norma;
           
            da1    = djacob1(1,:);
            da2    = djacob1(2,:);
            da3    = cross(da1,da2); 
            dnorma = norm(da3);
            da3    = da3/dnorma;
            
            
            da1_norm = da1/norm(da1);
            da1_perp = cross(da1_norm,da3);
            
            dE1 = da1_norm;
            dE2 = da1_perp;
            
            dJ1 = dnorma;
            
            J = dJ1/J1;
            
            % Restriction of the metric tensor to the tangent plane 
            % also called the 1st fundamental form of the surface, 
            % given by its components
            a11   = jacob2(1,:);
            a22   = jacob2(2,:);
            a12   = jacob2(3,:);
            
            da11   = djacob2(1,:);
            da22   = djacob2(2,:);
            da12   = djacob2(3,:);
        
            % dot products of ai and ei    
        
            % R_I,2*a1 + R_I,1*a2 for all shape functions
            
            [alpha,beta] = computeMetricC(a1,a2,a11,a12,a22,a3);
            [dalpha,dbeta] = computeMetricC(da1,da2,da11,da12,da22,da3);
            
            %Calculate Strains
          
            %Metric Coefficient
     
            e11 = 1/2*(dalpha(1,1) - alpha(1,1));
            e22 = 1/2*(dalpha(2,2) - alpha(2,2));
            e12 = 1/2*(dalpha(1,2) - alpha(1,2));
            e21 = 1/2*(dalpha(2,1) - alpha(2,1));
            
            eHat = [e11 e12; e21 e22];        
         
            k11 = 1/2*(beta(1,1) - dbeta(1,1));
            k22 = 1/2*(beta(2,2) - dbeta(2,2));
            k12 = 1/2*(beta(1,2) - dbeta(1,2));
            k21 = k12;

            kHat = [k11 k12; k21 k22];
            
            % Map eHat and kHat to local unitized cartesian basis E1xE2
            
            e_ortho = shiftBase2D(eHat,a1,a2,E1,E2,1);
            
            k_ortho = shiftBase2D(kHat,a1,a2,E1,E2,1);
            
            C_ortho = [1 nu 0;
                       nu 1 0;
                       0 0 (1-nu)/2];
                   
            e_ortho_v = [e_ortho(1,1); e_ortho(2,2); e_ortho(1,2)+e_ortho(2,1)];
            
            k_ortho_v = [k_ortho(1,1); k_ortho(2,2); k_ortho(1,2)+k_ortho(2,1)];

            S_n = memStiff*(C_ortho * e_ortho_v);
            S_m = benStiff*(C_ortho * k_ortho_v);
            
            S_n_tensor = [S_n(1) S_n(3)*0.5; S_n(3)*0.5 S_n(2)];
            S_m_tensor = [S_m(1) S_m(3)*0.5; S_m(3)*0.5 S_m(2)];  
            
            % Scale to correct value
            
            scaled_S_n = S_n_tensor * 1/J;
            scaled_S_m = S_m_tensor * 1/J;
            % Transform to actual configuration
            
            sigma_n = shiftBase2D(scaled_S_n,E1,E2,dE1,dE2,1);
            sigma_m = shiftBase2D(scaled_S_m,E1,E2,dE1,dE2,1);
                        
            n = [sigma_n(1,1) sigma_n(2,2) 2*sigma_n(1,2)];
            m = [sigma_m(1,1) sigma_m(2,2) 2*sigma_m(1,2)];

            

end