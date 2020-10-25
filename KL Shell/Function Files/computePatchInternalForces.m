function [fi_current] = computePatchInternalForces(pts,refPts,sctr,element,ngp,ep,BasisFun,GpW,GpJ,t_arr)

fi_current = zeros(3*length(pts(:,1)),1);

numOfEl = size(element,1);

if nargin == 9

   t = ep(3);
    
end


 for e = 1 : numOfEl
     
    conn   = element(e,:);          
    defEx     = pts(conn,1);
    defEy     = pts(conn,2);
    defEz     = pts(conn,3);
    
    Ex     = refPts(conn,1);
    Ey     = refPts(conn,2);
    Ez     = refPts(conn,3);
     
    E = ep(1);
    nu = ep(2);
   
    elPts = [Ex Ey Ez];
    defElPts = [defEx defEy defEz];
    fi_m = zeros(3*length(Ex),1);
    fi_b = zeros(3*length(Ex),1);
    gp = 0;
    
    for i = 1:ngp
        
        for j = 1:ngp
            
            gp = gp+1;
            
            if nargin == 10
            t = t_arr(gp,e);
            end
            
            memStiff = E*t/(1-nu^2);
            benStiff = E*t^3/12/(1-nu^2);
            
            R = BasisFun(:,1,gp,e)';
            dRdxi = BasisFun(:,2,gp,e)';
            dRdeta = BasisFun(:,3,gp,e)';
            dR2dxi = BasisFun(:,4,gp,e)';
            dR2det = BasisFun(:,5,gp,e)';
            dR2dxe = BasisFun(:,6,gp,e)';
            
            w12 = GpW(gp);
            J1 = GpJ(gp,e);
            
            jacob1  = [dRdxi; dRdeta] * elPts;                % 2x3(?) matrix
            jacob2 = [dR2dxi; dR2det; dR2dxe] * elPts;        % 3x3(?) matrix

            djacob1  = [dRdxi; dRdeta] * defElPts;                % 2x3(?) matrix
            djacob2 = [dR2dxi; dR2det; dR2dxe] * defElPts;        % 3x3(?) matrix
            
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
            
            
            a1_perp = a2-(dot(a2,a1_norm)*a1_norm);
            norm_a1perp = norm(a1_perp);
            a1_norm_perp = a1_perp/norm_a1perp;
            
            E1 = a1_norm;
            E2 = a1_norm_perp;
            
            J2 = norma;
            
            da1    = djacob1(1,:);
            da2    = djacob1(2,:);
            da3    = cross(da1,da2); 
            dnorma = norm(da3);
            da3    = da3/dnorma;
         
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
        
            da1e1  = da1(1); da1e2  = da1(2); da1e3  = da1(3);
            da2e1  = da2(1); da2e2  = da2(2); da2e3  = da2(3);
        
            % R_I,2*a1 + R_I,1*a2 for all shape functions
            
            % Mapping for voight notation
            
            [a1_c,a2_c] = computeContravariant(a1,a2);
    
            a1b1 = dot(a1_c,E1);
            a1b2 = dot(a1_c,E2);
            a2b1 = dot(a2_c,E1);
            a2b2 = dot(a2_c,E2);
            
            ab = [a1b1 a1b2;a2b1 a2b2];
            
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
     
             % R_I,2*a1 + R_I,1*a2 for all shape functions
             
            % Cross Products
            
            xda11da2 = cross(da11,da2);
            dotda3da11 = dot(da3,da11);
            xda1da11 = cross(da1,da11);
            
            xda22da2 = cross(da22,da2);
            dotda3da22 = dot(da3,da22);
            xda1da22 = cross(da1,da22);
            
            xda12da2 = cross(da12,da2);
            dotda3da12 = dot(da3,da12);
            xda1da12 = cross(da1,da12);
            
            xda2da3 = cross(da2,da3);
            xda3da3 = cross(da3,da3);
            
            noBasis = length(R);
            dRIa    = zeros(3,noBasis);
            for jj=1:noBasis
                dRIa(:,jj) = dRdeta(jj)*da1 + dRdxi(jj)*da2;
            end
            
            Bm = zeros(3,noBasis*3);
            Bb = zeros(3,noBasis*3);
            for k = 1:noBasis
                dRIdx = dRdxi (k);
                dRIdy = dRdeta(k);

                id    = (k-1)*3+1:3*k;
                
                Bm1_1 = [dRIdx*da1e1; dRIdy*da2e1; dRIa(1,k)];
                Bm2_1 = [dRIdx*da1e2; dRIdy*da2e2; dRIa(2,k)];
                Bm3_1 = [dRIdx*da1e3; dRIdy*da2e3; dRIa(3,k)];
                
                Bm1 = shiftBase2DVoigth(Bm1_1,ab,1);
                Bm2 = shiftBase2DVoigth(Bm2_1,ab,1);
                Bm3 = shiftBase2DVoigth(Bm3_1,ab,1);         
                        
                Bm(:,id) = [Bm1' Bm2' Bm3'];             

                BI1_alt = -dR2dxi(k)*da3 + 1/norma*(dRIdx*xda11da2 + dRIdy*xda1da11 + ...
                              dotda3da11*(dRIdx*xda2da3 + dRIdy*xda3da3)); 

                BI2_alt = -dR2det(k)*da3 + 1/norma*(dRIdx*xda22da2 + dRIdy*xda1da22 + ...
                              dotda3da22*(dRIdx*xda2da3 + dRIdy*xda3da3)); 

                BI3_alt = -dR2dxe(k)*da3 + 1/norma*(dRIdx*xda12da2 + dRIdy*xda1da12 + ...
                              dotda3da12*(dRIdx*xda2da3 + dRIdy*xda3da3)); 

                          
                B1_alt =   [BI1_alt(1); BI2_alt(1); 2*BI3_alt(1)];
                B2_alt =   [BI1_alt(2); BI2_alt(2); 2*BI3_alt(2)];
                B3_alt =   [BI1_alt(3); BI2_alt(3); 2*BI3_alt(3)];
                
                Bb1 = shiftBase2DVoigth(B1_alt,ab,1);
                Bb2 = shiftBase2DVoigth(B2_alt,ab,1);
                Bb3 = shiftBase2DVoigth(B3_alt,ab,1);
               
                Bb(:,id) = [Bb1' Bb2' Bb3'];
                
            end
            
            fi_m = fi_m + Bm' * S_n * J1 * J2 * w12;  
            fi_b = fi_b + Bb' * S_m * J1 * J2 * w12;
            
        end
        
    end

    %Assemble
    
    fi_e = fi_m + fi_b;
    
    % Assemble
    fi_current(sctr(e,:),1) = fi_current(sctr(e,:),1) + fi_e;
    
 end
end