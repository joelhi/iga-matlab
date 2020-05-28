function [stressVal, Sx, Sy ,Sz] = getPatchStress(a, pts, Xi, Eta, Xi_store, Eta_store, n, m, p, q,ep,weights,element,index,stressToPlot,layer,t_arr)
    
    % New control points
    at = a(1:3:end);
    au = a(2:3:end);
    an = a(3:3:end);
    
    defpts(:,1) = pts(:,1)+at;
    defpts(:,2) = pts(:,2)+au;
    defpts(:,3) = pts(:,3)+an;

    all_pts = [defpts(:,1),defpts(:,2),defpts(:,3)];
    all_refPts = [pts(:,1), pts(:,2), pts(:,3)];
    
    stressVal = zeros(length(Xi_store),length(Eta_store));
    
    Sx = zeros(length(Xi_store),length(Eta_store));
    Sy = zeros(length(Xi_store),length(Eta_store));
    Sz = zeros(length(Xi_store),length(Eta_store));
    
    B1 = zeros(4,n,m);
    
    B1(1:3,:) = defpts';
    B1(4,:) = weights;
    
    kk = 1;
    % Construct NURBS surface
    for i = 1:length(Xi_store)
        for j = 1:length(Eta_store)

            xi  = Xi_store(i);  % Get current xi
            eta = Eta_store(j); % Get current eta

            if nargin == 17
               
                ep(3) = t_arr(i,j);
                
            end
            
            % Find xi span and get b-spline basis
            xiSpan = FindSpanPS(n,p,xi,Xi);       % Find knot span
            N = BasisFun(xiSpan,xi,p,Xi);       % Get only the non-zero basis functions (there are p + 1 non-zero basis functions).

            % Find eta span and get b-spline basis
            etaSpan = FindSpanPS(m,q,eta,Eta);    % Find knot span
            M = BasisFun(etaSpan,eta,q,Eta);      % Get only the non-zero basis functions (there are q + 1 non-zero basis functions).
            
            % Construct the bi-variate b-spline basis using the outer product
            NM=N'*M;

            % Get relevant control points for current non-zero basis funcitons.
            X=reshape(B1(1,xiSpan-p+1:xiSpan+1,etaSpan-q+1:etaSpan+1),p+1,q+1);
            Y=reshape(B1(2,xiSpan-p+1:xiSpan+1,etaSpan-q+1:etaSpan+1),p+1,q+1);
            Z=reshape(B1(3,xiSpan-p+1:xiSpan+1,etaSpan-q+1:etaSpan+1),p+1,q+1);
            w=reshape(B1(4,xiSpan-p+1:xiSpan+1,etaSpan-q+1:etaSpan+1),p+1,q+1);

            R = NM .* w;            % Multiply Bspline basis by weights
            R = R / sum(sum(R));    % Divide by total W -> NURBS basis
            
            % find element
            
            for k = 1:length(index)
               
                if index(k,1) == xiSpan-(p-1) && index(k,2) == etaSpan-(q-1)  
                    
                e = k;
                break;
                end
                
            end

            
            conn   = element(e,:);          
            defpts     = all_pts(conn,:);
    
            refPts     = all_refPts(conn,:);

                        [R2, dRdxi, dRdeta, dR2dxi, dR2det, dR2dxe] = ...
            NURBS2DBasis2ndDers([xi; eta],p,q,Xi,Eta,weights);
        
            BasisFunc = [R2; dRdxi; dRdeta; dR2dxi; dR2det; dR2dxe];
            [S_n,S_m,E1,E2] = computePatchStress(BasisFunc,refPts,defpts,ep); 

            Stress = [S_n(1:2)/ep(3)+S_m(1:2)*layer/(ep(3)^2/6) S_n(3)/ep(3) S_m(3)]'/1e6;
            
            % Current point
            x = sum(sum(X.*R));
            y = sum(sum(Y.*R));
            z = sum(sum(Z.*R));

            % Save current surface point
            Sx(i,j) = x;
            Sy(i,j) = y;
            Sz(i,j) = z;

            stressVal(i,j) = Stress(stressToPlot);
            
            kk = kk + 1;
        end
    end

    


end