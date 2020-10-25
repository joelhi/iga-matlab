function [BasisFun,GpW,GpJ,u,v] = BasisFunctionsEvaluate(p,q,ngp,numOfEl,index,elRangeXi,elRangeEta,Xi,Eta,weights)

BasisFun = zeros((p+1)*(q+1),6,ngp*ngp,numOfEl);
GpW = zeros(ngp*ngp);
GpJ = zeros(ngp*ngp,numOfEl);
u = zeros(ngp*ngp,numOfEl);
v = zeros(ngp*ngp,numOfEl);

   for e = 1 : numOfEl

    xiID   = index(e,1);
    etaID  = index(e,2);
    xiE    = elRangeXi(xiID,:);
    etaE   = elRangeEta(etaID,:); 
    
    [Q,W] = GetIntegrationConstants(ngp);
    
    gp = 1;
    for i = 1:ngp
        Q1 = Q(i);
        w1 = W(i);
        
        for j = 1:ngp
            
            Q2 = Q(j);
            w2 = W(j);
            xiVal      = parent2ParametricSpace(xiE, Q1);
            etaVal     = parent2ParametricSpace(etaE,Q2);
            J1      = jacobianPaPaMapping(xiE,etaE);
            
            u(gp,e) = xiVal;
            v(gp,e) = etaVal;
            
            [R, dRdxi, dRdeta, dR2dxi, dR2det, dR2dxe] = ...
            NURBS2DBasis2ndDers([xiVal; etaVal],p,q,Xi,Eta,weights');    
            
            BasisFun(:,:,gp,e) = [R' dRdxi' dRdeta' dR2dxi' dR2det' dR2dxe'];
            GpW(gp) =  w1*w2;
            GpJ(gp,e) = J1;
            gp = gp + 1;           
        
        end        
    end
   end

end