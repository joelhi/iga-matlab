function [K,f_e,f_i] = applyKinematicConstraints(B,K,f_e,f_i,u_c)

    
        [B,keepCols1] = removeZeroCols(B);
    
     if nargin == 4
         u_c = zeros(size(B,2),1);
     end
     
        u_c = u_c(keepCols1);
        
        additionalSize = size(K,1) - size(B,1);
        
        if additionalSize == 1
        
        K = [K B;
         B' zeros(size(B,2))];
        else
        
            K = [K [B; zeros(additionalSize,size(B,2))];
                 [B; zeros(additionalSize,size(B,2))]' zeros(size(B,2))];
            
        end
        
        
            
        f_e = [f_e;
               u_c];
    
    
        f_i = [f_i;
               zeros(size(u_c,1),1)];



end