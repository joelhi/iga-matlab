function [a_tot,delta_a,delta_r,r_LG] = solveEquations(a_tot,numDofs,freeDofs,K,f_e,f_i,Bu,Bcouple,a_c,u_c)
    
    constrDofs = setdiff(1:numDofs,freeDofs);

    if length(freeDofs) == numDofs
        
        if Bu ~= -1
            [K,f_e,f_i] = applyKinematicConstraints(Bu,K,f_e,f_i,u_c);
        end
        
        if Bcouple ~= -1
            [K,f_e,f_i] = applyKinematicConstraints(Bcouple,K,f_e,f_i);
        end
    
    f = f_e - f_i;
    
    [temp_a] = solveq(K,f);
    
    delta_a = temp_a(1:numDofs);
    lambda = temp_a(numDofs+1:end);
    delta_r = -1;
    
    a_tot = a_tot + delta_a;
    
    r_LG = computeReactionForces(K,temp_a,freeDofs);
    
    else
    
    fe_f = f_e(freeDofs);
    fi_f = f_i(freeDofs);
    
    % In K, freeDofs, constrDofs
    [Kff,Kfc,Kcc] = partK(freeDofs,constrDofs,K);
    
    f_eq = Kfc*a_c(:,2);
    
        if Bu ~= -1
            [Kff,fe_f,fi_f] = applyKinematicConstraints(Bu(freeDofs,:),Kff,fe_f,fi_f,u_c);
        end
        
        if Bcouple ~= -1
            [Kff,fe_f,fi_f] = applyKinematicConstraints(Bcouple(freeDofs,:),Kff,fe_f,fi_f);
        end
        
        
    f = fe_f - fi_f;
    
    f_eq = [f_eq;
            zeros(length(f)-length(f_eq),1)];
            
    
    [temp_a] = solveq(Kff,f-f_eq);

    
    
    delta_a = temp_a(1:length(freeDofs));
    
    delta_r = Kfc'*temp_a(1:length(freeDofs)) + Kcc*a_c(:,2);
    
    a_tot(freeDofs) = a_tot(freeDofs) + delta_a;
    a_tot(constrDofs) = a_tot(constrDofs) + a_c(:,2);
    
    r_LG = computeReactionForces(Kff,temp_a,freeDofs);
    
    end


end
