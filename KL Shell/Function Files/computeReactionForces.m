function [r_LG] = computeReactionForces(Kff,delta_a_f,freeDofs)

    if size(Kff,1) == freeDofs
        
        r_LG = zeros(length(freeDofs),1);

    else
    Bk = Kff(length(freeDofs):end,1:length(freeDofs));
    
    r_LG = Bk' * delta_a_f(length(freeDofs):end);
  
    end
end