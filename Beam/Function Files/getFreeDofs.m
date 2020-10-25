function [freeDofs] = getFreeDofs(numDofs,zero_bc,nonzero_bc)

freeDofs = zeros(numDofs-length(zero_bc)-length(nonzero_bc),1);
i = 1;
for dof = 1:numDofs
    dofConstrained = false;  
    for zeroDofId = 1:length(zero_bc)     
        zeroDof = zero_bc(zeroDofId);
        if dof == zeroDof
            dofConstrained = true;
        end
    end
    for nonzeroDofId = 1:length(nonzero_bc)        
        nonzeroDof = nonzero_bc(nonzeroDofId);
        if dof == nonzeroDof
            dofConstrained = true;
        end       
    end
    if dofConstrained
    else
    freeDofs(i) = dof;
    i = i+1;
    end
end

end