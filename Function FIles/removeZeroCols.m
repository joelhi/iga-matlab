function [Bmod,keepCols] = removeZeroCols(B)

    keepCols = [size(B,2)];

    for i = 1:size(B,2)
        
        if max(abs(B(:,i))) > 1e-8
            
            keepCols(i) = i;
        
        end
        
    end

    keepCols = keepCols(keepCols ~= 0);
    
    Bmod = B(:,keepCols);

end