function [elRange,elConn,index,numElems] = getConnectivity(Xi,p) 
    
uniqueXiKnots       = unique(Xi);
numElems            = length(uniqueXiKnots)-1;

[elRange,elConn] = buildConnectivityPS(p,Xi,numElems);

index  = zeros(numElems,1);

for i = 1:numElems  
    index(i,1) = i;
end


end