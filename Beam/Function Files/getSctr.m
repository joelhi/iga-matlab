function [numCpts,numDofsPerNode,numDofs,nodesPerElem,sctr] = getSctr(elConn,pts,numElems)



numCpts = length(pts);
numDofsPerNode = length(pts(1,:));

numDofs = numCpts*numDofsPerNode;

nodesPerElem = length(elConn(1,:));

sctr = zeros(numElems,nodesPerElem*numDofsPerNode);


for i = 1 : numElems
    
for j = 1 : nodesPerElem
   
    sctr(i,(2*j)-1) = (2*elConn(i,j))-1;
    sctr(i,(2*j)) = (2*elConn(i,j));
    
end
    
end

end