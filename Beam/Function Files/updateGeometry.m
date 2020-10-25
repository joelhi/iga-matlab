function [defPts] = updateGeometry(pts,a, scaleFact)

numCpts = length(pts);
dim = size(pts,2);

defPts = zeros(numCpts,dim);

for i = 1:numCpts
    if dim == 2
    defPts(i,:) = pts(i,:) + scaleFact*[a((2*i)-1) a(2*i)];
    elseif dim == 3
    defPts(i,:) = pts(i,:) + scaleFact*[a((3*i)-2) a((3*i)-1) a(3*i)];    
    end
end

end