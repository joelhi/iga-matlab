function [K] = assembleSequentialK(K1,K2)

size1 = size(K1,1);
size2 = size(K2,1);

zeroK = zeros(size1,size2);

K = [K1 zeroK;
     zeroK' K2];

end

