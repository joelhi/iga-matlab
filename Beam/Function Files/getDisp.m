function [a,ax,ay] = getDispForce(a_temp,r_temp,f_temp,numDofs)

a = a_temp(1:numDofs);
a = a_temp(1:numDofs);

numCpts = numDofs/2; 

ax = zeros(numCpts,1);
ay = zeros(numCpts,1);

j = 1;
for i = 1:2:length(a)
   
    ax(j) = a(i);
    ay(j) = a(i+1);
    
    j = j + 1;
end


end