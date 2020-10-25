function [K,Kmm,Kbb,f,fi] = initK(numDofs)

K = zeros(numDofs, numDofs);
Kmm = zeros(numDofs, numDofs);
Kbb = zeros(numDofs, numDofs);
f = zeros(numDofs, 1);
fi = zeros(numDofs, 1);

end