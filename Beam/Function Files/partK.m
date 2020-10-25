function [Kff,Kfc,Kcc] = partK(freeDofs,constrDofs,K)

    Kff = K(freeDofs,freeDofs);
    Kfc = K(freeDofs,constrDofs);
    Kcc = K(constrDofs,constrDofs);

end