function [alpha,beta] = computeMetricC(a1,a2,a11,a12,a22,a3)

    alpha = zeros(2,2);
    
    alpha(1,1) = dot(a1,a1);
    alpha(1,2) = dot(a1,a2);
    alpha(2,1) = dot(a2,a1);
    alpha(2,2) = dot(a2,a2);

    beta = zeros(2,2);
    
    beta(1,1) = dot(a11,a3);
    beta(1,2) = dot(a12,a3);
    beta(2,1) = dot(a12,a3);
    beta(2,2) = dot(a22,a3);
    
    
end