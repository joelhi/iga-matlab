function [e_trans] = shiftBase2D(e,a1,a2,b1,b2,format)

% This function shifts the 2d tensor e from basis a1xa2 to b1xb2
%   
%
%
    
    e_trans = zeros(2,2);
    
    if format == 0
    % Getr
    
    a1b1 = dot(a1,b1);
    a1b2 = dot(a1,b2);
    a2b1 = dot(a2,b1);
    a2b2 = dot(a2,b2);
    
    ab = [a1b1 a1b2;a2b1 a2b2];
    
    
    for l = 1:2
        for k = 1:2
            for j = 1:2
                for i = 1:2
                    e_trans(k,l) = e_trans(k,l) + e(i,j)*ab(i,k)*ab(j,l); 
                end
            end
        end
    end
    
    end
    
    if format == 1
    % Getr
    
    [a1_c,a2_c] = computeContravariant(a1,a2);
    
    a1b1 = dot(a1_c,b1);
    a1b2 = dot(a1_c,b2);
    a2b1 = dot(a2_c,b1);
    a2b2 = dot(a2_c,b2);
    
    ab = [a1b1 a1b2;a2b1 a2b2];
    
    
    for l = 1:2
        for k = 1:2
            for j = 1:2
                for i = 1:2
                    e_trans(k,l) = e_trans(k,l) + e(i,j)*ab(i,k)*ab(j,l); 
                end
            end
        end
    end
    
    
    
    
    end
    
    if format == 2
    % Getr
    
    [b1_c,b2_c] = computeContravariant(b1,b2);
    
    a1b1 = dot(a1,b1_c);
    a1b2 = dot(a1,b2_c);
    a2b1 = dot(a2,b1_c);
    a2b2 = dot(a2,b2_c);
    
    ab = [a1b1 a1b2;a2b1 a2b2];
    
    
    for l = 1:2
        for k = 1:2
            for j = 1:2
                for i = 1:2
                    e_trans(k,l) = e_trans(k,l) + e(i,j)*ab(i,k)*ab(j,l); 
                end
            end
        end
    end
    
    end
    
    if format == 3
    % Getr
    
    [b1_c,b2_c] = computeContravariant(b1,b2);
    [a1_c,a2_c] = computeContravariant(a1,a2);
    
    a1b1 = dot(a1_c,b1_c);
    a1b2 = dot(a1_c,b2_c);
    a2b1 = dot(a2_c,b1_c);
    a2b2 = dot(a2_c,b2_c);
    
    ab = [a1b1 a1b2;a2b1 a2b2];
    
    
    for l = 1:2
        for k = 1:2
            for j = 1:2
                for i = 1:2
                    e_trans(k,l) = e_trans(k,l) + e(i,j)*ab(i,k)*ab(j,l); 
                end
            end
        end
    end
    
    end
    
end