function [a1_c,a2_c] = computeContravariant(a1,a2)

            a1a2 = [a1; a2]; 
            
            a1_c = a1a2 \ [1;0];
            a2_c = a1a2 \ [0;1]; 

end