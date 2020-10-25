function [maxDiff,fDiff] = checkEq(freeDofs,f_i,f_e,r_LG)

        fi_f = f_i(freeDofs);
        fe_f = f_e(freeDofs);
    
        fi_res = abs(fi_f-fe_f);
    
        fDiff = abs(fi_res) - abs(r_LG);
        
        maxDiff = max(abs(fi_res) - abs(r_LG));
        
end