function [K_mod,f_mod] = applyLaGrangeBC(Xi_b,Xi_dof,Xi_b_val,p,Xi,Eta,Weights,K,f)

N = zeros(length(f),length(Xi_b));

for i = 1:length(Xi_b)

    Xi_b_temp = Xi_b(i);
    Xi_dof_temp = Xi_dof(i);
    

%spanIndex = FindSpanPS(numCpts-1,p,Xi_b_temp,Xi)+1;

for j = 1:length(Weights)
    
  if Xi_dof_temp == 1    
  N(2*j-1,i) = NURBSbasis(j, p, Xi_b_temp, Xi, Weights);
  elseif Xi_dof_temp == 2
  N(2*j,i) = NURBSbasis(j, p, Xi_b_temp, Xi, Weights);
  elseif Xi_dof_temp == 3
  
  end
  
end

end

  K_mod = [K N;
           N' zeros(length(Xi_b))];
  f_mod = [f;
           Xi_b_val'];
      



end