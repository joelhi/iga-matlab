function [B,lambda] = applyLaGrangeBC2D(param,Xi_dof,Xi_b_val,n,m,p,q,Xi,Eta,index,element,weights,numDofs)

  
  B = zeros(numDofs,size(param,1));
  

  
  for i = 1:length(param)
  
      tempXi = param(i,1);
      tempEta = param(i,2);
      tempDof = Xi_dof(i);
      
      % Find element of span
      % Xi/Eta span
      spanXi = FindSpanPS(n,p,tempXi,Xi)-(p-1);
      spanEta = FindSpanPS(m,q,tempEta,Eta)-(q-1);
      
      for j = 1:length(index)     
          if index(j,1) == spanXi && index(j,2) == spanEta         
              e = j;
              break
          end
      end
      
      ptsInd = element(e,:);
      
      
      % Get dofs/locations for values in stiffness matrix
      
      currentDofs = zeros(length(ptsInd),1);
      
      for j = 1:size(element,2)
      
          if tempDof == 1
                currentDofs(j) = 3*ptsInd(j)-2;
          elseif tempDof == 2
                currentDofs(j) = 3*ptsInd(j)-1;
          elseif tempDof == 3
                currentDofs(j) = 3*ptsInd(j);
          end      
              
      end
      
      
      [R, ~, ~, ~, ~, ~] = NURBS2DBasis2ndDers(...
                                [tempXi,tempEta], p, q, Xi,Eta, weights);
                            
      B(currentDofs,i) = R';
      lambda = Xi_b_val;
      
      
  end

      
end