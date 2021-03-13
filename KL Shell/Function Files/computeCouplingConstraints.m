function [B] = computeCouplingConstraints(u,v,dof,n,m,p,q,Xi1,Xi2,Eta1,Eta2,index1,index2,elementGlobal1,elementGlobal2,weights1,weights2,numDofs)



        B = zeros(numDofs,1);
    % Return column matrix for one couple.
    % Compute shape functions for each node.
    
      spanXi1 = FindSpanPS(n(1),p(1),u(1),Xi1)-(p(1)-1);
      spanEta1 = FindSpanPS(m(1),q(1),v(1),Eta1)-(q(1)-1);
      
      spanXi2 = FindSpanPS(n(2),p(2),u(2),Xi2)-(p(2)-1);
      spanEta2 = FindSpanPS(m(2),q(2),v(2),Eta2)-(q(2)-1);
      
      for j = 1:length(index1)     
          if index1(j,1) == spanXi1 && index1(j,2) == spanEta1         
              e1 = j;
              break
          end
      end
      
      for j = 1:length(index2)     
          if index1(j,1) == spanXi2 && index2(j,2) == spanEta2         
              e2 = j;
              break
          end
      end
     
      ptsInd1 = elementGlobal1(e1,:);
      ptsInd2 = elementGlobal2(e2,:);
      
      % Get dofs/locations for values in stiffness matrix
      
      currentDofs1 = zeros(length(ptsInd1),1);
      currentDofs2 = zeros(length(ptsInd2),1);
      
      for j = 1:size(elementGlobal1,2)
      
          if dof == 1
                currentDofs1(j) = 3*ptsInd1(j)-2;
          elseif dof == 2
                currentDofs1(j) = 3*ptsInd1(j)-1;
          elseif dof == 3
                currentDofs1(j) = 3*ptsInd1(j);
          end      
              
      end
      
      for j = 1:size(elementGlobal2,2)
      
          if dof == 1
                currentDofs2(j) = 3*ptsInd2(j)-2;
          elseif dof == 2
                currentDofs2(j) = 3*ptsInd2(j)-1;
          elseif dof == 3
                currentDofs2(j) = 3*ptsInd2(j);
          end      
              
      end
      
      
      [R1, ~, ~, ~, ~, ~] = NURBS2DBasis2ndDers([u(1),v(1)], p(1), q(1), Xi1,Eta1, weights1);
      
                            
      [R2, ~, ~, ~, ~, ~] = NURBS2DBasis2ndDers([u(2),v(2)], p(2), q(2), Xi2,Eta2, weights2);
      B(currentDofs1) = B(currentDofs1) + R1';
      B(currentDofs2) = B(currentDofs2)-R2';
      
    
    
    
    


end
