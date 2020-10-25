function [fi_current] = getFinternal(ngp,pts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,numElems,elConn,elRange,sctr,t)

         %Prepare Height Map.
        
        t_val = cell2mat(t(1));
        t_knot = cell2mat(t(2));
        t_p = cell2mat(t(3));
        
        t_pts = zeros(length(t_val),2);
        
        t_step = 1/(length(t_val)-1);
                
        for j = 1:length(t_val)
           
            t_pts(j,:) = [(j-1)*t_step t_val(j)];
                      
        end

numDofs = numDofsPerNode*length(Weights);

fi_current = zeros(numDofs,1);

for e = 1:numElems
   
    xiE = elRange(e,:);
    
    conn = elConn(e,:);
    
    Ex = refPts(conn,1);
    Ey = refPts(conn,2);
    
    defEx = pts(conn,1);
    defEy = pts(conn,2);
    
    [fm,fb] = computeInternalF(ngp,Ex,Ey,defEx,defEy,numDofsPerNode,xiE,Xi,p,Weights,E,A,I,t_pts,t_knot,t_p);
   
    fie = fm + fb;
    fi_current(sctr(e,:)) = fi_current(sctr(e,:)) + fie;
end

end


