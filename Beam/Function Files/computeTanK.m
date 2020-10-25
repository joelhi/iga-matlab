function [K,KNL] = computeTanK(ngp,pts,refPts,numDofsPerNode,Xi,p,Weights,E,A,I,elRange,elConn,numElems,sctr,K,KNL,t)

         %Prepare Height Map.
        
        t_val = cell2mat(t(1));
        t_knot = cell2mat(t(2));
        t_p = cell2mat(t(3));
        
        t_pts = zeros(length(t_val),2);
        
        t_step = 1/(length(t_val)-1);
                
        for j = 1:length(t_val)
           
            t_pts(j,:) = [(j-1)*t_step t_val(j)];
                      
        end


for e = 1 : numElems

    xiE = elRange(e,:);
    
    conn = elConn(e,:);
    Ex = pts(conn,1);
    Ey = pts(conn,2);
    
    refEx = refPts(conn,1);
    refEy = refPts(conn,2);
    
    %Compute Linear K
    
    [Km,Kb] = computeK_Lin(ngp,Ex,Ey,numDofsPerNode,xiE,Xi,p,Weights,E,A,I,t_pts,t_knot,t_p);
    
    % Compute Non Linear K.
    
    [KmNL,KbNL] = computeNonLinearK(ngp,Ex,Ey,refEx,refEy,numDofsPerNode,xiE,Xi,p,Weights,E,A,I,t_pts,t_knot,t_p);
    
    KeL = Km + Kb;
    KeNL = KmNL + KbNL;
    
    Ke = KeL;% + KeNL;
    
    K(sctr(e,:),sctr(e,:)) = K(sctr(e,:),sctr(e,:)) + Ke;
    KNL(sctr(e,:),sctr(e,:)) = KNL(sctr(e,:),sctr(e,:)) + KeNL;
    
end