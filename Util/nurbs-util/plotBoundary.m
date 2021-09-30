%script plotboundary.m
%1d plots of the boundary of the solution vs. exact solution (if available)


NumPoints = 31;  %number of evaluation points on each element


for i=1:size(dirichlet,1)
    localindex = 0;
    elementind = dirichlet(i,1);
        
    curnod = element_nod(elementind, nument);
    ni = coord_ij(curnod, 1);
    nj = coord_ij(curnod, 2);
    
    corient = dirichlet(i, 5);
    
    EvalPoints = linspace(-1,1,NumPoints);
    PlotPointsX = zeros(1, NumPoints);
    PlotPointsY = zeros(1, NumPoints);
    PhysPoints = zeros(1, NumPoints);
    AnalSolPointsX = zeros(1, NumPoints);
    AnalSolPointsY = zeros(1, NumPoints);
    ErrPointsX = zeros(1, NumPoints);
    ErrPointsY = zeros(1, NumPoints);
       
    for iEval = 1:NumPoints
        
        [R, coord, normal, J] = nurbedge(EvalPoints(iEval), ni, nj, knotu, knotv, b_net, p, q, lenu, lenv, corient);                

        [displacement, stress] = gbeam(coord(1), coord(2), P, E, nu, W, L);        
        
        AnalSolPointsX(iEval) = displacement(1);
        AnalSolPointsY(iEval) = displacement(2);

        if (corient == 1) || (corient == 3)
            PhysPoints(iEval) = coord(1);
        else
            PhysPoints(iEval) = coord(2);
        end

       
        for j=1:nument
            globnum = element_nod(elementind,j);     
          
            locnum = nument+1-j;
            cR = R(locnum);
            
            xcontrolpt = coordinates(globnum,1);
            ycontrolpt = coordinates(globnum,2);            
                                    
            globindx = 2*globnum-1;
            globindy = 2*globnum;
                        
            PlotPointsX(iEval) = PlotPointsX(iEval) + cR.*u(globindx);
            PlotPointsY(iEval) = PlotPointsY(iEval) + cR.*u(globindy);                                             
        end                
    end
    plot(PhysPoints, PlotPointsX, '-b', PhysPoints, AnalSolPointsX, '-g');
    
    hold on
    plot(PhysPoints, PlotPointsY, '-b', PhysPoints, AnalSolPointsY, '-g');
    
    hold on    
end
title('Computed and analytical displacements on the Dirichlet boundary')
