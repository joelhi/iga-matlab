%Displaying the displacements
element4 = zeros(elementcounter, 4);
physcoord = zeros(4*elementcounter, 2);
dispcoord = zeros(4*elementcounter, 2);
straincoord = zeros(4*elementcounter, 3);
sigmacoord = zeros(4*elementcounter, 4);

disp('Phys. Coord.         Analytic Sol.                   Computed Sol.')

maxerr = 0;
for i=1:elementcounter
    localindex = 0;    
    element4(i, :) = [(i-1)*4+1:(i-1)*4+4];
    %for each node in the element, count the number of shape functions
    
        
    nx = 2;
    ny = 2;
    px = linspace(-1, 1, nx);
    py = linspace(-1, 1, ny);
    R = cell(ny, nx);
    dRdx = cell(ny, nx);
    J = cell(ny, nx);
    coord = cell(ny,nx);
   

    %for each shape function in shapelist compute the matrices B, E
    dispmatx = zeros(ny,nx);
    dispmaty = zeros(ny,nx);
    exactdispx = zeros(ny,nx);
    exactdispy = zeros(ny, nx);
    strain11 = zeros(ny,nx);
    strain12 = zeros(ny,nx);
    strain22 = zeros(ny,nx);
    stressvect = cell(ny,nx);
    
    sigmad = zeros(ny,nx);
    
    for ii=1:nx
        for jj=1:ny                 
            
            [R, dRdx, coord{jj,ii}, J] = nurbshaped(i, px(ii), py(jj), knotu, knotv, b_net, p, q, lenu, lenv, element_nod, coord_ij); 
            [eldisp, estress0] = gbeam(coord{jj,ii}(1), coord{jj,ii}(2), P, E, nu, W, L);
            exactdispx(jj,ii) = eldisp(1);
            exactdispy(jj,ii) = eldisp(2);
            

            for j=1:nument
                globnum = element_nod(i,j);                                
                cR = R(j);
                cdRdx = dRdx(1,j);
                cdRdy = dRdx(2,j);

                globindx = 2*globnum-1;
                globindy = 2*globnum;

                cR = R(j);
                cdRdx = dRdx(1,j);
                cdRdy = dRdx(2,j);

                %calculate the value of enrichment and derivatives
                xcontrolpt = coordinates(globnum,1);
                ycontrolpt = coordinates(globnum,2);

                %calculate displacement values
                dispmatx(jj,ii) = dispmatx(jj,ii) + cR*u(globindx);
                dispmaty(jj,ii) = dispmaty(jj,ii) + cR*u(globindy);

                %calculate strain values
                strain11(jj,ii) = strain11(jj,ii) + (cdRdx)*u(globindx);
                strain12(jj,ii) = strain12(jj,ii) + 1/2*(cdRdy)*u(globindx) + 1/2*(cdRdx)*u(globindy);
                strain22(jj,ii) = strain22(jj,ii) + (cdRdy)*u(globindy);                                

                stressvect{jj,ii} = C*[strain11(jj,ii); strain22(jj,ii); strain12(jj,ii)];
                sigma11 = stressvect{jj,ii}(1);
                sigma12 = stressvect{jj,ii}(3);
                sigma22 = stressvect{jj,ii}(2);
               
            end
    
        end     
    end
    
    maxerr = max([maxerr, max(max(abs(dispmatx - exactdispx))), max(max(abs(dispmaty - exactdispy)))]);
    
    
    physcoord((i-1)*4+1, :) = coord{1,1};
    physcoord((i-1)*4+2, :) = coord{1,2};
    physcoord((i-1)*4+3, :) = coord{2,2};
    physcoord((i-1)*4+4, :) = coord{2,1};
    
    %print the displacements at the corners of the elemnt
    if ((i/numix)==ceil(i/numix)) && (i>(numix*(numiy-1))) %if in the top right corner display values at all four corners
        disp(sprintf('(%1.1f  %1.1f)   %1.10f   %1.10f    %1.10f   %1.10f',  coord{1,1}(1),  coord{1,1}(2),  exactdispx(1,1), exactdispy(1,1), dispmatx(1,1), dispmaty(1,1)));
        disp(sprintf('(%1.1f  %1.1f)   %1.10f   %1.10f    %1.10f   %1.10f',  coord{1,2}(1),  coord{1,2}(2),  exactdispx(1,2), exactdispy(1,2), dispmatx(1,2), dispmaty(1,2)));
        disp(sprintf('(%1.1f  %1.1f)   %1.10f   %1.10f    %1.10f   %1.10f',  coord{2,1}(1),  coord{2,1}(2),  exactdispx(2,1), exactdispy(2,1), dispmatx(2,1), dispmaty(2,1)));
        disp(sprintf('(%1.1f  %1.1f)   %1.10f   %1.10f    %1.10f   %1.10f',  coord{2,2}(1),  coord{2,2}(2),  exactdispx(2,2), exactdispy(2,2), dispmatx(2,2), dispmaty(2,2)));
    elseif ((i/numix)==ceil(i/numix)) %if on the right edge display values from the bottom two nodes
        disp(sprintf('(%1.1f  %1.1f)   %1.10f   %1.10f    %1.10f   %1.10f',  coord{1,1}(1),  coord{1,1}(2),  exactdispx(1,1), exactdispy(1,1), dispmatx(1,1), dispmaty(1,1)));
        disp(sprintf('(%1.1f  %1.1f)   %1.10f   %1.10f    %1.10f   %1.10f',  coord{1,2}(1),  coord{1,2}(2),  exactdispx(1,2), exactdispy(1,2), dispmatx(1,2), dispmaty(1,2)));
    elseif (i>(numix*(numiy-1))) %if on the top edge, display values on the left two nodes
        disp(sprintf('(%1.1f  %1.1f)   %1.10f   %1.10f    %1.10f   %1.10f',  coord{1,1}(1),  coord{1,1}(2),  exactdispx(1,1), exactdispy(1,1), dispmatx(1,1), dispmaty(1,1)));
        disp(sprintf('(%1.1f  %1.1f)   %1.10f   %1.10f    %1.10f   %1.10f',  coord{2,1}(1),  coord{2,1}(2),  exactdispx(2,1), exactdispy(2,1), dispmatx(2,1), dispmaty(2,1)));
    else
        disp(sprintf('(%1.1f  %1.1f)   %1.10f   %1.10f    %1.10f   %1.10f',  coord{1,1}(1),  coord{1,1}(2),  exactdispx(1,1), exactdispy(1,1), dispmatx(1,1), dispmaty(1,1)));
        
    end
    
    
    
    dispcoord((i-1)*4+1, :) = [dispmatx(1,1) dispmaty(1,1)];
    dispcoord((i-1)*4+2, :) = [dispmatx(1,2) dispmaty(1,2)];
    dispcoord((i-1)*4+3, :) = [dispmatx(2,2) dispmaty(2,2)];
    dispcoord((i-1)*4+4, :) = [dispmatx(2,1) dispmaty(2,1)];
    straincoord((i-1)*4+1, :) = [strain11(1,1) strain12(1,1) strain22(1,1)];
    straincoord((i-1)*4+2, :) = [strain11(1,2) strain12(1,2) strain22(1,2)];
    straincoord((i-1)*4+3, :) = [strain11(2,2) strain12(2,2) strain22(2,2)];
    straincoord((i-1)*4+4, :) = [strain11(2,1) strain12(2,1) strain22(2,1)];
      
    
    sigmacoord((i-1)*4+1, :) = [stressvect{1,1}' sigmad(1,1)];
    sigmacoord((i-1)*4+2, :) = [stressvect{1,2}' sigmad(1,2)];
    sigmacoord((i-1)*4+3, :) = [stressvect{2,2}' sigmad(2,2)];
    sigmacoord((i-1)*4+4, :) = [stressvect{2,1}' sigmad(2,1)];
    
end

maxerr

factor = 10;
%factor = 0;
%close all
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), straincoord(:,1), 'facecolor','interp')
view(0,90)
title('Displacements and \epsilon_{11}')
colorbar('vert')

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), straincoord(:,2), 'facecolor','interp')
view(0,90)
colorbar('vert')
title('Displacements and \epsilon_{12}')


figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), straincoord(:,3), 'facecolor','interp')
view(0,90)
title('Displacements and \epsilon_{22}')
colorbar('vert')

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,1), 'facecolor','interp')
view(0,90)
title('Displacements and \sigma_{11}')
colorbar('vert')

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,3), 'facecolor','interp')
view(0,90)
title('Displacements and \sigma_{12}')
colorbar('vert')

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,2), 'facecolor','interp')
view(0,90)
title('Displacements and \sigma_{22}')
colorbar('vert')




    