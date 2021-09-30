%compute l2-error and energy norm

l2errvar = 0; %l2 norm of the error
l2norm = 0; %l2 norm of the solution
energynorm = 0; %energy norm of the computed solution
energyerrnorm = 0; %energy norm of the error
invC = inv(C);

for i=1:elementcounter
    localindex = 0;
    
    ximin = element_int(i,1);
    etamin = element_int(i,2);
    ximax = element_int(i,3);
    etamax = element_int(i,4);            

    scalefac = (ximax - ximin)*(etamax - etamin)/4;
    [gpx,gwx]=genGP_GW(ngaussx);
    [gpy,gwy]=genGP_GW(ngaussy);

    temp = 0;

    dispmatx = zeros(ngaussy,ngaussx);
    dispmaty = zeros(ngaussy,ngaussx);
    exdispmatx = zeros(ngaussy, ngaussx);
    exdispmaty = zeros(ngaussy, ngaussx);

    %for each shape function in shapelist compute the matrices B, E
    for ii=1:ngaussx
        for jj=1:ngaussy
            [R, dRdx, coord, J] = nurbshaped(i, gpx(ii), gpy(jj), knotu, knotv, b_net, p, q, lenu, lenv, element_nod, coord_ij);                 
            B = zeros(localindex, 3);
            physx = coord(1);
            physy = coord(2);
            [eldisp, estress0] = gbeam(physx, physy, P, E, nu, W, L);                       
            exdispmatx(jj, ii) = eldisp(1);
            exdispmaty(jj, ii) = eldisp(2);
            scrtx = zeros(dim*nument, 1);
            for j=1:nument

                globnum = element_nod(i,j);                                
                cR = R(j);
                cdRdx = dRdx(1,j);
                cdRdy = dRdx(2,j);
                
                B(2*j-1, :) = [cdRdx, 0, cdRdy];
                B(2*j, :) = [0, cdRdy, cdRdx];
                
                globindx = 2*globnum-1;
                globindy = 2*globnum;
                
                scrtx(2*j-1) = 2*globnum-1;
                scrtx(2*j) = 2*globnum;

                dispmatx(jj,ii) = dispmatx(jj,ii) + cR*u(globindx);
                dispmaty(jj,ii) = dispmaty(jj,ii) + cR*u(globindy);                                                                        

            end

            
            
            estress = C*B'*u(scrtx);

            energyerrnorm = energyerrnorm + (estress0-estress)'*invC*(estress0-estress)*J.*scalefac.*gwx(ii).*gwy(jj);
            energynorm = energynorm + estress'*invC*estress.*J.*scalefac.*gwx(ii).*gwy(jj);

            l2errvar = l2errvar + ((exdispmatx(jj,ii)-dispmatx(jj,ii))^2+(exdispmaty(jj,ii)-dispmaty(jj,ii))^2).*J.*scalefac.*gwx(ii).*gwy(jj);
            l2norm = l2norm + ((exdispmatx(jj,ii))^2+(exdispmaty(jj,ii))^2).*J.*scalefac.*gwx(ii).*gwy(jj);
        end
    end
    
end

l2errvar = sqrt(l2errvar);
l2errrel = l2errvar/sqrt(l2norm)
energyerrnorm = sqrt(0.5*energyerrnorm);
energynorm = sqrt(0.5*energynorm);
energyerrrel = energyerrnorm/energynorm