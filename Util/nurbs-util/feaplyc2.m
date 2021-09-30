function [kk,ff]=feaplyc2(kk,ff,bcdof,bcvall)

%----------------------------------------------------------
%  Purpose:
%     Apply constraints to matrix equation [kk]{x}={ff}
%
%  Synopsis:
%     [kk,ff]=feaplybc(kk,ff,bcdof,bcval)
%
%  Variable Description:
%     kk - system matrix before applying constraints 
%     ff - system vector before applying constraints
%     bcdof - a vector containging constrained d.o.f
%     bcval - a vector containing contained value 
%
%     For example, there are constraints at d.o.f=2 and 10
%     and their constrained values are 0.0 and 2.5, 
%     respectively.  Then, bcdof(1)=2 and bcdof(2)=10; and
%     bcval(1)=1.0 and bcval(2)=2.5.
%-----------------------------------------------------------
 
 n=length(bcdof);
 sdof=length(kk);

 for i=1:n
    c=bcdof(i);
    for j=1:sdof
       kk(c,j)=0;
    end
    kk(c,c)=1;
    ff(c,1)=bcvall(i);
 end

