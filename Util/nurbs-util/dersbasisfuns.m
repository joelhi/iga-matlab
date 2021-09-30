%  This Function evaluates the basis functions and first derivatives 
%               functions at a given parameter value u.
% 
%  Algorithm from Piegl, Les. "The NURBS Book". Springer-Verlag: 
%     Berlin 1995; pp. 72-73.
% 
%  June 17, 2003
%  J. Austin Cottrell
%  CES Graduate Student
%  Texas Institute for Computational Engineering Science
%  University of Texas at Austin
%
%  Modified to codes Matlab by :
%  Chien Thai Hoang & Hung Nguyen Xuan
%
%   Faculty of Mathematics & Informatics, University of Natural Sciences
%   Vietnam   National University–HCM


function [ders]=dersbasisfuns(i,pl,ml,u,u_knotl)

%     --------------variable declarations--------------------------------
% Input: i,pl,ml   %knot span, degree of curve, number of control points
%        u , u_knotl  %parameter value, vector of knots
% Out:   ders(1,:)  % shape matrix
%        ders(2,:)  % derivative matrix

left=zeros(1,pl+1);
right=zeros(1,pl+1);
ndu=zeros(pl+1,pl+1);
a=zeros(2,pl+1);
ders=zeros(2,pl+1);

%     -------------------------------------------------------------------

ndu(1,1) = 1;
for j = 1:pl
left(j+1) = u - u_knotl(i+1-j);
right(j+1) = u_knotl(i+j) - u;
saved = 0;
   for r = 0:j-1
       ndu(j+1,r+1) = right(r+2) + left(j-r+1);
       temp = ndu(r+1,j)/ndu(j+1,r+1);%%%%% error here
       ndu(r+1,j+1) = saved + right(r+2)*temp;
       saved = left(j-r+1)*temp;
   end 
ndu(j+1,j+1) = saved;
end 

% load basis functions
for j = 0:pl
ders(1,j+1) = ndu(j+1,pl+1);
end 

% compute derivatives
for r = 0:pl % loop over function index
s1 = 0;
% alternate rows in array a
s2 = 1;
a(1,1) = 1;

   for k = 1:1; % loop to compute kth derivative
       d = 0d+0;
       rk = fix(r-k);
       pk = fix(pl-k);
       if(r >= k)
          a(s2+1,1) = a(s1+1,1)./ndu(pk+2,rk+1);
          d = a(s2+1,1).*ndu(rk+1,pk+1);
       end
       if(rk >= -1)
          j1 = 1;
       else
          j1 = fix(-rk);
       end
       if((r-1) <= pk)
          j2 = fix(k-1);
       else
          j2 = fix(pl-r);
       end
    for j = j1:j2;
        a(s2+1,j+1) =(a(s1+1,j+1) - a(s1+1,j))./ndu(pk+2,rk+j+1);
        d = d + a(s2+1,j+1).*ndu(rk+j+1,pk+1);
    end % j =fix(j2+1);
    if(r <= pk)
       a(s2+1,k+1) = -a(s1+1,k)./ndu(pk+2,r+1);
       d = d + a(s2+1,k+1).*ndu(r+1,pk+1);
    end
ders(k+1,r+1) = d;
j = fix(s1);
s1 = fix(s2);
% switch rows
s2 = fix(j);
   end 
end    

%   multiply through by the correct factors;
r = fix(pl);
for k = 1:1;
   for j = 0:pl;
       ders(k+1,j+1) = ders(k+1,j+1).*r;
   end
r = fix(r.*(pl-k));
end

return
