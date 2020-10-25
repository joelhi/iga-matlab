
clc, close all, clear all
e = [0.2 0;0 0]

a1 = [10 0];
a2 = [3 3];

b1 = [1 0];
b2 = [cos(0.5*pi) sin(0.5*pi)];

[e_trans] = shiftBase2D(e,a1,a2,b1,b2,1)
[e_trans2] = shiftBase2D(e,a1,a2,b1,b2,2)


[a1_c,a2_c] = computeContravariant(a1,a2);

nu = 0.1;
t = 1;
E = 334456;

memStiff = E*t/(1-nu^2);

C = [1 nu 0; nu 1 0;0 0 (1-nu)/2];

s_v = C * [e(1,1); e(2,2); e(1,2)+e(2,1)]; 

s = [s_v(1) 0.5*s_v(3);0.5*s_v(3) s_v(2)]

[s_trans1] = shiftBase2D(s,a1,a2,b1,b2,1)
[s_trans2] = shiftBase2D(s,a1,a2,b1,b2,2)

hold on
plot([0 a1(1)],[0 a1(2)],'k--o')
plot([0 a2(1)],[0 a2(2)],'k--o')
plot([0 b1(1)],[0 b1(2)],'r-x')
plot([0 b2(1)],[0 b2(2)],'r-x')
axis equal
hold off


edots1 = sum(sum(e.*s))
edots2 = sum(dot(e_trans,s_trans1))
edots3 = sum(dot(e_trans,s_trans2))