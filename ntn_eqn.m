clc; clear;close all
syms m1 mb m2 l k c 

syms x t a r r_f dx dt da dr ddx ddt dda ddr g
q=[x t a r].'
dq=[dx dt da dr].'
ddq=[ddx ddt dda ddr].'
s=[q;dq];ds=[dq;ddq];

syms Ix_b Iy_b Iz_b

I_b=[Ix_b    0    0
        0 Iy_b    0
        0    0 Iz_b]%xyz

w_1=[dt;0;0]%XYZ
w_2=[0;0;da]%xyz

R_1=[1       0      0
     0  cos(t)  sin(t)
     0 -sin(t)  cos(t)]

R_2=[cos(a)  sin(a) 0
    -sin(a)  cos(a) 0
          0       0 1]
R=R_2*R_1
w=R*w_1+w_2%xyz
H_G=I_b*w%xyz
H_G=R.'*H_G%XYZ

r1=[x;0;0]%XYZ

rb_1=[l/2;0;0]%xyz
rf_1=[l;0;0]%xyz
r2_1=[r;0;0]%xyz

rb=R*r1+rb_1%xyz
rf=R*r1+rf_1%xyz
r2=R*r1+r2_1%xyz

rb=simplify(R.'*rb)%XYZ
rf=simplify(R.'*rf)%XYZ
r2=simplify(R.'*r2)%XYZ

v1=jacobian(r1,q)*dq%XYZ
vb=jacobian(rb,q)*dq%XYZ
v2=jacobian(r2,q)*dq%XYZ


a1=jacobian(v1,s)*ds%XYZ
ab=jacobian(vb,s)*ds%XYZ
a2=jacobian(v2,s)*ds%XYZ


syms To1 To2 f alpha
F=f*[-cos(alpha);-sin(alpha);0]%xyz
F=R.'*F%XYZ
eqn(1)=-F(1)+m1*a1(1)+m2*a2(1)+mb*ab(1)       %F_x=ma_x
Moment=cross(rb-r1,[0;-mb*g;0])+cross(r2-r1,[0;-m2*g;0])-[To2;0;0]+cross(rf-r1,F)%XYZ
dH_G=jacobian(H_G,s)*ds%XYZ
eqn(2)=-Moment(1)   +cross(r2-r1,m2*a2).'*[1;0;0]  +dH_G(1)    +cross(rb-r1,mb*ab).'*[1;0;0]

Moment=R*Moment%xyz
eqn(3)=-Moment(3) +To1 +(R*cross(r2-r1,m2*a2)).'*[0;0;1] +(R*dH_G).'*[0;0;1]+(cross(rb-r1,mb*ab)).'*[0;0;1] 

eqn(4)=(R*[0;-m2*g;0]).'*[1;0;0]-k*(r-r_f)-c*dr-(R*m2*a2).'*[1;0;0]
eqn=simplify(eqn).'



%eqn=M*ddq+A=0
M=jacobian(eqn,ddq)
A=subs(eqn,ddq,[0;0;0;0])
M=simplify(M)
A=simplify(A)
%ddq=-M\A

