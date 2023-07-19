clc; clear;close all
syms m1 mb m2 l k c
syms x t a r r_f dx dt da dr g
q=[x t a r].'
dq=[dx dt da dr].'
syms Ix_b Iy_b Iz_b

I_b=[Ix_b    0    0
        0 Iy_b    0
        0    0 Iz_b]%xyz

w_1=[dt;0;0]%XYZ
w_2=[0;0;da]%xyz

R_1=[1       0      0
    0  cos(t) sin(t)
    0 -sin(t) cos(t)]

R_2=[cos(a)  sin(a) 0
    -sin(a)  cos(a) 0
          0       0 1]
R=R_2*R_1
w=R*w_1+w_2%xyz
H_G=I_b*w%xyz

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
vf=jacobian(rf,q)*dq%XYZ
v2=jacobian(r2,q)*dq%XYZ

e1=simplify(1/2*m1*v1.'*v1)
e2=simplify(1/2*m2*v2.'*v2)
eb=simplify(1/2*mb*vb.'*vb+1/2*w.'*H_G)



T=1/2*(m1*v1.'*v1  +m2*v2.'*v2  +mb*vb.'*vb  +w.'*H_G)
V=m2*g*r2(2)+mb*g*rb(2)+1/2*k*(r-r_f)^2
L=T-V

syms To1 To2 f alpha
F=f*[-cos(alpha);-sin(alpha);0]%xyz
Po_i=-c*dr^2+F.'*(R*vf)-To2*dt-To1*da
Q=simplify(jacobian(Po_i,dq).')
Q(end)=Q(end)/2

syms ddx ddt dda ddr
ddq=[ddx ddt dda ddr].'
s=[q;dq];ds=[dq;ddq];
P=jacobian(L,dq)
term_1=jacobian(P,q)*dq+jacobian(P,dq)*ddq
term_2=-jacobian(L,q).'
eqn=term_1+term_2-Q

%eqn=M*ddq+A=0
M=jacobian(P,dq)
A=subs(eqn,ddq,[0;0;0;0])
M=simplify(M)
A=simplify(A)
%ddq=-M\A
E=simplify(T+V)
T
V
P_x=m1*v1(1)+m2*v2(1)+mb*vb(1)