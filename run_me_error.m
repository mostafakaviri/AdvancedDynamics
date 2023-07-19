clc; clear;close all
global m1 mb m2 l k c ro r_f g f To1 To2 alpha Ix_b Iy_b Iz_b
% r_f=r of free for stifness
m1=1;mb=1;m2=1;l=1;k=10;ro=0;r_f=1;g=9.81;alpha=pi/4;
Ix_b=0;Iy_b=1/12*mb*l^2;Iz_b=Iy_b;
c=.5;f=2; To1=0; To2=0;
x0=0
t0=0
a0=pi/2
r0=.75
dx0=0
dt0=0
da0=0
dr0=0


z0=[x0 t0 a0 r0 dx0 dt0 da0 dr0];
t_end=2;

for ss=1:100
    tic
    kk(ss)=(ss)/1000;
    h=kk(ss)
   z=z0;i=1;
    for t_=h:h:t_end
        i=i+1;
        k1=h*lag_fun(t_,z(i-1,:));
        k2=h*lag_fun(t_+h/2,z(i-1,:)+k1/2);
        k3=h*lag_fun(t_+h/2,z(i-1,:)+k2/2);
        k4=h*lag_fun(t_+h,z(i-1,:)+k3);
        z(i,:)=z(i-1,:)+1/6*(k1+2*k2+2*k3+k4);
    end
    time=(0:h:t_end);


for i=1:length(time)
    x=z(i,1); t=z(i,2);a=z(i,3); r=z(i,4);  dx=z(i,5); dt=z(i,6); da=z(i,7); dr=z(i,8);
    E(i) =(Iz_b*da^2)/2 + (Iy_b*dt^2)/2 + (dr^2*m2)/2 + (dx^2*m1)/2 + (dx^2*m2)/2 + (dx^2*mb)/2 + (k*r^2)/2 + (k*r_f^2)/2 + (da^2*l^2*mb)/8 + (dt^2*l^2*mb)/8 + (da^2*m2*r^2)/2 + (dt^2*m2*r^2)/2 + (Ix_b*dt^2*cos(a)^2)/2 - (Iy_b*dt^2*cos(a)^2)/2 - k*r*r_f + dr*dx*m2*cos(a) - (dt^2*l^2*mb*cos(a)^2)/8 - (dt^2*m2*r^2*cos(a)^2)/2 + (g*l*mb*sin(a)*cos(t))/2 + g*m2*r*sin(a)*cos(t) - (da*dx*l*mb*sin(a))/2 - da*dx*m2*r*sin(a);
    T(i) =(mb*((da*l*cos(a)*cos(t))/2 - (dt*l*sin(a)*sin(t))/2)^2)/2 + (mb*((da*l*cos(a)*sin(t))/2 + (dt*l*sin(a)*cos(t))/2)^2)/2 + (Iz_b*da^2)/2 + (m2*(dx + dr*cos(a) - da*r*sin(a))^2)/2 + (dx^2*m1)/2 + (m2*(dr*sin(a)*cos(t) + da*r*cos(a)*cos(t) - dt*r*sin(a)*sin(t))^2)/2 + (m2*(dr*sin(a)*sin(t) + da*r*cos(a)*sin(t) + dt*r*sin(a)*cos(t))^2)/2 + (mb*(dx - (da*l*sin(a))/2)^2)/2 + (Ix_b*dt^2*cos(a)^2)/2 + (Iy_b*dt^2*sin(a)^2)/2;
    V(i) =(k*(r - r_f)^2)/2 + (g*l*mb*sin(a)*cos(t))/2 + g*m2*r*sin(a)*cos(t);
    P_x(i) =m2*(dx + dr*cos(a) - da*r*sin(a)) + dx*m1 + mb*(dx - (da*l*sin(a))/2);
end
t_hal(ss)=toc;
erroe_E(ss)=max(E)-min(E);
erroe_P_x(ss)=max(P_x)-min(P_x);
end
figure
subplot(2,1,1)
plot(kk,erroe_E,kk,t_hal);xlabel('time step');legend('Error_E','sul time');grid on
subplot(2,1,2)
plot(kk,erroe_P_x,kk,t_hal);xlabel('time step');legend('Error_P','sul time');grid on