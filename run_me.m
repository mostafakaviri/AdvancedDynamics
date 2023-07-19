clc; clear;close all
global m1 mb m2 l k c ro r_f g f To1 To2 alpha Ix_b Iy_b Iz_b
% r_f=r of free for stifness
m1=1;mb=1;m2=1;l=1.2;k=10;ro=1;r_f=1;g=9.81;alpha=pi/4;
Ix_b=0;Iy_b=1/12*mb*l^2;Iz_b=Iy_b;
c=0;f=2; To1=0; To2=0;
x0=0
t0=0
a0=pi/2
r0=1
dx0=0
dt0=0
da0=0
dr0=0





z0=[x0 t0 a0 r0 dx0 dt0 da0 dr0];
t_end=2;
disp('lagrange(short rung)=1  Newton(short rung)=2   lagrange(ode45)=3   Newton(ode45)=4');
ord=input('ord=');

if ord==1;
    h=0.001;z=z0;i=1;
    for t_=h:h:t_end
        i=i+1;
        k1=h*lag_fun(t_,z(i-1,:));
        k2=h*lag_fun(t_+h/2,z(i-1,:)+k1/2);
        k3=h*lag_fun(t_+h/2,z(i-1,:)+k2/2);
        k4=h*lag_fun(t_+h,z(i-1,:)+k3);
        z(i,:)=z(i-1,:)+1/6*(k1+2*k2+2*k3+k4);
    end
    time=(0:h:t_end);
elseif ord==2
    h=0.001;z=z0;i=1;
    for t_=h:h:t_end
        i=i+1;
        k1=h*ntn_fun(t_,z(i-1,:));
        k2=h*ntn_fun(t_+h/2,z(i-1,:)+k1/2);
        k3=h*ntn_fun(t_+h/2,z(i-1,:)+k2/2);
        k4=h*ntn_fun(t_+h,z(i-1,:)+k3);
        z(i,:)=z(i-1,:)+1/6*(k1+2*k2+2*k3+k4);
    end
    time=(0:h:t_end);
    
elseif ord==3
    tspan=[0 t_end];
    [time,z]=ode45(@lag_fun,tspan,z0);
elseif ord==4
    tspan=[0 t_end];
    [time,z]=ode45(@ntn_fun,tspan,z0);
end
figure(1)
subplot(2,4,1);plot(time,z(:,1));xlabel('time');ylabel('x');grid on;hold on
subplot(2,4,2);plot(time,z(:,2));xlabel('time');ylabel('theta');grid on;hold on
subplot(2,4,3);plot(time,z(:,3));xlabel('time');ylabel('alpha');grid on;hold on
subplot(2,4,4);plot(time,z(:,4));xlabel('time');ylabel('r');grid on;hold on
subplot(2,4,5);plot(time,z(:,5));xlabel('time');ylabel('dx');grid on;hold on
subplot(2,4,6);plot(time,z(:,6));xlabel('time');ylabel('dtheta');grid on;hold on
subplot(2,4,7);plot(time,z(:,7));xlabel('time');ylabel('dalpha');grid on;hold on
subplot(2,4,8);plot(time,z(:,8));xlabel('time');ylabel('dr');grid on;hold on


for i=1:length(time)
    x=z(i,1); t=z(i,2);a=z(i,3); r=z(i,4);  dx=z(i,5); dt=z(i,6); da=z(i,7); dr=z(i,8);
    E(i) =(Iz_b*da^2)/2 + (Iy_b*dt^2)/2 + (dr^2*m2)/2 + (dx^2*m1)/2 + (dx^2*m2)/2 + (dx^2*mb)/2 + (k*r^2)/2 + (k*r_f^2)/2 + (da^2*l^2*mb)/8 + (dt^2*l^2*mb)/8 + (da^2*m2*r^2)/2 + (dt^2*m2*r^2)/2 + (Ix_b*dt^2*cos(a)^2)/2 - (Iy_b*dt^2*cos(a)^2)/2 - k*r*r_f + dr*dx*m2*cos(a) - (dt^2*l^2*mb*cos(a)^2)/8 - (dt^2*m2*r^2*cos(a)^2)/2 + (g*l*mb*sin(a)*cos(t))/2 + g*m2*r*sin(a)*cos(t) - (da*dx*l*mb*sin(a))/2 - da*dx*m2*r*sin(a);
    T(i) =(mb*((da*l*cos(a)*cos(t))/2 - (dt*l*sin(a)*sin(t))/2)^2)/2 + (mb*((da*l*cos(a)*sin(t))/2 + (dt*l*sin(a)*cos(t))/2)^2)/2 + (Iz_b*da^2)/2 + (m2*(dx + dr*cos(a) - da*r*sin(a))^2)/2 + (dx^2*m1)/2 + (m2*(dr*sin(a)*cos(t) + da*r*cos(a)*cos(t) - dt*r*sin(a)*sin(t))^2)/2 + (m2*(dr*sin(a)*sin(t) + da*r*cos(a)*sin(t) + dt*r*sin(a)*cos(t))^2)/2 + (mb*(dx - (da*l*sin(a))/2)^2)/2 + (Ix_b*dt^2*cos(a)^2)/2 + (Iy_b*dt^2*sin(a)^2)/2;
    V(i) =(k*(r - r_f)^2)/2 + (g*l*mb*sin(a)*cos(t))/2 + g*m2*r*sin(a)*cos(t);
    P_x(i) =m2*(dx + dr*cos(a) - da*r*sin(a)) + dx*m1 + mb*(dx - (da*l*sin(a))/2);
end
erroe_E=max(E)-min(E)
erroe_P_x=max(P_x)-min(P_x)
figure(2)
subplot(2,1,1);
plot(time,P_x);xlabel('time');ylabel('P_x');grid on; hold on
subplot(2,1,2);
plot(time,T,time,V,time,E);xlabel('time');ylabel('Energy');grid on; legend('kinetic Energy','Potential Energy','Energy');hold on
