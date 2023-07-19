function dz=ntn_fun(t,z)
global m1 mb m2 l k c ro r_f g f To1 To2 alpha Ix_b Iy_b Iz_b
dz=z;
x=z(1); t=z(2);a=z(3); r=z(4);  dx=z(5); dt=z(6); da=z(7); dr=z(8);

M =[                       m1 + m2 + mb,                                                                                                                             0,       -(sin(a)*(l*mb + 2*m2*r))/2, m2*cos(a)
                                  0, Ix_b/2 + Iy_b/2 + (l^2*mb)/8 + (m2*r^2)/2 + (Ix_b*cos(2*a))/2 - (Iy_b*cos(2*a))/2 - (l^2*mb*cos(2*a))/8 - (m2*r^2*cos(2*a))/2,                                 0,         0
 -(sin(a)*(2*m2*r + l*mb*cos(t)))/2,                                                                                              -(l^2*mb*cos(a)*sin(a)*sin(t))/4, (mb*cos(t)*l^2)/4 + m2*r^2 + Iz_b,         0
                         -m2*cos(a),                                                                                                                             0,                                 0,       -m2];
 
 
A =[                                                                                                                                                              f*cos(a + alpha) - 2*da*dr*m2*sin(a) - (da^2*l*mb*cos(a))/2 - da^2*m2*r*cos(a)
                                                     To2 + dr*dt*m2*r - Ix_b*da*dt*sin(2*a) + Iy_b*da*dt*sin(2*a) - dr*dt*m2*r*cos(2*a) - (g*l*mb*sin(a)*sin(t))/2 - g*m2*r*sin(a)*sin(t) + (da*dt*l^2*mb*sin(2*a))/4 + da*dt*m2*r^2*sin(2*a)
 To1 + f*l*sin(alpha) + (Ix_b*dt^2*sin(2*a))/2 - (Iy_b*dt^2*sin(2*a))/2 + 2*da*dr*m2*r - (dt^2*m2*r^2*sin(2*a))/2 + (g*l*mb*cos(a)*cos(t))/2 + g*m2*r*cos(a)*cos(t) - (da*dt*l^2*mb*cos(a)^2*sin(t))/2 - (dt^2*l^2*mb*cos(a)*sin(a)*cos(t))/4
                                                                                                                                                         k*r_f - k*r - c*dr + da^2*m2*r + dt^2*m2*r - dt^2*m2*r*cos(a)^2 - g*m2*sin(a)*cos(t)];
 
dz(1:4)=z(5:8);
dz(5:8)=-M\A;
% s=rank(M)
end