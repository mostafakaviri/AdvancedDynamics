function dz=lag_fun(t,z)
global m1 mb m2 l k c  r_f g f To1 To2 alpha Ix_b Iy_b Iz_b
dz=z;
x=z(1); t=z(2);a=z(3); r=z(4);  dx=z(5); dt=z(6); da=z(7); dr=z(8);

M =[                m1 + m2 + mb,                                                                            0, -(sin(a)*(l*mb + 2*m2*r))/2, m2*cos(a)
                           0, Ix_b - Ix_b*sin(a)^2 + Iy_b*sin(a)^2 + (l^2*mb*sin(a)^2)/4 + m2*r^2*sin(a)^2,                           0,         0
 -(sin(a)*(l*mb + 2*m2*r))/2,                                                                            0,  (mb*l^2)/4 + m2*r^2 + Iz_b,         0
                   m2*cos(a),                                                                            0,                           0,        m2];
 
 
A =[                                                                                                               f*cos(a + alpha) - 2*da*dr*m2*sin(a) - (da^2*l*mb*cos(a))/2 - da^2*m2*r*cos(a)
  To2 + 2*dr*dt*m2*r - Ix_b*da*dt*sin(2*a) + Iy_b*da*dt*sin(2*a) - 2*dr*dt*m2*r*cos(a)^2 - (g*l*mb*sin(a)*sin(t))/2 - g*m2*r*sin(a)*sin(t) + (da*dt*l^2*mb*sin(2*a))/4 + da*dt*m2*r^2*sin(2*a)
 To1 + f*l*sin(alpha) + (Ix_b*dt^2*sin(2*a))/2 - (Iy_b*dt^2*sin(2*a))/2 + 2*da*dr*m2*r - (dt^2*l^2*mb*sin(2*a))/8 - (dt^2*m2*r^2*sin(2*a))/2 + (g*l*mb*cos(a)*cos(t))/2 + g*m2*r*cos(a)*cos(t)
                                                                                                          c*dr + k*r - k*r_f - da^2*m2*r - dt^2*m2*r + dt^2*m2*r*cos(a)^2 + g*m2*sin(a)*cos(t)];
 
 
dz(1:4)=z(5:8);
dz(5:8)=-M\A;

end