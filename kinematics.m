%% Inertia computation

clc

q = [pi/2; pi/2];
q1 = q(1);
q2 = q(2);
l = 0.6;
m = 3.4;

% syms q1 q2 l m real

m_link = m/4;

C1 = [l/2; 0; 0];
C2 = [l+l/2*cos(q1); l/2*sin(q1); 0];
C3 = [l+l*cos(q1)+l/2*cos(q1+pi/2); l*sin(q1)+l/2*sin(q1+pi/2); 0];
C4 = [l+l*cos(q1)+l*cos(q1+pi/2)+l/2*cos(q1+q2+pi/2); l*sin(q1)+l*sin(q1+pi/2)+l/2*sin(q1+q2+pi/2); 0];
C = (C1+C2+C3+C4)/4;
r1 = C1-C;
r2 = C2-C;
r3 = C3-C;
r4 = C4-C;
R_12 = [cos(q1) -sin(q1) 0;
        sin(q1) cos(q1)  0;
        0       0        1];
R_23 = [cos(q1+pi/2) -sin(q1+pi/2) 0;
        sin(q1+pi/2) cos(q1+pi/2)  0;
        0            0             1];
R_34 = [cos(q1+pi/2+q2) -sin(q1+pi/2+q2) 0;
        sin(q1+pi/2+q2) cos(q1+pi/2+q2)  0;
        0               0                1];
I_link1 = [0 0               0;
           0 1/12*m_link*l^2 0;
           0 0               1/12*m_link*l^2];
I_link2 = R_12*I_link1*R_12';
I_link3 = R_23*I_link1*R_23';
I_link4 = R_34*I_link1*R_34';
I = I_link1+I_link2+I_link3+I_link4+m_link*(r1'*r1*eye(3)-r1*r1'+r2'*r2*eye(3)-r2*r2'+r3'*r3*eye(3)-r3*r3'+r4'*r4*eye(3)-r4*r4')

%% Allocation matrix computation

beta = 10*2*pi/360; 
ki = 0.01;
% syms beta ki real

R_F1C = [cos(-beta) 0 sin(-beta);
         0          1 0;
        -sin(-beta) 0 cos(-beta)];

R_F2C = [cos(q1) -sin(q1) 0;
         sin(q1) cos(q1) 0;
         0 0 1]*R_F1C;

R_F3C = [cos(q1+pi/2) -sin(q1+pi/2) 0;
    sin(q1+pi/2) cos(q1+pi/2) 0;
    0 0 1]*R_F1C;

R_F4C = [cos(q1+pi/2+q2) -sin(q1+pi/2+q2) 0;
    sin(q1+pi/2+q2) cos(q1+pi/2+q2) 0;
    0 0 1]*R_F1C;

Qp_tran = [R_F1C*[0; 0; 1], R_F2C*[0; 0; 1], R_F3C*[0; 0; 1], R_F4C*[0; 0; 1]];

Qp_rot = [(skew(C1-C)+ki*eye(3))*R_F1C*[0; 0; 1], (skew(C2-C)+ki*eye(3))*R_F2C*[0; 0; 1], (skew(C3-C)+ki*eye(3))*R_F3C*[0; 0; 1], (skew(C4-C)+ki*eye(3))*R_F4C*[0; 0; 1]];

up = ([Qp_tran(3,:); Qp_rot])\[1; 0; 0; 0];
us = m*9.81/norm(Qp_tran*up)*up;
fs = Qp_tran*us;
alfa_x = atan2(fs(2), fs(3));
alfa_y = atan2(-fs(1), sqrt(fs(2)^2+fs(3)^2));
R_C_CoG = [cos(alfa_y)  0 sin(alfa_y);
           0            1 0;
           -sin(alfa_y) 0 cos(alfa_y)]*...
           [1 0           0;
            0 cos(alfa_x) -sin(alfa_x);
            0 sin(alfa_x) cos(alfa_x)];
Q_tran = R_C_CoG*Qp_tran;
Q_rot = R_C_CoG*Qp_rot;

%% Plot drone configuration

figure(1)
axis square
axis([C(1)-1.2 C(1)+1.2 C(2)-1.2 C(2)+1.2])
hold on
xlabel("Drone configuration seen from above")
plot([0 0.6], [0 0], 'r-')
plot([0.6 0.6+l*cos(q1)], [0 l*sin(q1)], 'b-');
plot([0.6+l*cos(q1) 0.6+l*cos(q1)+l*cos(q1+pi/2)], [l*sin(q1) l*sin(q1)+l*sin(q1+pi/2)], 'r-');
plot([0.6+l*cos(q1)+l*cos(q1+pi/2) 0.6+l*cos(q1)+l*cos(q1+pi/2)+l*cos(q1+pi/2+q2)], [l*sin(q1)+l*sin(q1+pi/2) l*sin(q1)+l*sin(q1+pi/2)+l*sin(q1+pi/2+q2)], 'b-');
plot(C(1), C(2), 'kx')
hold off


function S = skew(v)
    if(numel(v)~= 1)
        S= [0 -v(3) v(2); 
            v(3) 0 -v(1);
            -v(2) v(1) 0];
    else
        S= zeros(3);
    end
end