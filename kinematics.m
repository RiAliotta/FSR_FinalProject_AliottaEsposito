%% Kinematics

clc

q = [0; 0];
l = 0.6;
m = 3.4;
m_link = m/4;

C1 = [l/2; 0; 0];
C2 = [l+l/2*cos(q(1)); l/2*sin(q(1)); 0];
C3 = [l+l*cos(q(1))+l/2*cos(q(1)+pi/2); l*sin(q(1))+l/2*sin(q(1)+pi/2); 0];
C4 = [l+l*cos(q(1))+l*cos(q(1)+pi/2)+l/2*cos(q(1)+q(2)+pi/2); l*sin(q(1))+l*sin(q(1)+pi/2)+l/2*sin(q(1)+q(2)+pi/2); 0];
C = (C1+C2+C3+C4)/4;
r1 = C1-C;
r2 = C2-C;
r3 = C3-C;
r4 = C4-C;
I_link1 = [0 0 0;
        0 1/12*m_link*l^2 0;
        0 0 1/12*m_link*l^2]
I_link2 = [cos(q(1)) -sin(q(1)) 0;
    sin(q(1)) cos(q(1)) 0;
    0 0 1]*I_link1*[cos(q(1)) -sin(q(1)) 0;
    sin(q(1)) cos(q(1)) 0;
    0 0 1]'
I_link3 = [cos(q(1)+pi/2) -sin(q(1)+pi/2) 0;
    sin(q(1)+pi/2) cos(q(1)+pi/2) 0;
    0 0 1]*I_link1*[cos(q(1)+pi/2) -sin(q(1)+pi/2) 0;
    sin(q(1)+pi/2) cos(q(1)+pi/2) 0;
    0 0 1]'
I_link4 = [cos(q(1)+pi/2+q(2)) -sin(q(1)+pi/2+q(2)) 0;
    sin(q(1)+pi/2+q(2)) cos(q(1)+pi/2+q(2)) 0;
    0 0 1]*I_link1*[cos(q(1)+pi/2+q(2)) -sin(q(1)+pi/2+q(2)) 0;
    sin(q(1)+pi/2+q(2)) cos(q(1)+pi/2+q(2)) 0;
    0 0 1]'
I = I_link1+I_link2+I_link3+I_link4+m_link*(r1'*r1*eye(3)-r1*r1'+r2'*r2*eye(3)-r2*r2'+r3'*r3*eye(3)-r3*r3'+r4'*r4*eye(3)-r4*r4')

figure(1)
axis square
axis([C(1)-1.2 C(1)+1.2 C(2)-1.2 C(2)+1.2])
hold on
xlabel("Drone configuration seen from above")
plot([0 0.6], [0 0], 'r-')
plot([0.6 0.6+l*cos(q(1))], [0 l*sin(q(1))], 'b-');
plot([0.6+l*cos(q(1)) 0.6+l*cos(q(1))+l*cos(q(1)+pi/2)], [l*sin(q(1)) l*sin(q(1))+l*sin(q(1)+pi/2)], 'r-');
plot([0.6+l*cos(q(1))+l*cos(q(1)+pi/2) 0.6+l*cos(q(1))+l*cos(q(1)+pi/2)+l*cos(q(1)+pi/2+q(2))], [l*sin(q(1))+l*sin(q(1)+pi/2) l*sin(q(1))+l*sin(q(1)+pi/2)+l*sin(q(1)+pi/2+q(2))], 'b-');
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
