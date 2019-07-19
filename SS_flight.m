close all
clear all
clc

format shortG

%%
%initial conditions
V = 35; %m/s airspeed
alpha = 2.1131e-001; %35
beta = -2.0667e-002; %35
z = 2000*0.3048; %35
x = 0; y = 0;

u = V*cos(alpha)*cos(beta);
v = V*sin(beta);
w = V*sin(alpha)*cos(beta);
phi = 0; theta = 1.9190e-001; psi = 0; %35
q = 0; r = 0; p = 0;

IC = [p; q; r; phi; theta; psi; u; v; w; x; y; z];
%%

tend = 100;

[t,sol] = ode45(@(t,sol) odefun(t,sol), [0 tend], IC);


%%
figure(1)
plot3(sol(:,10), sol(:,11), sol(:,12), 'k', 'LineWidth', 2)
xlabel('x'); ylabel('y'); zlabel('z');
title('position')
grid on

figure(2)
subplot(311)
plot(t, sol(:,4).*180/pi, 'k', 'LineWidth', 1.5)
title('phi (degree)'); xlabel('time'); grid on
subplot(312)
plot(t, sol(:,5).*180/pi, 'k', 'LineWidth', 1.5)
title('theta (degree)'); xlabel('time'); grid on
subplot(313)
plot(t, sol(:,6).*180/pi, 'k', 'LineWidth', 1.5)
title('psi (degree)'); xlabel('time'); grid on

figure(3)
subplot(211)
plot(t, atan2(sol(:,9),sol(:,7)).*180./pi, 'k', 'LineWidth', 1.5)
title('Alpha (degree)'); xlabel('time'); grid on
subplot(212)
V_ = sqrt(sol(:,7).^2 + sol(:,8).^2 + sol(:,9).^2);
plot(t, asin(sol(:,8)./V_).*180./pi, 'k', 'LineWidth', 1.5)
title('Beta (degree)'); xlabel('time'); grid on

figure(4)
subplot(311)
plot(t, sol(:,1), 'k', 'LineWidth', 1.5)
title('p'); xlabel('time'); grid on
subplot(312)
plot(t, sol(:,2), 'k', 'LineWidth', 1.5)
title('q'); xlabel('time'); grid on
subplot(313)
plot(t, sol(:,3), 'k', 'LineWidth', 1.5)
title('r'); xlabel('time'); grid on

figure(5)
subplot(311)
plot(t, sol(:,7), 'k', 'LineWidth', 1.5)
title('u'); xlabel('time'); grid on
subplot(312)
plot(t, sol(:,8), 'k', 'LineWidth', 1.5)
title('v'); xlabel('time'); grid on
subplot(313)
plot(t, sol(:,9), 'k', 'LineWidth', 1.5)
title('w'); xlabel('time'); grid on

figure(6)
plot(t, sqrt(sol(:,7).^2 + sol(:,8).^2 + sol(:,9).^2), 'k', 'LineWidth', 1.5)
title('V_{total}'); xlabel('time'); grid on

figure(7)
plot(t, sol(:,12), 'k', 'LineWidth', 1.5)
title('H'); xlabel('time'); grid on





