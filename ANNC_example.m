% 
% Author: Z-JC
% Data: 2022-10-16
clear
clc

%%
% states
x_1(:,1) = 2;
x_2(:,1) = 0.2;
h(5,:)=0;
fx(:,1) = 10 * x_1(:,1) * x_2(:,1);

% Control inputs
u(:,1) = 0;

% Desired 
x_d(:,1) = sin(0);
ddot_x_d(:,1) = cos(0);
dot_x_d(:,1) = -sin(0);

% Parameters
c = 1;
eta = 0.5;

% RBF
InpLayNum = 2;
HidLayNum = 5;
OutLayNum = 1;
gamma = 500;
eta = 0.5;
c_i = [-2 -1  0  1  2
       -2 -1  0  1  2];
b_i = 3;
hat_W = 0.1*ones(HidLayNum, OutLayNum);

%% Time state
t(1,1) = 0;
tBegin = 0;
tFinal = 20;
dT = 0.01;
times = (tFinal-tBegin)/dT;

% Iterations
for i=1:times
    % Record time
    t(:,i+1) = t(:,i) + dT;
    
    % error
    x_d(:,i+1) = sin(t(:,i+1));
    e = x_1(:,i) - x_d(:,i+1);
    
    dot_x_d(:,i+1) = cos(t(:,i+1));
    dot_e = x_2(:,i) - dot_x_d(:,i+1);
    
    s = c*e + dot_e;
    
    % RBF
    % calculate Gaussian basis function
    %norm 范数函数
    for j = 1:HidLayNum
        h(j,:) = exp( -norm([x_1(:,i); x_2(:,i)]-c_i(:,j))^2/(2*b_i^2) );
    end
    
    hat_fx(:,i+1) = hat_W' * h;
    
    % adaptive law
    hat_W = hat_W + dT * (gamma * s * h);
    
    % control input
    ddot_x_d(:,i+1) = -sin(t(:,i+1));
    u(:,i+1) = -c*dot_e + ddot_x_d(:,i+1) - eta * sign(s) - hat_fx(:,i+1);
    
    % update states
    x_2(:,i+1) = x_2(:,i) + dT * ( fx(:,i)+u(:,i+1) );
    x_1(:,i+1) = x_1(:,i) + dT * x_2(:,i+1);
    
    % update unknown
    fx(:,i+1) = 10 * x_1(:,i+1) * x_2(:,i+1);
end

%% Plot results
figure(1)
subplot(2,1,1)
plot(t,x_1, t,x_d, 'linewidth',1.5);
legend('$x_{1}$', '$x_{d}$', 'interpreter','latex');
grid on

subplot(2,1,2)
plot(t,x_2, t,dot_x_d, 'linewidth',1.5);
legend('$x_{2}$', '$\dot{x}_{d}$', 'interpreter', 'latex');
grid on

figure(2)
subplot(2,1,1)
plot(t,u, 'linewidth',1.5);
legend('$u$', 'interpreter','latex');
grid on

subplot(2,1,2)
plot(t,fx, t,hat_fx,'k:', 'linewidth',1.5);
legend('$f(x)$', '$\hat{f}(x)$', 'interpreter','latex');
grid on


function out = sat(s)
    out = sign(s) * min(10,abs(s));
end
