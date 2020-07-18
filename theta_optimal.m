function [theta] = theta_optimal(xp, x1, x2, x3, theta1, theta2, theta3)
% takes in 3 curves and returns theta s.t. theta minimizes the 3 lengths
% param = (xp, x1, x2, x3, theta1, theta2, theta3)
theta = -1;
theta_p = 0;
theta_tilda_min = 1000000;

theta_array = linspace(0, 2*pi, 5000);
theta_min = 30;

for i = 1:length(theta_array)
    theta_p = theta_array(i);
    theta_tilda_tot = 1000;
    
    [a0, a1, a2] = quadratic_curve(xp, x1, theta_p, theta1);
    theta_tilda1 = theta_tilda_calc([a0, a1, a2]);

    [a0, a1, a2] = quadratic_curve(xp, x2, theta_p+(2*pi/3), theta2);
    theta_tilda2 = theta_tilda_calc([a0, a1, a2]);

    [a0, a1, a2] = quadratic_curve(xp, x3, theta_p +(4*pi/3), theta3);
    theta_tilda3 = theta_tilda_calc([a0, a1, a2]);

    theta_tilda_tot = theta_tilda1 + theta_tilda2 + theta_tilda3;

    if theta_tilda_tot < theta_tilda_min
        theta_tilda_min = theta_tilda_tot;
        theta_min = theta_p;
    end
end

theta = theta_min;

