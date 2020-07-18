close all;
clear all;

x0 = [0; 0];
x1 = [1; 1];
x2 = [-1; -1];
x3 = [-1; 1];

%x0 = [0; 0];
%x1 = [cos(pi/6); sin(pi/6)];
%x2 = [cos(5*pi/6); sin(5*pi/6)];
%x3 = [0; -1];


points = [x0 x1 x2 x3];
[m n] = size(points); % m = number of points
theta_array = linspace(0, pi, 100);
theta_tilda_min = 10000;


% manual search optimization
for i = 1:length(theta_array)
    theta = theta_array(i);
    theta_tilda_tot = 0;
    
    for i = 1:n-1
        x_end = points(:,i+1);
        [a0 a1 a2] = quadratic_curve(x0, x_end, theta, 0);
        theta_tilda = theta_tilda_calc([a0, a1, a2]);
        theta_tilda_tot = theta_tilda_tot + theta_tilda;
    end

    if theta_tilda_tot < theta_tilda_min
        theta_tilda_min = theta_tilda;
        theta_min = theta;
    end
end


% plotting

theta = theta_min;
for i = 1:n-1
        
    mod(radtodeg(theta), 180)

    x_end = points(:,i+1);
    [a0 a1 a2] = quadratic_curve(x0, x_end, theta, 0);

    s_array = linspace(0,1,500);
    curve = [];
    for i=1:length(s_array)
        s = s_array(i);
        c = [a0(1) + a1(1)*s + a2(1)*s^2 ; a0(2) + a1(2)*s + a2(2)*s^2];
        curve = [curve c];
    end
      
    theta = theta + (2*pi/3);
    figure(1)
    hold on
    plot(curve(1,:), curve(2,:));
    plot(points(1,:), points(2,:), 'o');
    title('manual search - minimum length')
    axis([-2 2 -2 2])

end
    



     
