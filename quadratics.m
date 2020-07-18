close all;
clear all;
drawArrow = @(x,y) quiver( x(1),x(2),y(1)-x(1),y(2)-x(2),0 )    


X = [0 1 -1 0; 0 1 1 -1];
theta = pi/2;
Ti = [cos(theta); sin(theta)];
Tf = [cos(theta); sin(theta)];

Ni = [-sin(theta); cos(theta)];
Nf = [-sin(theta); cos(theta)];

for i = 1:3
   theta = pi/3*i;
   v = [cos(theta); sin(theta)];
   n = [-sin(theta); cos(theta)];
   Tf = [Tf v];
   Nf = [Nf n];
end

for i = 1:2
   theta = pi/3*i;
   v = [cos(theta); sin(theta)];
   n = [-sin(theta); cos(theta)];
   Ti = [Ti v];
   Ni = [Ni n];
end

for j = 1:3
    
    X0 = X(:,1);
    X1 = X(:,j+1);

    a0 = X0;

    T01 = (2*dot(X1-X0,Nf(:,j+1)))/dot(Ti(:,j),Nf(:,j+1));

    a1 = T01*Ti(:,j);
    a2 = X1-X0-(T01*Ti(:,j));

    s_array = linspace(0,1,500);
    curve = [X0];

    for i=1:length(s_array)
        s = s_array(i);
        c = [a0(1) + a1(1)*s + a2(1)*s^2 ; a0(2) + a1(2)*s + a2(2)*s^2];
        curve = [curve c];
    end

    figure(1)
    hold on
    plot(curve(1,:), curve(2,:));
     
end

hold on
plot(X(1,:),X(2,:),'o')

for i=1:4
    if i < 4
    x1 = X(:,1);
    y1 = x1+Ti(:,i)/5;
    drawArrow(x1, y1 );
    end

    x1 = X(:,i)
    y1 = X(:,i)+Tf(:,i)/5;

    drawArrow(x1, y1 );
end





 