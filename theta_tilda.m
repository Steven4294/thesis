function [theta] = theta_tilda_calc(param)

theta = 3.2;
% a0 = param(:,1);
% a1 = param(:,2);
% a2 = param(:,3);
% theta = 0;
% alpha = norm(a1)
% beta = 4*dot(a1,a2)
% gamma = 2*norm(a2)
% z = sqrt(alpha^2 + beta + gamma^2); %% what happens if this is complex??? 
% y = beta^2 - (4*alpha^2*gamma^2);
% 
% % taking the log of negatives fails....
% theta = (1/(8*gamma^3))*(2*gamma*(-alpha*beta+z*(beta+2*gamma^2))...
%     + y*log(beta+2*alpha*gamma)-y*log(beta+2*gamma*(gamma+z)));