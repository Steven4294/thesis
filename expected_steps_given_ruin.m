function [E_r] = expected_steps_given_ruin(x0, x1, delta_t)

T = norm(x1-x0);
n = delta_t;
m = T - n;

E_r = (T/m)*(((2/3)*n*T)-n^2+((1/(3*T))*n^3));