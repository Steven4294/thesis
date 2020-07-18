function [E_nr] = expected_steps_given_win(x0, x1, delta_t)

T = norm(x1-x0);
n = delta_t;
m = T - n;

E_nr = (T/n)*(((2/3)*m*T)-m^2+((1/(3*T))*m^3));


