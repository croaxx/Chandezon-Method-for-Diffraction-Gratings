function L=L_eval(gamma, m, StrucParam)
FF=@(x) (1/StrucParam.dx)*exp(1i*gamma*eval(StrucParam.a_x) - 1i*m*2*pi*x/StrucParam.dx);
L=quadgk(FF, 0, StrucParam.dx);
end