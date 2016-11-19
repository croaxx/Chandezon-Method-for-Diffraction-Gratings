function [m_set, n_set, alpha_m, beta1_m, beta2_m, n1_set, n1_set_ind, n2_set, n2_set_ind]=CalculateAlphaBetaAndMNSets(StrucParam)

K = 2*pi/StrucParam.dx; %inverse lattice vector
m1 = -((StrucParam.N_Tr-1)/2);
m2 = -m1;
m_set = m1:1:m2; % array of m according to Li
n_set = 1-StrucParam.N_Tr:1:StrucParam.N_Tr-1; % array of n according to Li
alpha_m = StrucParam.n1*(2*pi/StrucParam.lambda)*sin(StrucParam.theta)+K*m_set;
beta1_m = diag(sqrt(eye(StrucParam.N_Tr)*(StrucParam.n1*2*pi/StrucParam.lambda)^2-diag(alpha_m.^2)));
beta2_m = diag(sqrt(eye(StrucParam.N_Tr)*(StrucParam.n2*2*pi/StrucParam.lambda)^2-diag(alpha_m.^2)));

idx=find(~imag(beta1_m));
n1_set = m_set(idx);
n1_set_ind = idx.';

idx=find(~imag(beta2_m));
n2_set = m_set(idx);
n2_set_ind = idx.';

end