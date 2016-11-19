function [eig1,eig2,vect1,vect2]=SolveEigenvalues(a_diff,alpha_m,beta1_m,beta2_m,StrucParam)

a_diff_col=a_diff(StrucParam.N_Tr:1:2*StrucParam.N_Tr-1);
a_diff_row=a_diff(StrucParam.N_Tr:-1:1);
a_matr=toeplitz(a_diff_col,a_diff_row); % toeplitz matrix a_diff - equation (11) in the paper

Matr1=[-diag(1./(beta1_m.^2))*(diag(alpha_m)*a_matr+a_matr*diag(alpha_m)) diag(1./(beta1_m.^2))*(eye(StrucParam.N_Tr)+a_matr*a_matr);eye(StrucParam.N_Tr) zeros(StrucParam.N_Tr)];
Matr2=[-diag(1./(beta2_m.^2))*(diag(alpha_m)*a_matr+a_matr*diag(alpha_m)) diag(1./(beta2_m.^2))*(eye(StrucParam.N_Tr)+a_matr*a_matr);eye(StrucParam.N_Tr) zeros(StrucParam.N_Tr)];

[vect1,eig1]=eig(Matr1);
eig1=1./diag(eig1);
[vect2,eig2]=eig(Matr2);
eig2=1./diag(eig2);

end