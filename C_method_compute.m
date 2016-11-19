function [R_tot,T_tot]=C_method_compute(StrucParam)

[m_set, n_set, alpha_m, beta1_m, beta2_m, n1_set, n1_set_ind, n2_set, n2_set_ind] = CalculateAlphaBetaAndMNSets(StrucParam); % basic parameters
a_diff = ComputeDiff_of_a_array(n_set,StrucParam);

[eig1, eig2, vect1, vect2] = SolveEigenvalues(a_diff, alpha_m, beta1_m, beta2_m, StrucParam); %solve the eigenvalue problems in two media
[L_mn_beta1, L_mk_beta2, L_mo_beta1] = Compute_L_arrays(m_set, StrucParam, beta1_m, beta2_m, n1_set, n1_set_ind, n2_set, n2_set_ind); %computation of L_mn, L_mk, and L_mo matrices

[eig1_p, vect1_p, ~, ~] = SortEigenvaluesAndVectors(eig1, vect1, StrucParam.kVecImagMin);
[~, ~, eig2_m, vect2_m] = SortEigenvaluesAndVectors(eig2, vect2, StrucParam.kVecImagMin);

%% F_mn_R_p, F_mk_R_m, F_mo_R_in, F_mq_p, F_mr_m
F_mn_R_p = L_mn_beta1;
F_mk_R_m = L_mk_beta2;
F_mo_R_in = L_mo_beta1;
F_mq_p = vect1_p(1:StrucParam.N_Tr,length(n1_set)+1:end);
F_mr_m = vect2_m(1:StrucParam.N_Tr,length(n2_set)+1:end);

%% cut small elements, if the flag is set
if StrucParam.cut == 1,
    F_mn_R_p = CutSmallArrayElements(F_mn_R_p, StrucParam.accuracy);
    F_mk_R_m = CutSmallArrayElements(F_mk_R_m, StrucParam.accuracy);
    F_mo_R_in = CutSmallArrayElements(F_mo_R_in, StrucParam.accuracy);
    F_mq_p = CutSmallArrayElements(F_mq_p, StrucParam.accuracy);
    F_mr_m = CutSmallArrayElements(F_mr_m, StrucParam.accuracy);
end

[G_mn_R_p, G_mk_R_m, G_mo_R_in, G_mq_p, G_mr_m] = Compute_G_matrices(m_set, n1_set, n1_set_ind, n2_set, n2_set_ind, ...
                                                                     alpha_m, a_diff, beta1_m, beta2_m, L_mn_beta1, L_mk_beta2, ...
                                                                     L_mo_beta1, eig1_p, F_mq_p, eig2_m, F_mr_m, StrucParam);

GF_matrix = [F_mn_R_p F_mq_p -F_mk_R_m -F_mr_m;
             G_mn_R_p G_mq_p -G_mk_R_m -G_mr_m];

GF_col = -[F_mo_R_in; G_mo_R_in];

Amplitudes = GF_matrix\GF_col;

%% calculation of reflection efficiencies
R=zeros(1,length(n1_set));
for it=1:length(n1_set),
    R(it) = (beta1_m(n1_set_ind(it))/beta1_m((StrucParam.N_Tr+1)/2))*(abs(Amplitudes(it)))^2;
end
R_tot = sum(R);

%% calculation of transmission efficiencies
T_tot = 0;
T = zeros(1, length(n1_set));
if ~isempty(n2_set),
    for it = 1:length(n2_set),
        T(it) = (StrucParam.eps1*beta2_m(n2_set_ind(it))/(StrucParam.eps2*beta1_m((StrucParam.N_Tr+1)/2)))*(abs(Amplitudes(StrucParam.N_Tr+it)))^2;
    end
    T_tot = sum(T);
end

end