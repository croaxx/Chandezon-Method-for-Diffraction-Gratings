function [eig_p,vect_p,eig_m,vect_m]=SortEigenvaluesAndVectors(eig,vect, accuracyImag)
%% sortation of eigenvalues and corresp eigenvectros
idxRealAndPositive = abs(imag(eig)) < accuracyImag & real(eig) > 0;
idxRealAndNegative = abs(imag(eig)) < accuracyImag & real(eig) <= 0;
idxImagAndPositive = abs(imag(eig)) >= accuracyImag & imag(eig) > 0;
idxImagAndNegative = abs(imag(eig)) >= accuracyImag & imag(eig) <= 0;

eig_real_p = eig(idxRealAndPositive).';
vec_real_p = vect(:,idxRealAndPositive);
eig_real_m = eig(idxRealAndNegative).';
vec_real_m = vect(:,idxRealAndNegative);
eig_comp_p = eig(idxImagAndPositive).';
vec_comp_p = vect(:,idxImagAndPositive);
eig_comp_m = eig(idxImagAndNegative).';
vec_comp_m = vect(:,idxImagAndNegative);

%% sortation in descending order
if ~isempty(eig_real_p)
    [~,ind]=sort(real(eig_real_p),'descend');
    eig_real_p=eig_real_p(ind);
    vec_real_p=vec_real_p(:,ind);
end

if ~isempty(eig_comp_p)
    [~,ind]=sort(imag(eig_comp_p),'ascend');
    eig_comp_p=eig_comp_p(ind);
    vec_comp_p=vec_comp_p(:,ind);
end

if ~isempty(eig_real_m)
    [~,ind]=sort(real(eig_real_m),'ascend');
    eig_real_m=eig_real_m(ind);
    vec_real_m=vec_real_m(:,ind);
end

if ~isempty(eig_comp_m)
    [~,ind]=sort(imag(eig_comp_m),'descend');
    eig_comp_m=eig_comp_m(ind);
    vec_comp_m=vec_comp_m(:,ind);
end

eig_p = [eig_real_p, eig_comp_p];
vect_p = [vec_real_p, vec_comp_p];
eig_m = [eig_real_m, eig_comp_m];
vect_m = [vec_real_m, vec_comp_m];

end