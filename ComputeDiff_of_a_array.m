function [a_diff]=ComputeDiff_of_a_array(n_set,StrucParam)

a_diff=zeros(1,length(n_set));

it=1;
for n=n_set(1):1:n_set(end),
    F=@(x)(1/StrucParam.dx)*exp(-1i*n*2*pi*x/StrucParam.dx).*eval(StrucParam.diff_a_x);
    a_diff(it)=quadgk(F,0,StrucParam.dx);
    it=it+1;
end    

if StrucParam.cut==1; %cut small elements in the array
    a_diff=CutSmallArrayElements(a_diff,StrucParam.accuracy);
end

end