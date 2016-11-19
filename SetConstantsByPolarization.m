function opticalConstants=SetConstantsByPolarization(n1,n2,Pol)

eps0 = 8.8541878171e-12; % permittivity of vacuum
mu0 = 12.5663706141e-7;  % permeability of vacuum

% TM polarization
mu1 = 1;
mu2 = 1;
eps1 = n1^2/mu1;
eps2 = n2^2/mu2;

% for TE polarization the following transformations are done: E <--> H, eps0 <--> -mu0, and eps <--> -mu
if strcmp(Pol,'TE') == 1,
    
    [mu0, eps0] = deal(-eps0,-mu0);
    [mu1, eps1] = deal(eps1,-mu1);
    [mu2, eps2] = deal(eps2,-mu2);
    
end

Z0=sqrt(mu0/eps0);

opticalConstants = struct('n1', n1, 'n2', n2, ...
                          'eps0', eps0, 'eps1', eps1, 'eps2', eps2, ...
                          'mu0', mu0, 'mu1', mu1, 'mu2', mu2, ...
                          'Z0', Z0);

end