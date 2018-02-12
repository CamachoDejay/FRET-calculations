function [ fit_results ] = fitSinglePortrait( P )
%FITSINGLEPORTRAIT Calculation of the POLIM parameter based in the input
% portrait
%   TODO

assert(isa(P,'Portrait.pol_portrait'), 'input must be a pol_portrait object')

port   = P.I_ex_em;
ex_ang = P.ex_angles_rad;
em_ang = P.em_angles_rad;

% modulation in excitation %%
I_ex = sum(port,2);
I_ex = I_ex ./ sum(I_ex);
[Y, ~] = CosineFitFunc(ex_ang,I_ex,1);
m_ex           = Y(1);
fi0_ex         = Y(2)*pi/180;


% modulation in emission %%
I_em = sum(port,1);
I_em = I_em ./ sum(I_em);
[Y, ~] = CosineFitFunc(em_ang,I_em',1);
m_em           = Y(1);
fi0_em         = Y(2)*pi/180;

% anisotropy - classical:
tmp = abs(ex_ang - 0);
[val1,ex_par_ind] = min(tmp);

tmp   = abs(em_ang - 0);
[val2,em_par_ind] = min(tmp);

tmp   = abs(em_ang - pi/2);
[val3,em_per_ind] = min(tmp);


if all([val1, val2, val3] < pi/180)
    % then we can define anisotropy in the classical way
    Ipar = port(ex_par_ind,em_par_ind);
    Iper = port(ex_par_ind,em_per_ind);
    r = (Ipar - Iper) / (Ipar + 2*Iper);
else
    warning('current version of the code can NOT fit anisotropy, this will come later on')
    r = [];
    
end

% ansotropy POLIM
Fex_par = mod(fi0_ex,pi);
Fex_per = mod(fi0_ex+(pi/2),pi);

tmp = abs(ex_ang - Fex_par);
[val1,ex_par_ind] = min(tmp);

tmp   = abs(em_ang - Fex_par);
[val2,em_par_ind] = min(tmp);

tmp   = abs(em_ang - Fex_per);
[val3,em_per_ind] = min(tmp);

if all([val1, val2, val3] < pi/180)
    % then we can define anisotropy in the POLIM way
    Ipar = port(ex_par_ind,em_par_ind);
    Iper = port(ex_par_ind,em_per_ind);
    rPOLIM = (Ipar - Iper) / (Ipar + 2*Iper);
else
    warning('current version of the code can NOT fit anisotropy, this will come later on')
    rPOLIM = [];
    
end



% Luminescence shift
LS = fi0_ex - fi0_em; 
if LS>pi/2
    LS=LS-pi;
elseif LS<-pi/2
    LS=LS+pi;
end

% Average intensity
% AvInt = mean(port(:));
           
[SFAoutput, model_Mcurves] = polimETcalculation(P,m_ex,fi0_ex);
                                       
fit_results.Mex        = m_ex;
fit_results.Fex        = fi0_ex;
fit_results.Mem        = m_em;
fit_results.Fem        = fi0_em;
fit_results.LS         = LS;
fit_results.Anisotropy = r;
fit_results.rPOLIM     = rPOLIM;
fit_results.SFA        = SFAoutput;
end

