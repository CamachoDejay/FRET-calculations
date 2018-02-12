function [Y, Fit_Curve]=CosineFitFunc(Angles,DataToFit,phase_res)

    m   = length(Angles);
    tmp = min(size(Angles));
    
    assert(m > 1 && tmp ==1, 'Input Angles must be a vector')
    shapes=ones(m,2);

    if nargin == 2
        p_res = 1;

    elseif nargin == 3
        p_res = phase_res;
    end

    phases = 0:p_res:180;
    phases(phases==180) = [];
    y = NaN(length(phases),2);

    %%Fitting
     for i=1:length(phases)
        phase       = (phases(i)-90);
        shapes(:,2) = cos(2*(Angles-(phase*pi/180)));

        x       = lscov(shapes,DataToFit);
        eps     = abs(x);
        y_fit   = shapes*eps;
        y_error = DataToFit-y_fit;
        y(i,1)  = phase;
        y(i,2)  = sum(sum(y_error.*y_error));
     end
    %%Fitting output
    [Fit_Error, y_min_index] = min(y(:,2));
    Fit_Phase = y(y_min_index,1);

    shapes(:,2) = cos(2*(Angles-(Fit_Phase*pi/180)));
    x           = lscov(shapes,DataToFit);

    if ~all(x >= 0)
        warning('problems in the cosine fitting')
    end

    eps         = abs(x);
    Fit_Curve   = shapes*eps;

    % Imax = eps(1) + eps(2);
    % Imin = eps(1) - eps(2);
    % M = Imax-Imin / Imax+Imin;
    Fit_Modulation = eps(2) / eps(1);
    Y = [Fit_Modulation Fit_Phase Fit_Error];
end