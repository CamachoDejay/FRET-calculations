function [output, model_Mcurves] = polimETcalculation(P,M_ex,Phase_ex)
                                    
    % TODO: this function can be improved by fitting Epsilon via a least
    % square method and combine that with a f min search for the other
    % parameters. This is the way we do it in the python code. When I have
    % the time I will fix it. 

    assert(isa(P,'Portrait.pol_portrait'),...
           'input must be a pol_portrait object')
    
    portrait = P.I_ex_em;
    ex_ang_r = P.ex_angles_rad;
    em_ang_r = P.em_angles_rad;
    
        
    %% Calcualtion using single funnel approximation %%
%     ind = 1:10:180;
%     portrait = portrait(ind,ind);
%     ex_ang_r = ex_ang_r(ind);
%     em_ang_r = em_ang_r(ind);
    

    % get data to be used in fit - 2D representation of the portrait.
    [ExGrid,EmGrid]=meshgrid(ex_ang_r,em_ang_r);
    ex_ang = ExGrid(:);
    em_ang = EmGrid(:);
    % transponse is needed to have the right correspondance between angles
    % and intensity    
    exp_int = reshape(portrait',size(portrait,1)*size(portrait,2),1);
    % normalization of the intensity to its max value
    Ftotal = exp_int ./ max(exp_int);
    % this is probably not needed unless I fit only part of the portrait
    % which is probably a good idea. 
%     FitLeght = length(Ftotal);
    % I want angles to be a one row vector size = [1xn].
    ex_ang = ex_ang';
    em_ang = em_ang';

    %%%%%%%%% inputs for the fit routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   initial guess for fit
    a0 = [M_ex 0 1];   % Mf, thetaf, X
    
    % fit is based in a bounded fminsearch thus I have to define the
    % boundary conditions: boundary 1 is for the modulation of the funnel,
    % boundary 2 is for the relative phase of the funnel and boundary 3 is
    % for the geometrical ratio
    
    %     M_f   P_f     X
    LB = [0.01 -pi/2    0]; % Lower Boundary 
    UB = [1     pi/2    2*(1+M_ex)/(1-M_ex)]; % Upper Boundary
    
    % Note that in geometrical model we allow for the side dipoles to have
    % a different size relative to the main dipole. In some cases, this
    % allows to the 3 dipole model to mimic a 2 dipole system.
    
    % parameters not used from the fmincon function
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];
    
    % optimization settings 
    opt = optimset('MaxFunEvals',10000,'MaxIter',4000,...
                   'TolFun',1.0e-12,'Display','off'); % old value of tol 1e-12
               
    % input for function to minimize

    ExpInput.Mex= M_ex;
    ExpInput.Pex = Phase_ex;
    ExpInput.ExAng = ex_ang;
    ExpInput.EmAng = em_ang;
    ExpInput.Ftot = Ftotal;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Actual fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [a, RMSD] = fmincon(@(a) fun2min2D(a,ExpInput),a0,A,b,Aeq,beq,LB,UB,nonlcon,opt);

    Mf     = a(1);
    thetaf = a(2);
    X      = a(3);

    TwoDmodel = POLIM.twoDmodel(M_ex, Phase_ex, Mf, thetaf, X, [], ex_ang, em_ang);
    [ ~, epsilon, ~, ~ ] = SFA_fit_lsqnonneg( TwoDmodel, ExpInput.Ftot(:));
       
    output.epsilon = epsilon;
    output.Mf = Mf;
    output.Pf = thetaf;
    output.X = X;
    output.RMSD = RMSD;
    
    % generate fitted function using full range of angles 

    TwoDmodel = POLIM.twoDmodel(M_ex, Phase_ex, Mf, thetaf, X, epsilon, ex_ang, em_ang);
    
    fitPlot = reshape(TwoDmodel.intensity,length(unique(TwoDmodel.ex_ang_rad)),length(unique(TwoDmodel.em_ang_rad)));
    fitPlot = fitPlot';
    
    Pfit = P.modify_intensity(fitPlot);
    output.portrait = Pfit;
    
    % modulation in excitation for fit
    I_ex = sum(fitPlot,2);
    model_I_ex = I_ex ./ sum(I_ex);
    [Y, Fit_model_I_ex] = CosineFitFunc(ex_ang_r,model_I_ex,1);
    model_m_ex          = Y(1);
    model_fi0_ex        = Y(2)*pi/180;
    
    output.M_ex = model_m_ex;
    output.P_ex = model_fi0_ex;

    % modulation in emission for fit
    I_em = sum(fitPlot,1);
    model_I_em = I_em ./ sum(I_em);
    [Y, Fit_model_I_em] = CosineFitFunc(em_ang_r,model_I_em',1);
    model_m_em          = Y(1);
    model_fi0_em        = Y(2)*pi/180;
    output.M_em = model_m_em;
    output.P_em = model_fi0_em;

 
    % saving model curves and fit
    model_Mcurves = [model_I_ex; Fit_model_I_ex; model_I_em'; Fit_model_I_em]; 
        