function RMSD = fun2min2D(a, ExpInput)

Mex = ExpInput.Mex;
Pex = ExpInput.Pex;
MF    = a(1);
thetaf = a(2);
X     = a(3);
ExAng = ExpInput.ExAng;
EmAng = ExpInput.EmAng;

TwoDmodel = POLIM.twoDmodel(Mex, Pex, MF, thetaf, X, [], ExAng, EmAng);

[ ~, ~, RMSD, ~ ] = SFA_fit_lsqnonneg( TwoDmodel, ExpInput.Ftot(:));