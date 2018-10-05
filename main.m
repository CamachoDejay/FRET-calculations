% this script performs the simulations presented in the article: 2D
% polarization imaging as a low-cost fluorescence method to detect
% alpha-synuclein aggregation ex vivo by Camacho et.al.
%
% In general this code simulates the fluorescence emission of an
% anisotropic ensemble of GFP molecules after excitation by polarized
% light. In the simulations we consider the presence of homo-FRET between
% GFP molecules within the clasical Forster Teory.

% clearing of workspace
clear
% closing all figures
close all
% clearing of the command window
clc

%%%%%% User inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distances: interchromophoric distances. chromophores are placed in a
% cubic lattice, distances is the distance between the lattice points in
% nanometers.
distances = [4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 10 11 12 14 16]; % [nm]

% dimer_probs: dimer probabilities. probability that a place of the lattice
% contains a dimer instead of a monomer.
dimerProbs = [0 1/3 2/3 1];

% nReplicates: number of times we repeat the simulation - Due to the
% complexity of the problem we can only provide a numerical solution to the
% homo-FRET process in a lattice of finite size. Therefore, it is expected
% that the simulation yields slightly different results when repeated.
nReplicates = 10;

% nIterations: number of iterations. In each iteration step the
% polarization portrait obtained has fully polarized excitation (single
% absorbing dipole, i.e. central dipole) and an emission polarization that
% depends on the interaction between the central dipole and its buffer
% (neighbouring dipoles in the lattice). Therefore, to simulate the
% response coming from a large ensemble of randomly oriented dipoles (e.g.
% GFP in solution) many (hundreds-thousands) dipole model iterations have
% to be done and sum together.
nIterations = 500;

% buffer properties: properties of the neighbouring chromophores, i.e.
% geometrical properties of the lattice.
%   max_size: max dimension of the considered buffer in nm.
buffer.max_size = 500; % [nm]
%   min_size: if due to numerical reasons we have to use a buffer smaller
%   than the min_size a warning is trigered.
buffer.min_size = 10; % [nm]
%   max_n_dipoles: max number of dipoles that can be inside of the buffer.
%   The larger this is the more time comsuming the calculation, as the
%   energy transfer matrix has a size [n_dipoles n_dipoles]
buffer.max_n_dipoles = 200; 
%%%%%% END OF USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variables related to the chromophored used, calculation of spec overlap
% properties related to the FRET process

% load 'GFPemerald' matrix containing the absorption and emission spectra
% of Emerald GFP. This data was obtained from Thermo Fisher and confirmed
% via experiments.
%   load
load('GFPemeraldSpecData.mat')
% abs spec of the acceptor 
absA = GFPemerald(:,2); 
% emission spec of the donoer
emD  = GFPemerald(:,3); 
% common wavelength axis
wavelength = GFPemerald(:,1); % [nm] 
% extinction coef of the acceptor at max abs
extCoefA = 57500;     % [mol^-1 cm^-1] 
% These numbers are very important for the Forster equation as we have to
% calculate the spectral overlap between dipoles.
% refractive index
ri = 1.4;       % aqueous solution
% quantum yield of the donor
QY    = 0.68;      % GFP emerald
% lifetime of the donor
lifeD  = 2.7;       % [ns] GFP 
% creating the homoFRET object containing all the relevant FRET properties
FRETprops = homoFRET(absA,emD,wavelength,extCoefA,ri,QY,lifeD);

% Forster Radius calculation
ForRad = FRETprops.getFRadius;

%% ATTO550 is not used in the symulations.Moreover, I calcualte the Förster
% radius for comparison
load('ATTO550SpecData.mat')
% abs spec of the acceptor 
absA = ATTO550(:,2); 
% emission spec of the donoer
emD  = ATTO550(:,3); 
% common wavelength axis
wavelength = ATTO550(:,1); % [nm] 
% extinction coef of the acceptor at max abs
extCoefA = 1.2e5;     % [mol^-1 cm^-1] 
% These numbers are very important for the Forster equation as we have to
% calculate the spectral overlap between dipoles.
% refractive index
ri = 1.4;       % aqueous solution
% quantum yield of the donor
QY    = 0.8;      % GFP emerald
% lifetime of the donor
lifeD  = 3.6;       % [ns] GFP 
% creating the homoFRET object containing all the relevant FRET properties
FRETpropsATTO = homoFRET(absA,emD,wavelength,extCoefA,ri,QY,lifeD);

% Forster Radius calculation
ForRadATTO = FRETpropsATTO.getFRadius;

% print message
fprintf('Förster Radius of GFP emerald in aqueous solution: %0.2f\n',ForRad)
fprintf('Förster Radius of ATTO550 in aqueous solution: %0.2f\n',ForRadATTO)
%% simulating polarization portraits
% Properties of the polarization portrait to simulate
%   excitation angles in degrees
ex_angles_d = 0:1:180;
%   0 and 180 is the same angle so I remove it
ex_angles_d(ex_angles_d == 180) = [];
%   emission angles in degrees
em_angles_d = 0:1:180;
%   0 and 180 is the same angle so I remove it
em_angles_d(em_angles_d == 180) = [];
%   populate the portait porperties
portrait_prop.ex_angles = ex_angles_d;
portrait_prop.em_angles = em_angles_d;

% init variables used in the for loops for prealocation of memory.
%   number of distances considered
nDistances = length(distances);
%   number of dimer probabilities considered
nDimProb  = length(dimerProbs);
%   init of simulated polarization portraits
portraits = zeros(length(ex_angles_d),length(em_angles_d),nDistances,nDimProb,nReplicates);
%   variable to keep track of the actual buffer used
buffUsed = zeros(nDistances,nDimProb,nReplicates);

% timing of the code to give an estimate of the run time left to the user
tic
% t1: time point 1
t1 = toc;
% iteration over the number of replicates
for n = 1:nReplicates
    % t3: timpoint 3;
    t3 = toc;
    % iteration over the considered distances (dimentions of the lattice)
    for d = 1:nDistances
        % interCdist: inter chromophoric distance - distance between the
        % points of the squared lattice.
        interCdist = distances(d);
        % initializing dipole model object. Dinit: As of this moment it
        % only knows about the general properties of the lattice buffer to
        % use.
        Dinit = Dipoles.dipole_model('3D',buffer);
        % t4: time point 4
        t4 = toc;
        % iteration over the dimer probabilities considered
        for dp = 1:nDimProb
            % dimer probability to use
            dimerP = dimerProbs(dp);
            % pass information about dimer probability to the dipole model.
            % As Dinit can be reused we store the output in a new object: D
            D = Dinit.get_positions(interCdist,dimerP);
            % get the orientation of all dipoles in the lattice
            D = D.get_orientations;
            % once we know the position and angles of all dipoles we can
            % calculate the energy transfer matrix of the whole system
            D = D.get_et_matrix(FRETprops);
            % we are now able to calculate the polarization portrait if we
            % selectively excite the central dipole and observe the
            % emission of the whole system after the energy trasnfer
            % process has taken place.
            P = Portrait.pol_portrait(D,portrait_prop);
            % as mentioned above to simulate the properties of a large
            % ensemble we must simualate many dipole models (P) and sum
            % the results.
            %   We first init the total polarization portrait: totI to the
            %   results of our first iteration.
            totI = P.I_ex_em;
            %   We also store the buffer used as is common to all D
            buffUsed(d,dp,n) = D.buffer_used;
            % then we iterate nIterations times and sum the retuls of each
            % run to totI. As this is the most time consuming step we do a
            % parfor here.
            parfor i = 1:nIterations 
                % we re run the dipole model simulation many times and keep
                % track of the total result
                D = Dinit.get_positions(interCdist,dimerP);
                D = D.get_orientations;
                D = D.get_et_matrix(FRETprops);
                Ptmp = Portrait.pol_portrait(D,portrait_prop);
                totI = totI + Ptmp.I_ex_em;

            end
            % store output
            portraits(:,:,d,dp,n) = totI;
            % inform the user that all goes well
            disp(['Done for Distance: ' num2str(interCdist)...
                  '; Dimer_prob: ' num2str(dimerP)])
        end
        % t5: timepoint 5
        t5 = toc;
        % telling the user how long it is all taking
        disp(['time taken for distance ' num2str(interCdist) ': ' num2str(round(t5-t4,2)) ' sec'])

    end
    % t6: timepoint 6
    t6 = toc;
    % telling the user how long it is all taking
    disp(['time taken for replicate ' num2str(n) ': ' num2str(round((t6-t3)/60,2)) ' min'])
end
% t7: timepoint 7
t7 = toc;
% telling the user how long it is all taking
disp(['time taken for sim: ' num2str(round((t7)/60,2)) ' min'])
% we have now finished with the polariation portraits calculation
disp('Done with portrait calculation')
% save the output
run_prop_name = ['Ndip' num2str(buffer.max_n_dipoles)...
                 '-Nite' num2str(nIterations) ...
                 '-Nrep' num2str(nReplicates)];             
results_dir = [cd filesep 'Results' filesep run_prop_name filesep];
mkdir(results_dir)
save([results_dir 'portraits.mat'],'portraits')
save([results_dir 'distances.mat'],'distances')

% clearing tmp vars
clear t1 t3 t4 t5 t6 t7 d dimer_p dp inter_c_dist n   
clear n_iterations 
%% calculating polim params
% now that the portraits are done we must calculate the POLIM output for
% each portrait and store the resuts for later plotting and report.

% init of variables used to store the polim output: outParams
%   init structure outParams
outParams = [];
%   nCond: number of simulated conditions
nCond = nDistances * nDimProb;
%   init the fields of outParams
outParams(nCond).Distance      = [];
outParams(nCond).Dimer_prob    = [];
outParams(nCond).Mex           = [];
outParams(nCond).Mex_dev       = [];
outParams(nCond).Mem           = [];
outParams(nCond).Mem_dev       = [];
outParams(nCond).r_polim       = [];
outParams(nCond).r_polim_dev   = [];
outParams(nCond).r_normal      = [];
outParams(nCond).r_normal_dev  = [];
outParams(nCond).Epsilon       = [];
outParams(nCond).Epsilon_dev   = [];

% angles used for the anisotropy calculation - radians
r_ex_deg = round(P.ex_angles_rad * 180 / pi);
r_em_deg = round(P.em_angles_rad * 180 / pi);

% by choosing the portrait angles properly this calculation becomes very
% easy as the excitation angle 00 exists and the emission angles 00 and 90
% deg also exist
ex_00_deg = r_ex_deg == 0;
test_1 = sum(ex_00_deg) == 1;
em_00_deg = r_em_deg == 0;
test_2 = sum(em_00_deg) == 1;
em_90_deg = r_em_deg == 90;
test_3 = sum(em_90_deg) == 1;
assert(test_1,'problmnes with classical r calculation: can not find ex = 0 degrees');
assert(test_2,'problmnes with classical r calculation: can not find em = 0 degrees');
assert(test_3,'problmnes with classical r calculation: can not find em = 90 degrees');
clear test_1 test_2 test_3

% iteration over the different distances
% a simple counter used for indexing
counter = 0;
for i = 1:nDistances
    % iteration over the different dimer probabilities
    for k = 1:nDimProb
        % increasing the counter used for indexing
        counter = counter + 1;
        % fishing out the results of interest
        results_ik = outParams(counter);
        % init of structure array
        results_j(nReplicates) = results_ik;

        % iteration over the replicates
        for j = 1:nReplicates
            % indexing to the portrait of interest
            portrait_ikj = portraits(:,:,i,k,j);
            % calculatinf the anisotropy
            I_par = portrait_ikj(ex_00_deg,em_00_deg);
            I_per = portrait_ikj(ex_00_deg,em_90_deg);
            results_j(j).r_normal = (I_par - I_per) / (I_par + 2*I_per);
            % calculating the POLIM params
            P = P.modify_intensity(portrait_ikj);
            polimOUT  = POLIM.fitSinglePortrait( P );
            % storing POLIM results
            results_j(j).Mex = polimOUT.Mex;
            results_j(j).Mem = polimOUT.Mem;
            results_j(j).r_polim = polimOUT.Anisotropy;
            results_j(j).Epsilon = polimOUT.SFA.epsilon;
        end
        
        % storing results
        %   Distance
        results_ik.Distance   = distances(i);
        results_ik.Dimer_prob = dimerProbs(k);
        %   excitation modulation depth 
        results_ik.Mex     = nanmean([results_j.Mex]);
        results_ik.Mex_dev = nanstd([results_j.Mex]);
        %   emission modulation depth
        results_ik.Mem     = nanmean([results_j.Mem]);
        results_ik.Mem_dev = nanstd([results_j.Mem]);
        %   fluorescence anisotropy - classical calculation
        results_ik.r_normal = nanmean([results_j.r_normal]);
        results_ik.r_normal_dev = nanstd([results_j.r_normal]);
        %   fluorescence anisotropy - POLIM optimized for orientational
        %   artifacts
        results_ik.r_polim = nanmean([results_j.r_polim]);
        results_ik.r_polim_dev = nanstd([results_j.r_polim]);
        %   energy funneling efficiency
        results_ik.Epsilon = nanmean([results_j.Epsilon]);
        results_ik.Epsilon_dev = nanstd([results_j.Epsilon]);
        %   placing it at the right index
        outParams(counter) = results_ik;
        % informing the user that all goes well
        disp(['Done for distance: ' num2str(results_ik.Distance)...
              '; dimer_prob: ' num2str(results_ik.Dimer_prob)])
        % clearing out tmp variable 
        clear results_j  
    end
    
    
end

% saving results
save([results_dir 'outParams.mat'],'outParams')
% clearing tmp vars
clear i I_par I_per j k n_cond n_replicates portrait_ikj r_em_deg r_ex_deg
clear results_ik em_00_deg em_90_deg ex_00_deg ex_angles_d em_angles_d
clear n_distances

%% generating figure
figure(1)
subplot(2,1,1)
legend_str = cell(nDimProb,1);
for i = 1:nDimProb
    dp = [outParams.Dimer_prob];
    di = dimerProbs(i);
    idx = dp==di;
    
    errorbar([outParams(idx).Distance],[outParams(idx).r_normal],[outParams(idx).r_normal_dev])
    hold on
    legend_str{i} = ['Dimer sites: ' num2str(round(di*100)) '%'];
end
hold off
xlabel('distance in nm')
ylabel('fluorescence anisotrypy')
ylim([0 0.45])
xlim([3.5 16.5])
title({run_prop_name, 'FA'})
set(gca,'YTick',0:0.1:4)
grid on
xlabel('distance in nm')
ylabel('fluorescence anisotrypy')
legend(legend_str, 'Location', 'eastoutside' )

subplot(2,1,2)

for i = 1:nDimProb
    dp = [outParams.Dimer_prob];
    di = dimerProbs(i);
    idx = dp==di;
    
    
    errorbar([outParams(idx).Distance],[outParams(idx).Epsilon],[outParams(idx).Epsilon_dev])
    hold on

end
plot([4 16],[.5 .5],'--r','linewidth',1)
hold off
ylim([0 1])
xlim([3.5 16.5])
set(gca,'YTick',0:0.25:1)
grid on
xlabel('distance in nm')
ylabel('Epsilon')
legend(legend_str, 'Location', 'eastoutside' )

title({run_prop_name, 'Eplison'})