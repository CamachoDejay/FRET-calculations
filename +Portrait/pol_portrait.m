classdef pol_portrait
    %POL_PORTRAIT Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = protected)
        ex_angles_rad
        em_angles_rad
        et_steps
        res_ener
        I_ex_em
        em_after_et
    end

    methods
        function obj = pol_portrait (dipole_model, portrait_prop) % object constructor

            % check that inputs are as expected
            assert(isa(dipole_model,'Dipoles.dipole_model'),...
                   'Input must be a sipole_model object')

            % check that we have some positions
            if isempty( dipole_model.positions )
              warning('dipole model positions missing! returning empty pol_portrait')
              return
            end

            % check that we have some orientations
            if isempty( dipole_model.orientations )
              warning('dipole model orientations missing! returning empty pol_portrait')
              return
            end

            % check that we have a transfer_matrix
            if isempty( dipole_model.transfer_matrix )
              warning('dipole model transfer matrix missing! returning empty pol_portrait')
              return
            end

            if isa(portrait_prop,'struct')
                expected_fields = {'em_angles';'ex_angles'};
                if all(strcmp(sort(fieldnames(portrait_prop)),expected_fields))
                    assert(size(portrait_prop.ex_angles,1)==1 &&...
                           size(portrait_prop.ex_angles,2)>1,...
                           'Excitation angles must be a 1xn array')

                    assert(size(portrait_prop.em_angles,1)==1 &&...
                           size(portrait_prop.em_angles,2)>1,...
                           'Emission angles must be a 1xn array')
%                     disp('Portrait props are ok')
                else
                    error (['portrait_prop fields must be em_angles,'...
                        'ex_angles both in degrees'])
                end
            else
                error(['portrait_prop must be a structure with ',...
                       'fields: ex_angles and em_angles, boths in degrees'])
            end

            test1 = portrait_prop.ex_angles > pi*2;
            assert(any(test1),'Excitation angles must be in degrees')
            test2 = portrait_prop.em_angles > pi*2;
            assert(any(test2),'Emission angles must be in degrees')
            
            % look for central dipole for which the buffer was created
            ind_int = dipole_model.central_dipole;
            assert(sum(ind_int) == 1, 'Unexpected error')
            
            

            % change input angles into radians
            ex_angles_r = portrait_prop.ex_angles .* (pi/180);
            em_angles_r = portrait_prop.em_angles .* (pi/180);

            % load dipole orientation, orientation and ET matrix
%             pos       = dipole_model.positions;
            dip_ori   = dipole_model.orientations;
            ET_matrix = dipole_model.transfer_matrix;
            
            dipole_magnitude = (sum(dip_ori.*dip_ori)).^(.5);
            if any(dipole_magnitude < 0.99999)
                warning('Unexpected behaviour of the dipole magnitude, have a look')
            end


            %%%%%%%%%%%% EXCITATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the electric field vector of the excitation light
            E = get_ex_elect_field(ex_angles_r) ;

            % Once defined the polarization properties of the excitation
            % light we calculate the excitation probabilities. Here, the
            % dot product can be used because it nicely gives:
            % a.b = ||a|| ||b|| cos(theta);
            % where theta is the angle between a and b. Calculations are
            % done for light traveling along the x axis and
            % the electric flield lies on the y,z plane.
            abs_p = (E * dip_ori).^2;
            % in abs_p the first index (row) represents the excitation
            % angle and the second index (col) represents the dipole index

            %%%%%%%% EET   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % What is the effect of EET? so far we have calculated what
            % happens to the excitation energy once it reaches a dipole.
            % This info is stored in dipole_model.transfer_matrix. However,
            % this is a static picture, it does not consider the posibility
            % of ET from acceptor to donor. Thus,we still have to
            % calculate how the emission of the system is after all the
            % excitation energy has decayed via fluorescence - steady state
            % emission after EET.
            [ obj.em_after_et, obj.et_steps, obj.res_ener] = get_em_after_ET( ET_matrix, ind_int );
%             disp(['ET steps: ' num2str(et_steps)])
%             disp(['Residual energy: ' num2str(res_energy)])

            %%%%%%%% EMISSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Polarizer = get_em_pol( em_angles_r );

            % Once defined the observation axis (Polarizer), meaning the
            % orientation of the emission analyser in space, we can
            % calculate the probability of observing light from each dipole

            % Emission depends on the cos(angle between dipole and analizer)^2
            % to calculate the angle between the dipole and the analizer I use:
            % cos(theta) = (eal(a.b) / ( ||a|| ||b|| )
            % I divide by the magnitude of the dipole to get out of the calculation
            % only the angle. No magnitude for the analyser is needed as I set it to
            % one by definition [1 0 0] [0 1 0] or [0 0 1]

            em_p = ((Polarizer * dip_ori)./ (ones(length(em_angles_r),1) * dipole_magnitude)).^2;
            % em_p = ((Polarizer * dip_ori)./dipole_magnitude ).^2;
            assert ( max (em_p(:)) < 1.0001,...
                                          'look at emission probabilities')

            %%%%%%% INTENSITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Intensity is defined as:
            %  abs_p                           *  em_p                                =
            % ( E . d(i) )^2                   * ( E . d(i) / (  ||E|| ||d(i)|| ) )^2 =
            % ||E||^2 ||d(i)||^2 cos(theta)^2  *  cos(theta)^2;
            % where theta is the angle between the electric field E and the
            % dipole d(i) Intensity I_12 is observed from axis 1 and
            % excited though axis 2. this is the intensity for each dipole

            % do note that I only excite the dipole of interest therefore
            abs_int = zeros(size(abs_p));
            abs_int(:,ind_int) = abs_p(:,ind_int);
            % Thus, as mentioned before I do not need to calculate the
            % steady state transfer matrix for all dipole, only for the
            % dipole of interest, as it is the only one that has any energy
            % to redistribute.
            tm      = zeros(size(ET_matrix));
            tm(ind_int,:) = obj.em_after_et;
%             tm(ind_int,:) = round(obj.em_after_et); % small test to get 0.4 at
%             large distances
            % calculate how the excited state is redistributed after EET.
            ex_ET = abs_int * tm;

            I  = ex_ET(:,:) * em_p(:,:)';          % Intensity of all dipoles

            obj.ex_angles_rad = ex_angles_r;
            obj.em_angles_rad = em_angles_r;
            obj.I_ex_em = I;

        end

        function display_portrait(obj)
            figure()
            contourf(obj.ex_angles_rad, obj.em_angles_rad, obj.I_ex_em')
            colorbar
            xlabel('Excitation angle')
            ylabel('Emission angle')
            axis equal
            title('Polarization portrait')
        end
        
        function obj = modify_intensity(obj, new_int)
            obj.I_ex_em = new_int;
        end



    end

end
