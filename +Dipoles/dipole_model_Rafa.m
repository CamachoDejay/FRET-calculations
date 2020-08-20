classdef dipole_model
    % DIPOLE_MODEL Object that contains the positions, orientations and
    % transfer matrix of all dipoles in the model system. The list of
    % positions contains a central dipole in the origin surrounded by a
    % large number of buffer dipoles in a cubic lattice. All dipoles are
    % randomly oriented. There is the option of replacing the monomer sites
    % by a dimer.
    %   TODO

    properties (SetAccess = protected)
        distance_model
        buffer
        buffer_used
        positions
        orientations
        transfer_matrix
        central_dipole
    end

    methods
% object constructor
        function obj = dipole_model(dist_model, buffer)
            % distance model must be a string: 1D model 2D model 3D model
            if isa(dist_model,'char')
                % input is a string
                switch dist_model
                    case '1D'
                        disp('1D model')
                    case '2D'
                        disp('2D model')
                    case '3D'
                        disp('3D model')
                    otherwise
                        error(['Unknown distance model, '...
                               'input must be 1D 2D or 3D, ',...
                               'you entered ' dist_model])
                end
            else
                error('Distance model must be a string')
            end
            % if dist_model is ok we store it
            obj.distance_model = dist_model;

            % buffer must be a structure with fields: max_n_dipoles,
            % max_size and min_size.
            if isa(buffer,'struct')
                expected_fields = {'max_n_dipoles';'max_size';'min_size'};
                if all(strcmp(sort(fieldnames(buffer)),expected_fields))
                    disp('buffer is ok')
                else
                    error ('buffer fields must be max_n_dipoles, max_size and min_size')
                end
            else
                error(['Buffer must be a structure with ',...
                       'fields max_n_dipoles, max_size and min_size in units of nm'])
            end
            % storing buffer
            obj.buffer = buffer;
            % end of object construction
        end
% methods of the object:

        function obj = get_positions(obj, inter_dist, dimer_prob)
        % get_positions: defines the positions of all dipoles in the system
        % assuming a squared lattice

            % inter_dist is the interchromophorix distance in nm
            assert(isa(inter_dist,'float'),['Inter chromophoric distance '...
                                        'must be a float in units of nm']);

            assert(dimer_prob >= 0 && dimer_prob <= 1,...
                   ['The probability that a site in the cubic lattice'...
                    'contains a dimer must be bwteeen [0-1]']);

            % Due to computational times the max number of dipoles that we
            % can consider is finite
            max_n_dipoles = obj.buffer.max_n_dipoles;

            % init the lattice (cube) dimentions
            cube.x_radius = 0; % [nm]
            cube.y_radius = 0; % [nm]
            cube.z_radius = 0; % [nm]

            % now the dimentions of the lattice are changed depending on
            % the distance_model
            switch obj.distance_model
                case '1D'
                    d_sugested =  (max_n_dipoles-1)* inter_dist;
                    d_sugested = round((d_sugested / 2)*1.1,1);
                    d2use = min([obj.buffer.max_size,d_sugested]);
                    cube.x_radius = d2use;
                case '2D'
                    d_sugested = (((max_n_dipoles*4/pi)-1)^(1/2)) * inter_dist;
                    d_sugested = round((d_sugested / 2)*1.1,1);
                    d2use = min([obj.buffer.max_size,d_sugested]);
                    cube.x_radius = d2use;
                    cube.y_radius = d2use;
                case '3D'
                    d_sugested = (((max_n_dipoles*3*2^3/(4*pi))-1)^(1/3)) * inter_dist;
                    d_sugested = round((d_sugested / 2)*1.1,1);
                    d2use = min([obj.buffer.max_size,d_sugested]);
                    cube.x_radius = d2use;
                    cube.y_radius = d2use;
                    cube.z_radius = d2use;
                otherwise
                    error('Unexpected error')
            end

            % get positions in a cubic lattice: equidistant and centered
            % around 0
            pos = get_pos_cubic_lattice( cube, inter_dist );
            
            % now we replace some of the monomer sites by dimers
            %   distance between the chromophores in the dimer
            dimer_distance   = 3.5; % [nm] GFP barrel size 3x4 nm
            %   modify positions to account for dimers
            pos = set_dimers( pos, dimer_prob, dimer_distance );
            
            % find the distance between the dipoles and the origin (x:0,
            % y:0)
            pos_magnitude = (sum(pos.*pos)).^(.5);

            % remove positions that are larger than the desired radious
            good_pos  = pos_magnitude <= cube.x_radius;
            pos = pos(:,good_pos);
            pos_magnitude = pos_magnitude (good_pos);
            
            % the generation of dimers might displace the central dipole by
            % spliting it into a dimmer. Thus I have to move the reference
            % system to make it 0 0 0 again;
            n_dipoles = size(pos,2);
            [~,min_idx] = min(pos_magnitude);
            pos_cent_dip = pos(:,min_idx);
            pos = pos - repmat(pos_cent_dip, 1,n_dipoles);
            pos_magnitude = (sum(pos.*pos)).^(.5);

            % sort pos by their distance to the center (dipole of interest)
            [pos_magnitude, i] = sort(pos_magnitude);
            pos = pos(:,i);

            if n_dipoles > max_n_dipoles
                % then we have to decrease the buffer size
                approx_buff_size = pos_magnitude(max_n_dipoles);
                pos = pos(:,1:max_n_dipoles);
                pos_magnitude = pos_magnitude(1:max_n_dipoles);

                if approx_buff_size < obj.buffer.min_size
                    warning(['We need to use a buffer of '...
                              num2str(approx_buff_size,3)...
                              '[nm] instead of ' ...
                              num2str(obj.buffer.min_size,3) '[nm]'] )
                end
                % store the actuall bufer used
                obj.buffer_used = approx_buff_size;
            else
                obj.buffer_used = cube.x_radius;

            end
            obj.positions = pos;
            c_dip = pos_magnitude==0;
            assert(c_dip(1),'Unexpected central dipole')
            obj.central_dipole = c_dip;
%             disp(['Number of dipoles: ' num2str(size(pos,2))])

        end

        function obj = get_orientations(obj,sigma)
            % please note that for the calculations done later on it is
            % much better to define the orientation of the dipoles via
            % unitary vectors.
            
            % we first check that we have some positions
            if isempty( obj.positions )
               warning('You must define the positions first! returning empty orientations')
               return
            end
        
            % Its important that dipoles are randomly oriented in the whole
            % sphere. the direction vector of dipole i is given by: vi =
            % v(:,i); in x y z
            n_dipoles = size(obj.positions,2);
            
            switch nargin
                case 1
                    % uniform distribution
                    [ obj.orientations ] = get_random_orient_vec( n_dipoles );
                case 2
                    % random distribution
                    if sigma > 5
                        % this is so flat that uniform dist works best
                        [ obj.orientations ] = get_random_orient_vec( n_dipoles );
                    else
                        % random/gaussian like distribution
                        [ obj.orientations ] = get_orient_vec( n_dipoles, sigma );
                    end
                otherwise
                    error('??')
            end
                    
                

        end

        function obj = get_et_matrix(obj, FRETprops)
            % preconditions:
            assert(isa(FRETprops,'homoFRET'), 'FRETprops must be of class homoFRET')

            % we first check that we have some positions
            if isempty( obj.positions )
             warning('You must define the positions first! returning empty et_matrix')
             return
            end

            % we first check that we have some orientations
            if isempty( obj.orientations )
             warning('You must define the orientations first! returning empty et_matrix')
             return
            end
          
            % calculation of the orientation factor
            [ k2 ] = get_ori_factor( obj.positions, obj.orientations );

            % calculation of the eucli distance between all dipoles in nm
            [ r ] = get_euc_dist_matrix( obj.positions );

            % calculation of the spectral overlap and other physical
            % constants that will be use for the FRET efficiency
            % calculation
            [ k_t ] = get_transfer_rate_matrix( FRETprops, k2, r );

            %transfer efficiency
            [ obj.transfer_matrix ] = get_transfer_prob_matrix( k_t, FRETprops.lifetimeD);


            % check that the porbability transfer matrix is ok
            test   = sum(obj.transfer_matrix,2);
            assert (all(test > (1e5/(1e5+1))),...
                                   'Unexpected problem in transfer matrix')


        end



    end

end
