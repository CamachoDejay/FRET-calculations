classdef twoDmodel
    %TWODMODEL object that contains the data to calculate the 2D sinfle
    %funnel approximation model.
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        M_ex
        P_ex
        M_f
        P_f
        X
        epsilon
        ex_ang_rad
        em_ang_rad
        intensity
    end
    
    methods
        function obj = twoDmodel(mex, pex, mf, pf, X, epsilon, ex_ang_rad, em_ang_rad) % object constructor
            assert(mex >= 0 && mex <= 1, 'mex is defined out of bounds [0,1]')
            obj.M_ex = mex;
            assert(pex >= -pi/2 && pex <= pi/2, 'pex is defined out of bounds [-pi/2, pi/2]')
            obj.P_ex = pex;
            assert(mf >= 0 && mf <= pi, 'mf is defined out of bounds [0,1]')
            obj.M_f = mf;
            assert(pf >= -pi/2 && pf <= pi/2, 'pf is defined out of bounds [-pi/2, pi/2]')
            obj.P_f = pf;
            Xlb = 0;
            Xub = 2*(1+mex)/(1-mex);
            assert(X >= Xlb && X <= Xub, 'X, geometrical ratio defined out of bounds')
            obj.X = X;
            % I want angles to be a one row vector size = [1xn].
            assert(size(ex_ang_rad,1) == 1, 'excitation angles must be a row vector size [1xn]')
            test = all(ex_ang_rad(:)<= pi+0.00001);
            assert(test, 'excitation angles must be in radians and [0 pi]')
            test = all(ex_ang_rad(:)>= 0);
            assert(test, 'excitation angles must be in radians and [0 pi]')
            obj.ex_ang_rad = ex_ang_rad;
            
            % I want angles to be a one row vector size = [1xn].
            assert(size(em_ang_rad,1) == 1, 'emission angles must be a row vector size [1xn]')
            test = all(em_ang_rad(:)<= pi+0.00001);
            assert(test, 'emission angles must be in radians and [0 pi]')
            test = all(em_ang_rad(:)>= 0);
            assert(test, 'emission angles must be in radians and [0 pi]')
            obj.em_ang_rad = em_ang_rad;
            
            if isempty(epsilon)
                % this can happen
                obj.epsilon = [];
            else
                assert(epsilon >=0 && epsilon <= 1, 'epsilon defined out of bounds')
                obj.epsilon = epsilon;
                [Fnoet, Fet] = get_noet_et(obj);
                % to get full model I still have to do:
                F = (1-epsilon)*Fnoet(:) + epsilon*Fet(:);
                F = F./max(F);
                obj.intensity = F;                
            end
            
        end
        
        function [noet, et] = get_noet_et(obj)
            
            % information for the geometrical model, which depends on Mex
            % and geometrical ratio X
            alpha = 0.5*acos(0.5*(((obj.X+2)*obj.M_ex)-obj.X));

            % given alpha I can then calculate the orientation of the dipoles used
            % to mimic the excitation of the system.
            fii   = [(obj.P_ex-alpha) obj.P_ex (obj.P_ex+alpha)];

            % get the No ET component for each dipole
            EnNoTT(:,1) =       cos(obj.ex_ang_rad-fii(1,1)).^2  .*cos(obj.em_ang_rad-fii(1)).^2;
            EnNoTT(:,2) = obj.X*cos(obj.ex_ang_rad-fii(1,2)).^2  .*cos(obj.em_ang_rad-fii(2)).^2;
            EnNoTT(:,3) =       cos(obj.ex_ang_rad-fii(1,3)).^2  .*cos(obj.em_ang_rad-fii(3)).^2;

            % (2+X) is a normalization factor
            noet(1,:) = (sum(EnNoTT,2))./(2+obj.X);    

            % factor 0.25 is for normalization
            et = 0.25 * (1+obj.M_ex*cos(2*(obj.ex_ang_rad-obj.P_ex)))  .*(1+obj.M_f*cos(2*(obj.em_ang_rad-obj.P_f-obj.P_ex)));
            % note that the components come already normalized.
            
        end
    end
    
end

