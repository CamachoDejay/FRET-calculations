classdef homoFRET
    %HOMOFRET contains all the properties that are relevant to the FRET
    %process, which we will use to calculate the FRET effciency
    
    properties
        absSpec
        emiSpec
        wavelength
        extCoefA {mustBePositive}
        lifetimeD{mustBePositive}
        quantumYieldD{mustBePositive}
        refIndex {mustBePositive}
        J
        R0 {mustBePositive}
    end
    
    methods
        function obj = homoFRET(absA,emD,wavelength,extCoefA,refIndex,QYD,ltD)
            %HOMOFRET Construct an instance of this class, homoFRET
            %contains all the properties that are relevant to the FRET
            %process, which we will use to calculate the FRET effciency
           
            %   wavelength: wavelength axis
            assert(iscolumn(absA),'Absorption of the acceptor must be a colum vector')
            assert(iscolumn(emD),'Emission of the donor must be a colum vector')
            assert(iscolumn(wavelength),'Wavelength axis must be a colum vector')
            assert(and(length(absA)==length(wavelength),length(emD)==length(wavelength)),...
                   'abs and emission spec must be on the same wavelength axis');
            % extCoefA must be a positive real number
            obj.extCoefA = extCoefA;
            % refIndex must be a positive real number
            obj.refIndex = refIndex;
            % quantumYieldD must be a positive real number 0-1
            obj.quantumYieldD = QYD;
            % lifetimeD must be a positive real number 0-1
            obj.lifetimeD = ltD;
            
            % delta in the wavelength axis
            dWavelength = unique(diff(wavelength)); %[nm]
            assert(length(dWavelength)==1,'wavelength axis must be equaly spaced, or code modified')
            % maximun normalized absorption of acceptor
            norAbsA = absA ./ max(absA);
            % epsilon from Forster equation
            epsilon      = norAbsA .*extCoefA; % [mol^-1 cm^-1]
            % sum normalized fluorescence spectra
            fluoNorm    = (emD.*dWavelength) ./ (sum(emD.*dWavelength)); %[dimensionless]
            %   Spectral overlap parameter - J in Forster theory
            obj.J   = sum(fluoNorm .* epsilon .* (wavelength.^4)); % [mol^-1 cm^-1 nm^4]
           
        end
        
        function R = getFRadius(obj)
            %%%%%%% Physical constants %%%%%%%
            %avogadros number                %
            N_a = 6.022e23; %[mol^-1]        %
            %index of refraction             %
            n   = obj.refIndex;%  1.4 aqueous sol
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % contant for Forster integral calculation
            f1 = 9000 * log(10) / (128 * pi^5 * N_a); % [mol cm^3]

            %%%%%%% Chromophore properties %%%
            %QY of the donor                 %
            Q_d = obj.quantumYieldD;   % 0.68 for GFP emerald
            % Spectral overlap
            Jfactor = obj.J; % [mol^-1 cm^-1 nm^4]
            f2 = 1e14; % [nm^2 * cm^-2] This factor takes into account the units of J.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % orientation factor
            k2 = 2/3;
            % Forster radius
            R = (f1 * f2 * ((Q_d * Jfactor * k2) ./ (n^4)))^(1/6);% [nm];
            
        end
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

