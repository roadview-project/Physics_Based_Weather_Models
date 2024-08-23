% Time taken to apply snow (seconds) for 20 frames: 608.6816
% Frame Rate (per second) = 0.032858 (1 frame every 30 second)

%%% For single input and output:
% Inputs and Output Definition (s): 7.24e-05
% Setting Weather Paramaters (s): 4.8e-05
% Setting LiDAR Paramaters (s): 0.0001471
% Obtaining Energy Threshold (s): 0.0010715
% Obtaining Rain DSD (s): 0.0005694
% Obtaining Extinction Coefficient Alpha (s): 28.1954 %% Can potentially reduce this by appriximating Q_EXT
% Calculating Emission Angles (s): 0.21276
% Setting Zero-Range Emission Angles (s): 0.24366
% Appling Snow Model (s): 2.1042
% Saving Simuated PCs (s): 5.26e-05


function [PC_SnowSIM, Intensity_SnowSIM]= LiDAR_SnowModel_Ouster(Clear_PC, Clear_Intensity, Snow_Rate)
 
    %% Define Number of Input pointclouds and Output pointclouds

    pointCloudNum = 1; % Number of clear pointclouds inputted
    Repetitions = 1; % Number of simulated Snowy pointclouds to be outputted

    PC = Clear_PC; % Clear PC Coordinates (Cartesian) - stored in cell
    Intensity_value = Clear_Intensity; % Clear PC Intensity Values - stored in cell
    
    %% Set Parameters of the weather
    Snow_Intensity = Snow_Rate; % mm/h (snow rates used in validation = 0.25 and 0.835 mm/h)
    
    % Maximum and Minumum diamater size for snow flakes - taken from Gunn et
    % al: https://journals.ametsoc.org/view/journals/atsc/15/5/1520-0469_1958_015_0452_tdwsoa_2_0_co_2.xml
    min_FlakeSize = 0.1; % mm 
    max_FlakeSize = 4; % mm
    D_c = 1; % flakelet diameter of max probability (mm) - i.e. most probable rain flake diameter (see section 4.5 in: https://www.mdpi.com/1424-8220/23/15/6891)
    
    % Ice refractive index for 865 nm wavelength was taken from this link: https://refractiveindex.info/?shelf=3d&book=crystals&page=ice#:~:text=Ice%20%28H2O%29%20Wavelength%3A%20%C2%B5m%20%284.430E-002%E2%80%932.000E%2B006%29%20Complex%20refractive%20index,Refractive%20index%20%5B%20i%20%5D%20n%20%3D%201.3098
    % reference = . G. Warren, and R. E. Brandt. Optical constants of ice from the ultraviolet to the microwave: A revised compilation, J. Geophys. Res. 113, D14220 (2008) Optical constants available online at https://atmos.washington.edu/ice_optical_constants/
    Ice_RefractIndex = complex(1.3038,2.4000*(10^-7)); % This needs to be the complex refractive index of the flake
   
    %% Set Parameters of the LiDAR

    % Datasheet link - https://data.ouster.io/downloads/datasheets/datasheet-revd-v2p0-os1.pdf
    BeamRadius = 0.0095/2; % m
    LiDAR_Max_Range = 120; % m
    wavelength = 865*(10^-9); % m (The wavelength of Ouster LiDAR is 865nm)
    pulse_width = 10*(10^-15); % s (Estimated by Valentina and Jonathan)
    Receiver_Area = 1*(10^-6); % m^2 (Estimated by Valentina and Jonathan - original estimate was 1mm^2 = 1*10^-6m^2)
    Channels = 128; % number
    Horizontal_FOV = 360*(pi/180); % radians (360 degrees)
    min_az_angle = -Horizontal_FOV/2; % radians (Minimum Horizontal FOV)
    max_az_angle = Horizontal_FOV/2; % radians (Maximum Horizontal FOV)
    Vertical_FOV = 45*(pi/180); % radians (45 degrees)
    min_z_angle = -22.5*(pi/180); % radians (Minimum Vertcal FOV)
    max_z_angle = 22.5*(pi/180); % radians (Maximum Vertcal FOV)
    Horizontal_Res = Horizontal_FOV/2048; % radians (angular resolution)
    Vertical_Res = Vertical_FOV/Channels; % radians (angular resolution)
    
    Planck = 6.62607015*(10^-34); % js
    speedOfLight = 299792458 ; % m/s (Actual Value = 299792458 m/s, approximation = 3*10^8 m/s)
    
    Energy_value_Max = (2^16)*Planck*(speedOfLight/wavelength); % 16 bit unsigned int - Ouster (Max_energy_value = max intensity value*planck*speed_of_light)
    
    %Definitions for snow particle reflectivity (rho_s) and object
    %reflectivivity (rho_0) are taken from hadher et al: https://arxiv.org/pdf/2203.15118.pdf
    rho_s = 0.9; %snow particle reflectivity (rho_s) (see hadher et al: https://arxiv.org/pdf/2203.15118.pdf)
    rho_0 = (1*10^-6)/pi; %object reflectivity (rho_s) (see hadher et al: https://arxiv.org/pdf/2203.15118.pdf)

    %% Find Energy of LiDAR Returns

    % Ouster LiDAR declares intensity as number of photons up to the maximum of
    % 2^16 photons. Hence we simply mutliply the intensity value by the energy
    % of a photon to get the received energy
    
    Energy_value = cell(pointCloudNum,1);
    
    for e=1:pointCloudNum
        Energy_value{e} = Intensity_value{e}*Planck*(speedOfLight/wavelength);
    end
    
    Energy_Threshold = 0.000005*Energy_value_Max; % currently set to 0.0005% of maximum emitted energy - The reason why the threshold is very low is because most of the intensity values in the PCs are between 0 - 15, (almost all of them roughly 240000/260000)

    %% Distribution function

    % Define the flakesize distribution for the snow flake diamater size
    % distribution using the general form (exponential) provided in 
    % Rasshofer and the paramater values for N_0 and lambda.

    % Convert snow Intensity to Rain intensity equailivant 
    average_snow_density = 0.178*(D_c)^-0.992; % See Equation 7 from https://journals.ametsoc.org/view/journals/apme/46/5/jam2489.1.xml#i1558-8432-46-5-634-f13
    Rain_Intensity = ((Snow_Intensity/(487*average_snow_density*(D_c*10^-3)*1.6)).^3).^(1/4); % Using equation 2 from: https://light.princeton.edu/wp-content/uploads/2022/04/LiDAR_Snowfall_Sim_supplementary_CVPR_2022_compressed.pdf)

    % Distribution paramaters defined in: https://ieeexplore.ieee.org/document/9860868
    N_0 = (7.6*10^3)*Rain_Intensity^(-0.87); % See Gunn and Marshall in Table 3 of Rasshofer et al: https://ars.copernicus.org/articles/9/49/2011/ars-9-49-2011.pdf
    lambda = 2.55*Rain_Intensity^(-0.48); % See lambda Equation in https://light.princeton.edu/wp-content/uploads/2022/04/LiDAR_Snowfall_Sim_supplementary_CVPR_2022_compressed.pdf
    
    % snow flake sizes 'D':
    D = linspace(min_FlakeSize,max_FlakeSize,199);

    % Snow size distribution 
    % (REFERENCES: Exponential from Rasshofer and https://journals.ametsoc.org/view/journals/apme/46/5/jam2489.1.xml#i1558-8432-46-5-634-f13)
    % Hahner et al also use exponential distribution in their snow model
    N_D = N_0.*(exp(-lambda.*D)); 
    
    % Normalizing N_D and multiplying it by 600
    N_D = (N_D./max(N_D)).*600;
    
    % Obtaining the total number of snoflakes per unit area (m^2)
    int_sum = trapz(D,N_D); % Integrate DSD (N(D)) between min and max snowflake diamater size to get mean number of snowflakes per unit volume (per unit m^2)

    %% Extinction coefficient

    % Calulation of the extinction coefficient alpha. The equation must be
    % changed to have the correct DSD. This is in dB/km so must be converted to
    % dB/m. x_EXT is the size parameter. 
    x_EXT = (pi/wavelength)*(D*(10^-3)); %multiply D by 10^-3 to convert from mm to m
    m_EXT = real(Ice_RefractIndex);
    N_MAX = round((x_EXT + (4.*(x_EXT.^(1/3))) + 10), 0);
    
    Q_EXT = zeros(1,length(x_EXT));
    for idx = 1:length(x_EXT)
        Q_EXT(idx) = (2/(x_EXT(idx)^2))*C_EXT(x_EXT(idx),m_EXT,N_MAX(idx));
    end

    h = D(2) - D(1); % obtaining step size

    % Using Trapezium Rule to Calculate Extinction Coefficient:
    alpha_sum = 0;
    for D_index = 1:length(D)
        if D_index == 1
            alpha_sum = alpha_sum + ((D(D_index)^2) * Q_EXT(D_index) * N_D(D_index));
        elseif D_index == length(D)
            alpha_sum = alpha_sum + ((D(D_index)^2) * Q_EXT(D_index) * N_D(D_index));
        else
            alpha_sum = alpha_sum + (2*((D(D_index)^2) * Q_EXT(D_index) * N_D(D_index)));
        end
    end
    
    % Calculating alpha (Extinction Coefficient)
    alpha = (1/2)*h*alpha_sum;
    alpha = (pi/8)*alpha; % dB/km
    
    % Convert alpha to dB/m
    alpha = alpha/1000; %dB/m
 
    %% Calculate angles

    % Calculation of the angles of the points returned in azimuth and elevation 
    % in the clean PC. The x axis is set to be 0 degrees with anticlockwise 
    % negative up to -pi and with clockwise positive up to pi. If a point has 
    % zero range, the angles are set to 1000. angle = 2000 is set as an error 
    % catcher.
    
    Angles = cell(pointCloudNum,1);
    Range = cell(pointCloudNum,1);
    
    for pc_angles = 1:pointCloudNum
        Angles{pc_angles} = zeros(length(PC{pc_angles}(:,1)),2);
        Range{pc_angles} = zeros(length(PC{pc_angles}(:,1)),1);
        
        for k=1:length(Angles{pc_angles}(:,1))
            Range{pc_angles}(k,1) = sqrt(PC{pc_angles}(k,1)^2 + PC{pc_angles}(k,2)^2 + PC{pc_angles}(k,3)^2);
            if Range{pc_angles}(k,1) > 0
                % Calculating azimuth angles:
                if PC{pc_angles}(k,1) > 0 && PC{pc_angles}(k,2) >= 0
                    % If x and y are positive, angle is between 0 and pi/2
                    Angles{pc_angles}(k,1) = atan(PC{pc_angles}(k,2)/PC{pc_angles}(k,1));
                elseif PC{pc_angles}(k,1) < 0 && PC{pc_angles}(k,2) >= 0
                    % If x is negative and y is positive, angle is between pi/2 and pi
                    Angles{pc_angles}(k,1) = mod(atan(PC{pc_angles}(k,2)/PC{pc_angles}(k,1)),pi);
                elseif PC{pc_angles}(k,1) < 0 && PC{pc_angles}(k,2) < 0
                    % If x and y are negative, angle is between -pi/2 and -pi
                    Angles{pc_angles}(k,1) = atan(PC{pc_angles}(k,2)/PC{pc_angles}(k,1)) - pi;
                elseif PC{pc_angles}(k,1) > 0 && PC{pc_angles}(k,2) < 0
                    % If x is positive and y is negative, angle is between 0 and -pi/2
                    Angles{pc_angles}(k,1) = atan(PC{pc_angles}(k,2)/PC{pc_angles}(k,1));
                elseif PC{pc_angles}(k,1) == 0 && PC{pc_angles}(k,2) > 0
                    Angles{pc_angles}(k,1) = pi/2;
                elseif PC{pc_angles}(k,1) == 0 && PC{pc_angles}(k,2) < 0
                    Angles{pc_angles}(k,1) = -pi/2;
                else
                    Angles{pc_angles}(k,1) = 2000;
                end
                % Calculating elevation angles:
                Angles{pc_angles}(k,2) = asin(PC{pc_angles}(k,3)/Range{pc_angles}(k,1));
        
            else
                Angles{pc_angles}(k,1) = 1000;
                Angles{pc_angles}(k,2) = 1000;
            end
        end
    end

    %% Set all zero range points to an angle

    % As we know we have points with angle set to 1000 we can use this to find
    % the points with zero range. These points have no angle we can extract so
    % we need to predict an emission angle. This is such that we can still place
    % a snowflake point even if no point was present originally.
    
    for ang_k2 = 1:pointCloudNum
        for k2 = 1:length(Angles{ang_k2}(:,1))
            if Angles{ang_k2}(k2,1) == 1000
                if k2 == 1
                    Angles{ang_k2}(k2,1) = 0;
                    Angles{ang_k2}(k2,2) = 0;
                elseif k2 == 2
                    Angles{ang_k2}(k2,1) = Angles{ang_k2}(k2-1,1) - Horizontal_Res;
                    Angles{ang_k2}(k2,2) = Angles{ang_k2}(k2-1,2) - Vertical_Res;
                else
                    if Angles{ang_k2}(k2-1,1) + (Angles{ang_k2}(k2-1,1)-Angles{ang_k2}(k2-2,1)) > max_az_angle || ...
                            Angles{ang_k2}(k2-1,1) + (Angles{ang_k2}(k2-1,1)-Angles{ang_k2}(k2-2,1)) < min_az_angle
                        Angles{ang_k2}(k2,1) = Angles{ang_k2}(k2-1,1) - (Angles{ang_k2}(k2-1,1)-Angles{ang_k2}(k2-2,1));
                    else
                        Angles{ang_k2}(k2,1) = Angles{ang_k2}(k2-1,1) + (Angles{ang_k2}(k2-1,1)-Angles{ang_k2}(k2-2,1));
                    end
    
                    if Angles{ang_k2}(k2-1,2) + (Angles{ang_k2}(k2-1,2)-Angles{ang_k2}(k2-2,2)) > max_z_angle || ...
                            Angles{ang_k2}(k2-1,2) + (Angles{ang_k2}(k2-1,2)-Angles{ang_k2}(k2-2,2)) < min_z_angle
                        Angles{ang_k2}(k2,2) = Angles{ang_k2}(k2-1,2) - (Angles{ang_k2}(k2-1,2)-Angles{ang_k2}(k2-2,2));
                    else
                        Angles{ang_k2}(k2,2) = Angles{ang_k2}(k2-1,2) + (Angles{ang_k2}(k2-1,2)-Angles{ang_k2}(k2-2,2));
                    end
                end
            end
        end
    end

    %% Loop through points and change points for snow

    % These variables are to store the points in the snow model point cloud
    PC_SnowAdjusted = cell(pointCloudNum,Repetitions);
    Intensity_value_SnowAdjusted = cell(pointCloudNum,Repetitions);
    
    for pc_1 = 1:pointCloudNum % Loop through the point clouds (inputted pointclouds)
        for rep=1:Repetitions % Loop through the repetitions (number of output pointclouds)
            
            % These variables are to store the points in the snow model point cloud
            PC_SnowAdjusted = zeros(length(PC{pc_1}(:,1)),length(PC{pc_1}(1,:)));
            Intensity_value_SnowAdjusted = zeros(length(PC{pc_1}(:,1)),1);

            for i=1:length(PC{pc_1}(:,1)) % Loop through the points in the original clean PC
                DistanceToObject = Range{pc_1}(i,1); % Set the DistanceToObject variable to the Range for that point (DistanceToObject is used in the unreal code)

                % If Range is greater than zero, we adjust intensity of
                % object and intensity of flake. If not, we just calculate
                % intensity of snowflake.
                if DistanceToObject > 0 
        
                    LiDARYaw = Angles{pc_1}(i,1); % Set LiDARYaw to the azimuth angle (LiDARYaw is used in unreal code)
                    LiDARPitch = Angles{pc_1}(i,2); % Set LiDARPitch to the elevation angle (LiDARPitch is used in unreal code)
        
                    VBeam = 3.1416 * (BeamRadius^2)*DistanceToObject; % Calulate volume of beam using pi*(r^2)*h
                    Mean_V = VBeam*int_sum; % Multiply the above to get mean number of flakes in the beam
                    poisson_num = poissrnd(Mean_V); % Sample the sample of flakes from poisson distribution using the mean
                    if poisson_num > 0 % If at least one flake, calculate the point intensities
                        % This is a term in the literature for calculating
                        % power. The difference here is we are using energy
                        % rather than power so the equation for power
                        % received becomes energy received
                        C_A_Peak_Energy = Energy_value{pc_1}(i)*((DistanceToObject^2)/rho_0); % Equation taken from Hahner et al: https://arxiv.org/pdf/2203.15118.pdf
                        
                        % Calculate flakeDistances based on number of flakes
                        flakeDistances = ((rand(poisson_num,1) * 90) + 10) * 0.01 * DistanceToObject;
                        
                        % Start sampling flake sizes
                        probabilities = N_D./ sum(N_D);
                        % Number of samples you want to generate
                        num_samples = length(flakeDistances);
                        % Sample flake sizes using known frequencies of D
                        flakeSizes = randsample(D, num_samples, true, probabilities);
                      
                        % calculate flake energy equation (equation 11 in
                        % https://arxiv.org/pdf/2203.15118.pdf) -
                        % multiply this equation by the pulse-width to
                        % convert it from power to energy
                        flakeEnergies = pulse_width .* ((C_A_Peak_Energy*rho_s)./(flakeDistances.^2)) .* (sin((pi.*(DistanceToObject - flakeDistances))./(speedOfLight*pulse_width))).^2;
                        
                        % Adjusting returned flake energies based on ratio
                        % between snow flake size and beam diamater
                        ratios = (flakeSizes./(BeamRadius*1000*2)).^2;
                        ratios(ratios > 1) = 1;
                        flakeEnergies = flakeEnergies .* ratios';

                        % convert any Nan values to 0
                        flakeEnergies(isnan(flakeEnergies))=0;
                        
                        % Convert inf values to max numerical value
                        max_value = max(flakeEnergies(~isinf(flakeEnergies)));
                        flakeEnergies(isinf(flakeEnergies)) = max_value;
                        
                        % Find the index of the flake with the maximum
                        % energy and set Energy of soft target to max
                        % energy and Dist_to_flake to distance of flake with
                        % max energy
                        [flakeMaxEnergy, idx] = max(flakeEnergies);
                        Dist_to_flake = flakeDistances(idx);

                        % Calculate returned energy of the flake and the 
                        ReturnedEnergyObj = Energy_value{pc_1}(i)*exp(-2*alpha*DistanceToObject); %Reduce object energy based on rainy environment
                        ReturnedEnergyFlake = flakeMaxEnergy;

                    else % If no snow flakes, energy from object remains the same, energy from flake is 0
                        ReturnedEnergyObj = Energy_value{pc_1}(i);
                        ReturnedEnergyFlake = 0;
                    end

                    % Changing the position of the PC point based on
                    % returned energies
                    if ReturnedEnergyFlake > ReturnedEnergyObj
                        % If energy of snowflake is higher than the
                        % returned energy of the original point then
                        % change the position of that point to the
                        % position of the snow flake in the beam
 
                        xsnow = cos(LiDARPitch)*cos(LiDARYaw)*Dist_to_flake;
                        ysnow = cos(LiDARPitch)*sin(LiDARYaw)*Dist_to_flake;
                        zsnow = sin(LiDARPitch)*Dist_to_flake;
                
                        PC_SnowAdjusted(i,:) = [xsnow,ysnow,zsnow];
                        % Convert energy to intensity
                        Intensity_value_SnowAdjusted(i) = (ReturnedEnergyFlake/pulse_width)/Receiver_Area;
                        
                    elseif ReturnedEnergyFlake <= ReturnedEnergyObj
                        if ReturnedEnergyObj < Energy_Threshold
                            % If energy threshold is higher than 
                            % object energy, remove the point
                            % by setting to zero
                            PC_SnowAdjusted(i,:) = [0,0,0];
                            Intensity_value_SnowAdjusted(i) = 0;
                        elseif ReturnedEnergyObj >= Energy_Threshold
                            % If energy of original point is higher than the
                            % returned energy of the snowflake then
                            % keep the position the same
                            PC_SnowAdjusted(i,:) = PC{pc_1}(i,:);
                            % Convert energy to intensity
                            Intensity_value_SnowAdjusted(i) = (ReturnedEnergyObj/pulse_width)/Receiver_Area;
                        end
                    end
                
                % Modifying Pointclouds points set to [0,0,0] initially  
                else
                    
                    LiDARYaw = Angles{pc_1}(i,1); % Set LiDARYaw to the azimuth angle (LiDARYaw is used in unreal code)
                    LiDARPitch = Angles{pc_1}(i,2); % Set LiDARPitch to the azimuth angle (LiDARPitch is used in unreal code)
        
        
                    VBeam = 3.1416 * (BeamRadius^2)*DistanceToObject; %Calulate volume of beam using pi*(r^2)*h
                    Mean_V = VBeam*int_sum; % Multiply the above to get mean number of flakes in the beam
                    poisson_num = poissrnd(Mean_V); % Sample the sample of flakes from poisson distribution using the mean
        
                    if poisson_num > 0 % If more than one flake, calculate the flake intensities
                        % This is a term in the literature for calculating
                        % power. The difference here is we are using energy
                        % rather than power so the equation for power
                        % received becomes energy received
                        C_A_Peak_Energy = Energy_value{pc_1}(i)*((DistanceToObject^2)/rho_0); % equation taken from Hahner et al: https://arxiv.org/pdf/2203.15118.pdf
                        
                        % Calculate flakeDistances based on number of flakes
                        flakeDistances = ((rand(poisson_num,1) * 90) + 10) * 0.01 * DistanceToObject;
                        
                        % Start sampling flake sizes
                        probabilities = N_D./ sum(N_D);
                        % Number of samples you want to generate
                        num_samples = length(flakeDistances);
                        % Sample flake sizes using known frequencies of D
                        flakeSizes = randsample(D, num_samples, true, probabilities);
                      
                        % calculate flake energy equation (equation 11 in
                        % https://arxiv.org/pdf/2203.15118.pdf)
                        flakeEnergies = ((C_A_Peak_Energy*rho_s)./(flakeDistances.^2)) .* (sin((pi.*(DistanceToObject - flakeDistances))./(speedOfLight*pulse_width))).^2;
                        
                        % Adjusting returned flake energies based on ratio
                        % between snow flake size and beam diamater
                        ratios = (flakeSizes./(BeamRadius*1000*2)).^2;
                        ratios(ratios > 1) = 1;
                        flakeEnergies = flakeEnergies .* ratios';


                        % convert any Nan values to 0
                        flakeEnergies(isnan(flakeEnergies))=0;
                        
                        % Convert inf values to max numerical value
                        max_value = max(flakeEnergies(~isinf(flakeEnergies)));
                        flakeEnergies(isinf(flakeEnergies)) = max_value;
                        
                        % Find the index of the flake with the maximum
                        % energy and set Energy of soft target to max
                        % energy and dist_to_flake to distance of flake with
                        % max energy
                        [flakeMaxEnergy, idx] = max(flakeEnergies);
                        Dist_to_flake = flakeDistances(idx);

                        ReturnedEnergyFlake = flakeMaxEnergy;
                    else % If no snow flakes/flakes, energy from flake is 0
                        ReturnedEnergyFlake = 0;
                    end
        
                    
                
                    if ReturnedEnergyFlake < Energy_Threshold
                        % If energy of snowflake is less than 
                        % threshold, position the point is set to
                        % [0,0,0] and intensity = 0 (i.e. remove this
                        % point)
                        PC_SnowAdjusted(i,:) = [0,0,0];
                        Intensity_value_SnowAdjusted(i) = 0;
                        
                    else 
                        % If energy of snowflake is higher than threshold, position the point
                        % to represent the flake and update the intensity value
                        xsnow = cos(LiDARPitch)*cos(LiDARYaw)*Dist_to_flake;
                        ysnow = cos(LiDARPitch)*sin(LiDARYaw)*Dist_to_flake;
                        zsnow = sin(LiDARPitch)*Dist_to_flake;
                
                        PC_SnowAdjusted(i,:) = [xsnow,ysnow,zsnow];
                        % Convert energy to intensity
                        Intensity_value_SnowAdjusted(i) = (ReturnedEnergyFlake/pulse_width)/Receiver_Area;
                    end
                end
            end
        end
    end
   
    %% Export Variables

    PC_SnowSIM = PC_SnowAdjusted;
    Intensity_SnowSIM = Intensity_value_SnowAdjusted;

end