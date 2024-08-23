% Time taken to apply fog (seconds) for 20 frames: 357.9829 
% Frame Rate (per second) = 0.055869 (1 frame every 18 seconds)

%%% For single input and output:
% Inputs and Output Definition (s): 5.35e-05
% Setting Weather Paramaters (s): 4.13e-05
% Setting LiDAR Paramaters (s): 0.0001272
% Obtaining Energy Threshold (s): 0.0007662
% Obtaining Rain DSD (s): 0.0022265
% Obtaining Extinction Coefficient Alpha (s): 0.0002612
% Obtaining BackScatter Coefficient Beta (s): 0.0014437
% Calculating Emission Angles (s): 0.19376
% Setting Zero-Range Emission Angles (s): 0.28374
% Appling Fog Model (s): 17.4208
% Saving Simuated PCs (s): 2.51e-05

%% Fog Model Function (Takes cell containing single PC coordinates and cell containing single PC intensity values as input + fog rate)
function [PC_FogSIM, Intensity_FogSIM] = LiDAR_FogModel_Ouster(Clear_PC, Clear_Intensity, Fog_Visibility_Distance)
    %% Define Number of Input pointclouds and Output pointclouds

    pointCloudNum = 1; % Number of clear pointclouds inputted
    Repetitions = 1; % Number of simulated Foggy pointclouds to be outputted

    PC = Clear_PC;
    Intensity_value = Clear_Intensity;
    
    %% Set Parameters of the weather

    Fog_visibility = Fog_Visibility_Distance; % m
    min_FogDrop_size = 0.2; % micro-m 
    max_FogDrop_size = 20; % micro-m
    D_c = 10; % droplet diameter of max probability Micrometers
    Water_RefractIndex = complex(1.328,4.86*(10^-7)); % This needs to be the complex refractive index of the drop

    %% Set Parameters of the LiDAR

    BeamRadius = 0.95/2; %cm
    LiDAR_Max_Range = 120; %m
    wavelength = 865*(10^-9); % The wavelength of Ouster LiDAR is 865nm
    pulse_width= 10*(10^-15); %s
    Receiver_Area = 1*(10^-6);
    Channels = 128;
    Horizontal_FOV = 360*(pi/180);
    min_az_angle = -Horizontal_FOV/2;
    max_az_angle = Horizontal_FOV/2;
    Vertical_FOV = 45*(pi/180);
    min_z_angle = -22.5*(pi/180);
    max_z_angle = 22.5*(pi/180);
    Horizontal_Res = Horizontal_FOV/2048;
    Vertical_Res = Vertical_FOV/Channels;
    
    Planck = 6.62607015*(10^-34);
    speedOfLight = 3*(10^8);
    
    Energy_value_Max = (2^16)*Planck*(speedOfLight/wavelength); %16 bit unsigned int - Ouster
    
    beta_zero = (1*10^-6)/pi; % Differential reflectivity (set to value found in Hahner et al Fog model paper: https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/517119/Hahner_Fog_Simulation_on_Real_LiDAR_Point_Clouds_for_3D_Object_ICCV_2021_paper.pdf)
    
    %% Find Energy of LiDAR Returns

    %Ouster LiDAR declares intensity as number of photons up to the maximum of
    %2^16 photons. Hence we simply mutliply the intensity value by the energy
    %of a photon to get the received energy

    Energy_value = cell(pointCloudNum,1);
    
    for e=1:pointCloudNum
        Energy_value{e} = Intensity_value{e}*Planck*(speedOfLight/wavelength);
    end

    Energy_Threshold = (3.7*10^-5)*Energy_value_Max; % 0.0037% of emitted energy

    %% Distribution function

    %Define the dropsize distribution for the fog as the Gamma distribution. We set gamma to
    %a value and then along with setting r_e we can calculate the otehr needed parameters.
    %Then we create a distribution object for use in sampling the size of the drops.
    %To understand the equations used please follow the documentation provided
    %as this makes it much clearer

    % Advection Fog
    Gamma = 1.0;
    r_e = 10; %micro-m effective radius
    LWC = 36.489*(Fog_visibility^-1.1407); %Liquid Water Content
    
    h = 0.1;
    D = min_FogDrop_size:h:max_FogDrop_size;

    alpha_fog_numerator = 3.*Gamma.*((D_c/2)^Gamma);
    alpha_fog_denominator = ((D./2).^Gamma).*r_e - (Gamma.*((D_c/2)^Gamma).*(D./2));
    alpha_fog = alpha_fog_numerator./alpha_fog_denominator;
    
    gamma_Caps = gamma((alpha_fog+1)./Gamma);
    
    rho_numerator = (3*(10^6)).*LWC.*(((alpha_fog+3)./r_e).^(alpha_fog+4)).*gamma((alpha_fog+1)./Gamma);
    rho_denominator = 4*pi.*gamma(alpha_fog+4).*Gamma.*(((alpha_fog) ./ (Gamma*((D_c/2)^Gamma))).^((alpha_fog+1)/Gamma));
    rho = rho_numerator./rho_denominator;
    
    b = alpha_fog./(Gamma*((D_c/2)^Gamma));
    
    ND_numerator = Gamma.*rho.*(b.^((alpha_fog+1)./Gamma)); 
    N_D = (ND_numerator./gamma_Caps).* ((D/2).^alpha_fog).* exp(-b.*((D/2).^Gamma));
    
    int_sum = trapz(D,N_D); % Integrate DSD (N(D)) between min and max drop size to get mean number of drops per unit volume (per unit cm^2)

    %% Extinction coefficient (OLD)
    
    % %Calulation of the extinction coefficient alpha. This is in dB/km so must 
    % %be converted to dB/m. x_EXT is the size parameter. 
    % % tic
    % 
    % 
    % x_EXT = (pi/wavelength)*(D*(10^-6));
    % m_EXT = real(Water_RefractIndex);
    % N_MAX = round((x_EXT + (4.*(x_EXT.^(1/3))) + 10), 0);
    % 
    % Q_EXT = zeros(1,length(x_EXT));
    % for idx = 1:length(x_EXT)
    %     Q_EXT(idx) = (2/(x_EXT(idx)^2))*C_EXT(x_EXT(idx),m_EXT,N_MAX(idx));
    % end
    % 
    % % Q_EXT = (2/(x_EXT^2))*C_EXT(x_EXT,m_EXT,N_MAX)
    % 
    % D = D.*(1*10^-3); % Convert drop size from microns to mm
    % %Trapezium Rule
    % alpha_sum = 0;
    % for D_index = 1:length(D)
    %     if D_index == 1
    %         alpha_sum = alpha_sum + ((D(D_index)^2) * Q_EXT(D_index) * N_D(D_index));
    %     elseif D_index == length(D)
    %         alpha_sum = alpha_sum + ((D(D_index)^2) * Q_EXT(D_index) * N_D(D_index));
    %     else
    %         alpha_sum = alpha_sum + (2*((D(D_index)^2) * Q_EXT(D_index) * N_D(D_index)));
    %     end
    % end
    % 
    % alpha = (1/2)*h*alpha_sum;
    % alpha = (pi/8)*alpha
    % % alpha = ((pi/8)*alpha)/1000 %dB/m
    % % 
    % % % alpha = alpha/2000;
    
    %% Extinction coefficient (NEW)
    % alpha definition below based on this paper: https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/517119/Hahner_Fog_Simulation_on_Real_LiDAR_Point_Clouds_for_3D_Object_ICCV_2021_paper.pdf
    % alpha = 0.06*(50/Fog_visibility);

    MOR = [inf, 600, 300, 150, 100, 50];
    alpha_values = [0, 0.005, 0.01, 0.02, 0.03, 0.06];
    % apply linear regresson (assuming relationship between alpha values and
    % MOR is (alpha_values = a* (1/MOR) + b)
    x = 1./MOR;
    X = [ones(length(x),1) x'];
    coefficients = X\alpha_values';
    
    alpha = coefficients(1) + coefficients(2)*(1/Fog_visibility);

    %% Backscattering coefficient

    %Calulation of the backscatter coefficient beta. This is in dB/km so must 
    %be converted to dB/m

    x_EXT = (pi/wavelength)*(D_c*(10^-6));
    m_EXT = real(Water_RefractIndex);
    N_MAX = x_EXT + (4.*(x_EXT.^(1/2))) + 10;
    Q_BACK = (1/(x_EXT^2))*C_BACK(x_EXT,m_EXT,N_MAX);
    
    %Trapezium Rule
    beta_sum = 0;
    for D_index = 1:length(D)
        if D_index == 1
            beta_sum = beta_sum + ((D(D_index)^2) * Q_BACK * N_D(D_index));
        elseif D_index == length(D)
            beta_sum = beta_sum + ((D(D_index)^2) * Q_BACK * N_D(D_index));
        else
            beta_sum = beta_sum + (2*((D(D_index)^2) * Q_BACK * N_D(D_index)));
        end
    end

    beta = (1/2)*h*beta_sum;
    beta = ((pi/8)*beta)/1000; % dB/m

    %% Calculate angles
    
    %Calculation of the angles of the points returned in azimuth and elevation 
    %in the clean PC. The x axis is set to be 0 degrees with anticlockwise 
    %negative up to -pi and with clockwise positive up to pi. If a point has 
    %zero range, the angles are set to 1000. angle = 2000 is set as an error 
    %catcher.

    Angles = cell(pointCloudNum,1);
    Range = cell(pointCloudNum,1);
    
    for pc_angles = 1:pointCloudNum
        Angles{pc_angles} = zeros(length(PC{pc_angles}(:,1)),2);
        Range{pc_angles} = zeros(length(PC{pc_angles}(:,1)),1);
        
        for k=1:length(Angles{pc_angles}(:,1))
            Range{pc_angles}(k,1) = sqrt(PC{pc_angles}(k,1)^2 + PC{pc_angles}(k,2)^2 + PC{pc_angles}(k,3)^2);
            if Range{pc_angles}(k,1) > 0
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
                Angles{pc_angles}(k,2) = asin(PC{pc_angles}(k,3)/Range{pc_angles}(k,1));
        
            else
                Angles{pc_angles}(k,1) = 1000;
                Angles{pc_angles}(k,2) = 1000;
            end
        end
    end
    
    %% Set all zero range points to an angle

    %As we know we have points with angle set to 1000 we can use this to find
    %the points with zero range. These points have no angle we can extract so
    %we need to predict an emission angle. This is such that we can still place
    %a fog droplet point even if no point was present originally.
    
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
    
    %% Loop through points and change points for Fog

    %These variables are to store the points in the fog model point cloud
    PC_FogAdjusted = cell(pointCloudNum,Repetitions);
    Intensity_value_FogAdjusted = cell(pointCloudNum,Repetitions);
    
    for pc_1 = 1:pointCloudNum %Loop through the point clouds
        for rep=1:Repetitions %Loop through the repetitions
            
            %Initialise the PC and intensity values for this point cloud and
            %repetition
            PC_FogAdjusted = zeros(length(PC{pc_1}(:,1)),length(PC{pc_1}(1,:)));
            Intensity_value_FogAdjusted = zeros(length(PC{pc_1}(:,1)),1);
           
            for i=1:length(PC{pc_1}(:,1)) %Loop through the points in the original clean PC
                DistanceToObject = Range{pc_1}(i,1); %Set the DistanceTo Object variable to the Range for that point (DistanceToObject is used in the unreal code)

                %If Range is greater than zero, we adjust intensity of
                %object and intensity of drop. If not, we just calculate
                %intensity of drop.
                if DistanceToObject > 0
    
                    LiDARYaw = Angles{pc_1}(i,1); %Set LiDARYaw to the azimuth angle (LiDARYaw is used in unreal code)
                    LiDARPitch = Angles{pc_1}(i,2); %Set LiDARPitch to the elevation angle (LiDARPitch is used in unreal code)

                    VBeam = 3.1416 * ((BeamRadius/100)^2)*DistanceToObject; %Calulate volume of beam using pi*(r^2)*h (divide r by 100 to convert from cm to m to obtain volume in m^3) 
                    % int_sum = integral(N_D,min_FogDrop_size,max_FogDrop_size); %Integrate DSD to get mean number of drops per unit volume (per cm^2)
                    Mean_V = VBeam*(int_sum*10000); %Multiply the above to get mean number of drops in the beam (multiply int_sum by 10000 to convert from drops per cm^2 to drops per m^2)
                    poisson_num = poissrnd(Mean_V); %Sample the sample of drops from poisson distribution using the mean
        
                    
                    if poisson_num > 0 %If more than one drop, calculate the drop intensities
                        %This is a term in the literature for calculating
                        %power. The difference here is we are using energy
                        %rather than power so the equation fro power
                        %received becomes energy received
                        C_A_Peak_Energy = Energy_value{pc_1}(i)*((DistanceToObject^2)/beta_zero);
                        
                        % Calculate dropDistances based on number of drops
                        dropDistances = ((rand(poisson_num,1) * 90) + 10) * 0.01 * DistanceToObject;

                        % Sample 100 drop distance (sort dropDistance and then divide dropDistance
                        % array into 100 subarrays and extarct mean od each
                        % sub array)
                        input = sort(dropDistances, 'ascend');
                        N = numel(input);
                        M = 100;
                        % Calculate the number of elements per chunk
                        elements_per_chunk = floor(N / M);
                        % Reshape the input array into a 2D array with the desired number of rows (chunks)
                        reshaped_input = reshape(input(1:elements_per_chunk*M), elements_per_chunk, []);
                        % Calculate the mean along the rows (dimension 1)
                        output = mean(reshaped_input, 1);
                        dropDistances = output';

                        % Define the array of linearly spaced values 't'
                        t = linspace(0, 2*pulse_width, length(dropDistances));
                                               
                        % Create a matrix 'y' using vectorization to
                        % calculate drop energy integral input equation
                        y = C_A_Peak_Energy * beta * (sin((pi.*t)./(2*pulse_width)).^2) .* ...
                            (exp(-2*alpha*(dropDistances - (speedOfLight.*t/2))) ./ ((dropDistances-(speedOfLight.*t/2)).^2));
                                                
                        % integrate y to obtain drop energies for each drop
                        dropEnergies = trapz(t,y,2);

                        % convert any Nan values to 0
                        dropEnergies(isnan(dropEnergies))=0;
                        
                        % Convert inf values to max numerical value
                        max_value = max(dropEnergies(~isinf(dropEnergies)));
                        dropEnergies(isinf(dropEnergies)) = max_value;
                        
                        % Find the index of the snowflake with the maximum
                        % energy and set returned drop energy to max
                        % energy and dist_to_drop to distance of drop with
                        % max energy
                        [Energy_SoftTarget, idx] = max(dropEnergies);
                        Dist_to_drop = dropDistances(idx);
        
                        ReturnedEnergyObj = exp(-2*alpha*DistanceToObject)*Energy_value{pc_1}(i); %Reduce object energy based on foggy environment
                        ReturnedEnergyDrop = Energy_SoftTarget;
                    else %If no fog drops, energy from object remains the same, energy from drop is 0
                        ReturnedEnergyObj = Energy_value{pc_1}(i);
                        ReturnedEnergyDrop = 0;
                    end

                    if ReturnedEnergyDrop > Energy_Threshold && ReturnedEnergyObj < ReturnedEnergyDrop
                        %If energy of drop is highest, reposition the point
                        %to represent the drop and update the intensity
                        %value
                        xfog = cos(LiDARPitch)*cos(LiDARYaw)*Dist_to_drop;
                        yfog = cos(LiDARPitch)*sin(LiDARYaw)*Dist_to_drop;
                        zfog = sin(LiDARPitch)*Dist_to_drop;
                
                        PC_FogAdjusted(i,:) = [xfog,yfog,zfog];
                        %Convert energy to intensity
                        Intensity_value_FogAdjusted(i) = (ReturnedEnergyDrop/pulse_width)/Receiver_Area;

                    elseif ReturnedEnergyDrop < Energy_Threshold && ReturnedEnergyObj < Energy_Threshold
                        %If energy threshold is highest, remove the point
                        %by setting to zero
                        PC_FogAdjusted(i,:) = [0,0,0];
                        Intensity_value_FogAdjusted(i) = 0;
                    else
                        %If neither of the above apply, take the original
                        %point
                        PC_FogAdjusted(i,:) = PC{pc_1}(i,:);
                        %Convert energy to intensity
                        Intensity_value_FogAdjusted(i) = (ReturnedEnergyObj/pulse_width)/Receiver_Area;
                    end
                
                else
                    LiDARYaw = Angles{pc_1}(i,1); %Set LiDARYaw to the azimuth angle (LiDARYaw is used in unreal code)
                    LiDARPitch = Angles{pc_1}(i,2); %Set LiDARPitch to the azimuth angle (LiDARPitch is used in unreal code)
        
                    VBeam = 3.1416 * (BeamRadius^2)*(DistanceToObject*100); %Calulate volume of beam using pi*(r^2)*h (h is multiplied by 100 to comvert into cm) 
                    % int_sum = integral(N_D,min_FogDrop_size,max_FogDrop_size); %Integrate DSD to get mean number of drops per unit volume (per cm^2)
                    Mean_V = VBeam*int_sum; %Multiply the above to get mean number of drops in the beam
                    poisson_num = poissrnd(Mean_V); %Sample the sample of drops from poisson distribution using the mean
        
                    
                    
                    if poisson_num > 0 %If more than one drop, calculate the drop intensities
                        %This is a term in the literature for calculating
                        %power. The difference here is we are using energy
                        %rather than power so the equation fro power
                        %received becomes energy received
                        C_A_Peak_Energy = Energy_value{pc_1}(i)*((DistanceToObject^2)/beta_zero);
        
                        % Calculate dropDistances based on number of drops
                        dropDistances = ((rand(poisson_num,1) * 90) + 10) * 0.01 * DistanceToObject;

                        % Sample 100 drop distance (sort dropDistance and then divide dropDistance
                        % array into 100 subarrays and extarct mean od each
                        % sub array)
                        input = sort(dropDistances, 'ascend');
                        N = numel(input);
                        M = 100;
                        % Calculate the number of elements per chunk
                        elements_per_chunk = floor(N / M);
                        % Reshape the input array into a 2D array with the desired number of rows (chunks)
                        reshaped_input = reshape(input(1:elements_per_chunk*M), elements_per_chunk, []);
                        % Calculate the mean along the rows (dimension 1)
                        output = mean(reshaped_input, 1);
                        dropDistances = output';

                        % Define the array of linearly spaced values 't'
                        t = linspace(0, 2*pulse_width, length(dropDistances));
                                               
                        % Create a matrix 'y' using vectorization to
                        % calculate drop energy integral input equation
                        y = C_A_Peak_Energy * beta * (sin((pi.*t)./(2*pulse_width)).^2) .* ...
                            (exp(-2*alpha*(dropDistances - (speedOfLight.*t/2))) ./ ((dropDistances-(speedOfLight.*t/2)).^2));
                                                
                        % integrate y to obtain drop energies for each drop
                        dropEnergies = trapz(t,y,2);

                        % convert any Nan values to 0
                        dropEnergies(isnan(dropEnergies))=0;
                        
                        % Convert inf values to max numerical value
                        max_value = max(dropEnergies(~isinf(dropEnergies)));
                        dropEnergies(isinf(dropEnergies)) = max_value;
                        
                        % Find the index of the drop with the maximum
                        % energy and set Energy of soft target to max
                        % energy and dist_to_drop to distance of drop with
                        % max energy
                        [Energy_SoftTarget, idx] = max(dropEnergies);
                        Dist_to_drop = dropDistances(idx);
        
                        ReturnedEnergyDrop = Energy_SoftTarget;
                    else %If no fog droplets, energy from drop is 0
                        ReturnedEnergyDrop = 0;
                    end
                    
        
            
                    if ReturnedEnergyDrop > Energy_Threshold
                        %If energy of drop is higher than threshold, position the point
                        %to represent the drop and update the intensity value
                        xfog = cos(LiDARPitch)*cos(LiDARYaw)*Dist_to_drop;
                        yfog = cos(LiDARPitch)*sin(LiDARYaw)*Dist_to_drop;
                        zfog = sin(LiDARPitch)*Dist_to_drop;
                
                        PC_FogAdjusted(i,:) = [xfog,yfog,zfog];
                        %Convert energy to intensity
                        Intensity_value_FogAdjusted(i) = (ReturnedEnergyDrop/pulse_width)/Receiver_Area;

                    else %If not, the values should be zero
                        PC_FogAdjusted(i,:) = [0,0,0];
                        Intensity_value_FogAdjusted(i) = 0;
                    end
                end
            end
        end
    end

    %% Export Variables

    PC_FogSIM = PC_FogAdjusted;
    Intensity_FogSIM = Intensity_value_FogAdjusted;

end