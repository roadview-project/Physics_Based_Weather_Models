% Time taken to apply rain (seconds) for 20 frames: 480.3364
% Frame Rate (per second) = 0.0417 (1 frame every 24 seconds)

% Inputs and Output Definition (s): 4.51e-05
% Setting Weather Paramaters (s): 6.09e-05
% Setting LiDAR Paramaters (s): 0.0001169
% Obtaining Energy Threshold (s): 0.0007903
% Obtaining Rain DSD (s): 0.20079
% Obtaining Extinction Coefficient Alpha (s): 0.087218
% Calculating Emission Angles (s): 0.20997
% Setting Zero-Range Emission Angles (s): 0.23951
% Appling Rain Model (s): 23.1518
% Saving Simuated PCs (s): 6.86e-05

%% Rain Model Function (Takes cell containing single PC coordinates and cell containing single PC intensity values as input + rain rate)
function [PC_RainSIM, Intensity_RainSIM] = LiDAR_RainModel_Ouster(Clear_PC, Clear_Intensity, Rain_Rate)
    %% Define Number of Input pointclouds and Output pointclouds

    pointCloudNum = 1; % Number of clear pointclouds inputted
    Repetitions = 1; % The number of point clouds with added rain noise that will be generated

    PC = Clear_PC;
    Intensity_value = Clear_Intensity;

    %% Set Parameters of the weather

    Rain_Intensity = Rain_Rate; %mm/h (for THI rain rates are 10, 25 and 50 mmh)
    min_DropSize = 0.5; %mm (from Met Office Factsheeet - https://www.metoffice.gov.uk/binaries/content/assets/metofficegovuk/pdf/research/library-and-archive/library/publications/factsheets/factsheet_3-water-in-the-atmosphere_2023.pdf)
    max_DropSize = 6; %mm
    D_c = 1; % droplet diameter of max probability (mm) - i.e. most probable rain drop diameter (see section 4.5 in: https://www.mdpi.com/1424-8220/23/15/6891)
    Water_RefractIndex = complex(1.328,4.86*(10^-7)); % This needs to be the complex refractive index of the drop

    %% Set Parameters of the LiDAR

    % Datasheet link - https://data.ouster.io/downloads/datasheets/datasheet-revd-v2p0-os1.pdf
    BeamRadius = 0.0095/2; %m
    LiDAR_Max_Range = 120; %m
    wavelength = 865*(10^-9); %m (The wavelength of Ouster LiDAR is 865nm)
    pulse_width = 10*(10^-15); %s (Estimated by Valentina and Jonathan)
    Receiver_Area = 1*(10^-6); %m^2 (Estimated by Valentina and Jonathan - original estimate was 1mm^2 = 1*10^-6m^2)
    Channels = 128; %number
    Horizontal_FOV = 360*(pi/180); %radians (360 degrees)
    min_az_angle = -Horizontal_FOV/2; %radians
    max_az_angle = Horizontal_FOV/2; %radians
    Vertical_FOV = 45*(pi/180); %radians (45 degrees)
    min_z_angle = -22.5*(pi/180); %radians
    max_z_angle = 22.5*(pi/180); %radians
    Horizontal_Res = Horizontal_FOV/2048; %radians (angular resolution)
    Vertical_Res = Vertical_FOV/Channels; %radians (angular resolution)
    
    Planck = 6.62607015*(10^-34); %js
    speedOfLight = 299792458; %m/s (Actual Value = 299792458 m/s, approximation = 3*10^8 m/s)
    
    Energy_value_Max = (2^16)*Planck*(speedOfLight/wavelength); %16 bit unsigned int - Ouster - Max_energy_value = max intensity value*planck*speed_of_light

    %% Find Energy of LiDAR Returns

    %Ouster LiDAR declares intensity as number of photons up to the maximum of
    %2^16 photons. Hence we simply mutliply the intensity value by the energy
    %of a photon to get the received energy
    
    Energy_value = cell(pointCloudNum,1);
    
    for e=1:pointCloudNum
        Energy_value{e} = Intensity_value{e}*Planck*(speedOfLight/wavelength);
    end
    
    Energy_Threshold = 0.00005*Energy_value_Max; % (currently set to 0.005%) - The reason why the threshold is very low is because most of the intensity values in the PCs are between 0 - 15, (almost all of them roughly 240000/260000)

    %% Distribution function

    %Define the dropsize distribution for the rain as Marshall Palmer. Create a
    %distribution object for use in sampling the size of the drops.
    
    % Marshall-Palmer DSD
    N_D = @(D) 8000*exp(-4.1*(Rain_Intensity^-0.21)*D); % Equation found in Rasshofer et al
    DSD = makedist('Exponential','mu',4.1*(Rain_Intensity^-0.21)); % This is the same as N(D) but contains numerical values which can be sampled
   
    %% Extinction coefficient

    %Calulation of the extinction coefficient alpha. The equation must be
    %changed to have the correct DSD. This is in dB/km so must be converted to
    %dB/m. x_EXT is the size parameter. 

    x_EXT = (pi/wavelength)*(D_c*(10^-3)); %This is approximate, should really use D
    m_EXT = real(Water_RefractIndex);
   
    %%% Marshall Palmer
    alpha = (pi/8)*integral(@(D) (D.^2).*(2./(x_EXT.^2)).*C_EXT(x_EXT,m_EXT,round(x_EXT + (4.*(x_EXT.^(1/3))) + 10,0)).*(8000*exp(-4.1*(Rain_Intensity^-0.21)*D)),min_DropSize,max_DropSize);

    % Convert to dB/m
    alpha = alpha/1000; %dB/m
    
    %% Backscattering coefficient (NOT USED IN RAIN MODEL)
    
    %Calulation of the extinction coefficient beta. The equation must be
    %changed to have the correct DSD. This is in dB/km so must be converted to
    %dB/m
    
    %Gamma
    % beta = (pi/8)*integral(@(D) (D.^2).*((1./(((pi.*D)./wavelength).^2))*C_BACK(x_EXT,m_EXT,round(((pi.*D)./wavelength) + (4.*(((pi.*D)./wavelength).^(1/3))) + 10,0))).*((gamma*rho*(b^((alpha_rain+1)/gamma)))/gamma_Caps).*((D/2).^alpha_rain).*exp(-b.*((D/2).^gamma)),min_DropSize,max_DropSize);
    
    %Marshall-Palmer
    % beta = (pi/8)*integral(@(D) (D.^2).*((1./(((pi.*D)./wavelength).^2))*C_BACK(x_EXT,m_EXT,round(((pi.*D)./wavelength) + (4.*(((pi.*D)./wavelength).^(1/3))) + 10,0))).*(8000*exp(-4.1*(Rain_Intensity^-0.21)*D)),min_DropSize,max_DropSize);
    
    %Measured Value
    % beta = (pi/8)*integral(@(D) (D.^2).*((1./(((pi.*D)./wavelength).^2))*C_BACK(x_EXT,m_EXT,round(((pi.*D)./wavelength) + (4.*(((pi.*D)./wavelength).^(1/3))) + 10,0))).*(N_0*exp(-lambda_Cap*D)),min_DropSize,max_DropSize);
    % 
    % beta = beta/1000; %dB/m
    
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

    %As we know we have points with angle set to 1000 we can use this to find
    %the points with zero range. These points have no angle we can extract so
    %we need to predict an emission angle. This is such that we can still place
    %a raindrop point even if no point was present originally.
    
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
    
    %% Loop through points and change points for rain

    %These variables are to store the points in the rain model point cloud
    PC_RainAdjusted = cell(pointCloudNum,Repetitions);
    Intensity_value_RainAdjusted = cell(pointCloudNum,Repetitions);
 
    
    for pc_1 = 1:pointCloudNum %Loop through the point clouds
        for rep=1:Repetitions %Loop through the repetitions
            
            %Initialise the PC and intensity values for this point cloud and
            %repetition
            PC_RainAdjusted = zeros(length(PC{pc_1}(:,1)),length(PC{pc_1}(1,:)));
            Intensity_value_RainAdjusted = zeros(length(PC{pc_1}(:,1)),1);
                      
            
            for i=1:length(PC{pc_1}(:,1)) %Loop through the points in the original clean PC
                DistanceToObject = Range{pc_1}(i,1); %Set the DistanceTo Object variable to the Range for that point (DistanceToObject is used in the unreal code)
             
                %If Range is greater than zero, we adjust intensity of
                %object and intensity of drop. If not, we just calculate
                %intensity of drop.
                if DistanceToObject > 0 
        
                    LiDARYaw = Angles{pc_1}(i,1); %Set LiDARYaw to the azimuth angle (LiDARYaw is used in unreal code)
                    LiDARPitch = Angles{pc_1}(i,2); %Set LiDARPitch to the elevation angle (LiDARPitch is used in unreal code)
        
                    VBeam = 3.1416 * (BeamRadius^2)*DistanceToObject; %Calulate volume of beam using pi*(r^2)*h
                    int_sum = integral(N_D,min_DropSize,max_DropSize); %Integrate DSD to get mean number of drops per unit volume
                    Mean_V = VBeam*int_sum; %Multiply the above to get mean number of drops in the beam
                    poisson_num = poissrnd(Mean_V); %Sample the sample of drops from poisson distribution using the mean
        
                    if poisson_num > 0 %If more than one drop, calculate the drop intensities
                        dropDistances = zeros(poisson_num,1);
                        dropSizes = zeros(poisson_num,1);
                        dropReflectivities = zeros(poisson_num,1);
                        dropReturnedEnergies = zeros(poisson_num,1);
                        dropMaxEnergy = 0;
                        Dist_to_drop = 0;
                        for d2d = 1:poisson_num %Loop through the drops
                            dropDistances(d2d) = ((rand*90)+10)*0.01*DistanceToObject; %Uniform distribution between 10% and 100% of distance to object
                            dropSizes(d2d) = random(DSD); %Sample drop size from DSD
                            dropReflectivities(d2d) = (real((Water_RefractIndex-1)/(Water_RefractIndex+1)))^2 + (imag((Water_RefractIndex-1)/(Water_RefractIndex+1)))^2; % Drop reflectivity calculated using Fresnel
                            %Calculate the returned energy from the drop (pulse width * returned power) = E = P*t
                            dropReturnedEnergies(d2d) = pulse_width*((dropReflectivities(d2d)*exp(-2*alpha*dropDistances(d2d)))/(dropDistances(d2d)^2))*(min([(dropSizes(d2d)/(BeamRadius*1000*2))^2,1])); % No beam divergence
                            %If the drop has the highest energy, set its
                            %distance and energy to the stored values
                            % disp(dropReturnedEnergies(d2d))
                            if dropReturnedEnergies(d2d) > dropMaxEnergy
                                dropMaxEnergy = dropReturnedEnergies(d2d);
                                Dist_to_drop = dropDistances(d2d);
                            elseif dropReturnedEnergies(d2d) == dropMaxEnergy
                                %If the energy is exactly the same, take
                                %the closest drop
                                if dropDistances(d2d) < Dist_to_drop
                                    Dist_to_drop = dropDistances(d2d);
                                end
                            end
                        end
                        ReturnedEnergyObj = Energy_value{pc_1}(i)*exp(-2*alpha*DistanceToObject); %Reduce object energy based on rainy environment
                        ReturnedEnergyDrop = dropMaxEnergy;
                    else %If no rain drops, energy from object remains the same, energy from drop is 0
                        ReturnedEnergyObj = Energy_value{pc_1}(i);
                        ReturnedEnergyDrop = 0;
                    end

                
                    if ReturnedEnergyDrop > Energy_Threshold && ReturnedEnergyObj < ReturnedEnergyDrop
                        %If energy of drop is highest, reposition the point
                        %to represent the drop and update the intensity
                        %value
                        xrain = cos(LiDARPitch)*cos(LiDARYaw)*Dist_to_drop;
                        yrain = cos(LiDARPitch)*sin(LiDARYaw)*Dist_to_drop;
                        zrain = sin(LiDARPitch)*Dist_to_drop;
                
                        PC_RainAdjusted(i,:) = [xrain,yrain,zrain];
                        %Convert energy to intensity
                        Intensity_value_RainAdjusted(i) = (ReturnedEnergyDrop/pulse_width)/Receiver_Area;
                    elseif ReturnedEnergyDrop < Energy_Threshold && ReturnedEnergyObj < Energy_Threshold
                        %If energy threshold is highest, remove the point
                        %by setting to zero
                        PC_RainAdjusted(i,:) = [0,0,0];
                        Intensity_value_RainAdjusted(i) = 0;
                    else
                        %If neither of the above apply, take the original
                        %point
                        PC_RainAdjusted(i,:) = PC{pc_1}(i,:);
                        %Convert energy to intensity
                        Intensity_value_RainAdjusted(i) = (ReturnedEnergyObj/pulse_width)/Receiver_Area;
                    end

                else
                    
                    LiDARYaw = Angles{pc_1}(i,1); %Set LiDARYaw to the azimuth angle (LiDARYaw is used in unreal code)
                    LiDARPitch = Angles{pc_1}(i,2); %Set LiDARPitch to the azimuth angle (LiDARPitch is used in unreal code)
        
        
                    VBeam = 3.1416 * (BeamRadius^2)*DistanceToObject; %Calulate volume of beam using pi*(r^2)*h
                    int_sum = integral(N_D,min_DropSize,max_DropSize); %Integrate DSD to get mean number of drops per unit volume
                    Mean_V = VBeam*int_sum; %Multiply the above to get mean number of drops in the beam
                    poisson_num = poissrnd(Mean_V); %Sample the sample of drops from poisson distribution using the mean
        
                    if poisson_num > 0 %If more than one drop, calculate the drop intensities
                        dropDistances = zeros(poisson_num,1);
                        dropSizes = zeros(poisson_num,1);
                        dropReflectivities = zeros(poisson_num,1);
                        dropReturnedEnergies = zeros(poisson_num,1);
                        dropMaxEnergy = 0;
                        Dist_to_drop = 0;
                        for d2d = 1:poisson_num %Loop through the drops
                            dropDistances(d2d) = ((rand*90)+10)*0.01*LiDAR_Max_Range; %Uniform distribution between 10% and 100% of distance to object
                            dropSizes(d2d) = random(DSD); %Sample drop size from DSD
                            dropReflectivities(d2d) = (real((Water_RefractIndex-1)/(Water_RefractIndex+1)))^2 + (imag((Water_RefractIndex-1)/(Water_RefractIndex+1)))^2; % Drop reflectivity calculated using Fresnel
                            %Calculate the returned energy from the drop
                            dropReturnedEnergies(d2d) = pulse_width*((dropReflectivities(d2d)*exp(-2*alpha*dropDistances(d2d)))/(dropDistances(d2d)^2))*(min([(dropSizes(d2d)/(BeamRadius*1000*2))^2,1])); % No beam divergence
                            %If the drop has the highest energy, set its
                            %distance and energy to the stored values
                            % disp(dropReturnedEnergies(d2d))
                            if dropReturnedEnergies(d2d) > dropMaxEnergy
                                dropMaxEnergy = dropReturnedEnergies(d2d);
                                Dist_to_drop = dropDistances(d2d);
                            elseif dropReturnedEnergies(d2d) == dropMaxEnergy
                                %If the energy is exactly the same, take
                                %the closest drop
                                if dropDistances(d2d) < Dist_to_drop
                                    Dist_to_drop = dropDistances(d2d);
                                end
                            end
                        end
                        ReturnedEnergyDrop = dropMaxEnergy;
                    else %If no rain drops, energy from drop is 0
                        ReturnedEnergyDrop = 0;
                    end
        
                    
                
                    if ReturnedEnergyDrop > Energy_Threshold
                        %If energy of drop is higher than threshold, position the point
                        %to represent the drop and update the intensity value
                        xrain = cos(LiDARPitch)*cos(LiDARYaw)*Dist_to_drop;
                        yrain = cos(LiDARPitch)*sin(LiDARYaw)*Dist_to_drop;
                        zrain = sin(LiDARPitch)*Dist_to_drop;
                
                        PC_RainAdjusted(i,:) = [xrain,yrain,zrain];
                        %Convert energy to intensity
                        Intensity_value_RainAdjusted(i) = (ReturnedEnergyDrop/pulse_width)/Receiver_Area;

                    else %If not, the values should be zero
                        PC_RainAdjusted(i,:) = [0,0,0];
                        Intensity_value_RainAdjusted(i) = 0;
                    end
                    
                end
            end
        end
    end
    
    %% Export Variables (Save PCs (Sim_Rain (coordinates and intensity values))

    PC_RainSIM = PC_RainAdjusted;
    Intensity_RainSIM = Intensity_value_RainAdjusted;

end