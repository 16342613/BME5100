% Author: Irish Senthilkumar
% ID: 16342613
% Course:
% Module: BME5100 Advanced Computational Biomechanics
% Assignment: Summary of experimental and computational investigation
%             anisotropic hyperelastic behaviour of arterial tissue. (3)

clc;
clear;

K = 1;
G = 1;
k1 = 1;
k2 = 1;
j = 1;
fibreAngle = 20; % For second test change gamma to 90 minus 30

simulatedAxialStressStrain = computeStressStrainCurve(1.52, 30, K, G, k1, k2, j, 90 - fibreAngle, 2);
simulatedCircumferentialStressStrain = computeStressStrainCurve(1.52, 30, K, G, k1, k2, j, fibreAngle, 1);

recordedAxialStressStrain = xlsread("Aniso_Aorta_2017.xlsx", "Axial");
recordedCircumferentialStressStrain = xlsread("Aniso_Aorta_2017.xlsx", "Circum");

experimentalAxialStrain = [];
experimentalAxialStress = [];
experimentalCircumferentialStrain = [];
experimentalCircumferentialStress = [];

% Now get the nominal stress and strain from the data
originalLength = 0.021;
originalWidth = 0.021; % WAS 0.01
thickness = 0.003;
originalAxialArea = originalLength * thickness;
originalCircumferentialArea = originalWidth * thickness;

for i = 1:length(recordedAxialStressStrain)
    lengthChange = recordedAxialStressStrain(i,1)/1000;

    experimentalAxialStrain = [experimentalAxialStrain lengthChange/originalLength];
    experimentalAxialStress = [experimentalAxialStress computeNominalStress(recordedAxialStressStrain(i,2), originalAxialArea, originalLength, lengthChange)];
end

for i = 1:length(recordedCircumferentialStressStrain)
    widthChange = recordedCircumferentialStressStrain(i,1)/1000;

    experimentalCircumferentialStrain = [experimentalCircumferentialStrain widthChange/originalWidth];
    experimentalCircumferentialStress = [experimentalCircumferentialStress computeNominalStress(recordedCircumferentialStressStrain(i,2), originalCircumferentialArea, originalWidth, widthChange)];
end

% Plot the results
plot(experimentalAxialStrain, experimentalAxialStress);
hold on
plot(simulatedAxialStressStrain(:,1) - 1, simulatedAxialStressStrain(:,2)*1000);
hold on
plot(experimentalCircumferentialStrain, experimentalCircumferentialStress);
hold on
plot(simulatedCircumferentialStressStrain(:,1) - 1, simulatedCircumferentialStressStrain(:,2)*1000);
hold off
legend('Axial - Experimental','Axial - Simulated', 'Circumferential - Experimental','Circumferential - Simulated');
%legend('Axial - Experimental','Circumferential - Experimental');
title('Axial and Circumferential Stress-Strain Curves');
xlabel('Strain (Proportional)');
ylabel('Stress (Pa)');

disp("done");

function sigma_nominal = computeNominalStress(currentForce, originalArea, originalLength, lengthChange)
    sigma_nominal = (currentForce / originalArea) * (1 + lengthChange / originalLength);
    % sigma_nominal = currentForce / originalArea;
end

function stressStrain = computeStressStrainCurve(maximumStretch, measurementCount, K, G, k1, k2, j, fibreAngle, axis)
    x = linspace(1, maximumStretch, measurementCount + 1);
    stressStrain = [1 0];

    for i = 1:length(x)
        deformationGradient = [x(i) 0 0; ...
                              0 1/sqrt(x(i)) 0;
                              0 0 1/sqrt(x(i))];
    
        % Compute C, the right cauchy green tensor
        cauchyGreen_right = (deformationGradient.') * deformationGradient;
        % Compute B, the left cauchy green tensor
        cauchyGreen_left = deformationGradient * (deformationGradient.');
    
        % Compute the unit vector indicating the direction of fibre 4
        m_4 = [cosd(fibreAngle); sind(fibreAngle); 0];
        % Compute the unit vector indicating the direction of fibre 6
        m_6 = [cosd(-fibreAngle); sind(-fibreAngle); 0];
    
        % Compute the 4th anisotropic invariant
        I_4 = (m_4.') * cauchyGreen_right * m_4;
        % Compute the 6th anisotropic invariant
        I_6 = (m_6.') * cauchyGreen_right * m_6;
    
        % Compute the volumetric stress of the matrix
        sigma_vol = (2/(2/K)) * (j - 1) * (eye(3, 3));
    
        % Compute the isochoric stress of the matrix
        I_1 = trace(cauchyGreen_left);
    
        sigma_iso = (2/j) * (G/2) * ((cauchyGreen_left / (j^(2/3))) - (I_1/(3*j^(2/3))) * eye(3, 3));
    
        % Compute the anisotropic stress of the fibre 4
        sigma_aniso_f4 = 2 * k1 * (I_4 - 1) * exp(k2 * (I_4 - 1)^2) * ((deformationGradient * m_4) * (deformationGradient * m_4).');
        % Compute the anisotropic stress of the fibre 6
        sigma_aniso_f6 = 2 * k1 * (I_6 - 1) * exp(k2 * (I_6 - 1)^2) * ((deformationGradient * m_6) * (deformationGradient * m_6).');
    
        % Sum up all of the stresses together to get the overall stress
        sigma_total = sigma_vol + sigma_iso + sigma_aniso_f4 + sigma_aniso_f6;
    
        % Compute the first Piola-Kirchoff stress tensor
        P = j * sigma_total * (deformationGradient.');
    
        stressStrain = [stressStrain ; [x(i) P(axis, axis)]];
    end
end