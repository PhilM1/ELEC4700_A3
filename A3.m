%ELEC 4700 Assignment 3
close all;
clear;
clc;
warning off;
tic;

% Question 1.a) E field?
% Equation: E = V/d
yVoltage = 0.1; %Volts
xDistance = 200e-9; %Meters
E_Field = yVoltage / xDistance;
disp(['Question 1.a) Electric field is: ', num2str(E_Field), ' V/m']);

%Question 1.b) Force on Electron?
% Equation: F = qE
electron_Charge = 1.602e-19; %coulombs
electron_Force = electron_Charge * E_Field;
disp(['Question 1.b) Force on each electron: ', num2str(electron_Force), ' N']);

%Question 1.c) Acceleration on Electron from E-field?
% Equation: A = F/m , Where m is the effective mass of the electron.
mass_electron = 9.10938356e-31; %electron rest mass [kg]
mass_effective = 0.26*mass_electron; %electron effective mass [kg]

electron_accelE = electron_Force/mass_effective;
disp(['Question 1.c) Acceleration on each electron due to E-field: ', num2str(electron_accelE), ' m/s^2']);

timestep = 1e-14; %10fs
region_size_x = 200e-9; %meter
region_size_y = 100e-9; %meter
numParticles = 10000;
simLength = 1000; %number of iterations of the simulation
temperature = 300; %Kelvin

const_boltzman = 1.38064852e-23; %Boltzman constant [m^2 kg / s^2 K]
velocity_thermal = sqrt(const_boltzman*temperature/mass_effective); %thermal velocity
collisionTime = 0.2e-12; %time between particle collisions

%particle array is in following format: 
%xPos, yPos, thermal Velocity Mag, xVel, yVel, TimeSinceLastScatter
particleArray = [];

tempArray = []; %array to store temperature history
particleXPos7Array = []; %array to store 7 X positions history
particleYPos7Array = []; %array to store 7 Y positions history
scatterPathSum = 0; %store the summed path between scatters
scatterTimeSum = 0; %store the summed time between scatters
numScatters = 0; %total number of particles scattered
xCurrent_Array = []; %Array to store the currents in X direction for every iteration
electron_Concentration = 10e15*10000; %m^-2
densityMap_Array = [];
temperatureMap_Array = [];
yResolution = 20;
xResolution = yResolution * 2;

scatterProbability = 1 - exp(-timestep/collisionTime); %probability of particle being scattered

%Specs for boxes
box_Left = 80e-9;
box_Right = 120e-9;
box_Top = 60e-9;
box_Bottom = 40e-9;

%initialize the particles
for i = 1:numParticles
    %init positions
    particleArray(i,1) = rand * region_size_x;
    particleArray(i,2) = rand * region_size_y;
    vth = randn*velocity_thermal + velocity_thermal;
    
    %init velocity
    particleArray(i,3) = vth; %thermal velocity magnitude
    particleArray(i,4) = ((rand * 2) - 1)*vth; %xVelocity
    particleArray(i,5) = sqrt(particleArray(i,3)^2 - particleArray(i,4)^2); %yVelocity
    if(rand > 0.5)
        particleArray(i,5) = particleArray(i,5) * -1;
    end
    
    particleArray(i,6) = 0; %set time since last scatter to 0
        
end

%run simulation
for simCount = 1:simLength
    currTime = simCount * timestep;
        
    %scatter particles
    for i = 1:numParticles
        %update time since last scatter
        particleArray(i,6) = particleArray(i,6) + timestep;
        
        if(rand <= scatterProbability) %scatter the particle
            scatterPathSum = scatterPathSum + (particleArray(i,6)*particleArray(i,3)); %store path between scatters
            scatterTimeSum = scatterTimeSum + particleArray(i,6); %store time between scatters
            particleArray(i,6) = 0; %reset time since last scatter            
            numScatters = numScatters + 1;
            particleArray(i,3) = randn*velocity_thermal + velocity_thermal; %randomize velocity
            particleArray(i,4) = ((rand * 2) - 1)*particleArray(i,3); %xVelocity
            particleArray(i,5) = sqrt(particleArray(i,3)^2 - particleArray(i,4)^2); %yVelocity
            if(rand > 0.5)
                particleArray(i,5) = particleArray(i,5) * -1;
            end
        end
    end
    
    %update x velocity due to E-field
    particleArray(:,4) = particleArray(:,4) + timestep * electron_accelE;
    
    %update particle positions
    new_xPos = particleArray(:,1) + timestep * particleArray(:,4);
    new_yPos = particleArray(:,2) + timestep * particleArray(:,5);    
    
    
    %check boundary conditions
    for i = 1:numParticles
        if(new_xPos(i) < 0) %pass through (x-dir)
            new_xPos(i) = new_xPos(i)+region_size_x;
        elseif(new_xPos(i) > region_size_x) %bounce off boundary
            new_xPos(i) = new_xPos(i) - region_size_x;
        end
        
        if(new_yPos(i) < 0) %bounce off boundary (y-dir)
            new_yPos(i) = abs(new_yPos(i));
            particleArray(i,5) = particleArray(i,5) * -1; %swap direction
        elseif(new_yPos(i) > region_size_y) %bounce off boundary
            new_yPos(i) = 2*region_size_y - new_yPos(i);
            particleArray(i,5) = particleArray(i,5) * -1; %swap direction
        end
       
    end
    
    %set new positions
    particleArray(:,1) = new_xPos;
    particleArray(:,2) = new_yPos;
    
    %Scatter plotting of the particles uncomment if you wanna see live plot
%     figure(4);
%     scatter(particleArray(:,1), particleArray(:,2), 5);
%     axis([0, region_size_x, 0, region_size_y]);
%     title(['Scatter Plot of Particles in a box. Iteration: ', num2str(simCount)]);
%     xlabel('X (m)');
%     ylabel('Y (m)');


    %save the particle position history to be plotted later
    particleXPos7Array(1:7,simCount) = particleArray(1:7,1);
    particleYPos7Array(1:7,simCount) = particleArray(1:7,2);

    %save the temperature history to be plotted later
    %tempArray(simCount) = mean(particleArray(:,3))^2*mass_effective/const_boltzman;

    %store the X current
    xCurrent_Array(simCount) = electron_Charge * electron_Concentration * mean(particleArray(:,4)) * region_size_x;

    %pause(0.0001); %quick pause so that the graphs can be displayed live
    if(mod(simCount,100) == 0) %Update the console with current simCount
        disp(['simcount: ', num2str(simCount), ' of ', num2str(simLength)]);
    end
    
    if(simCount == simLength) %last iteration, populate Maps
        
        xDiv = region_size_x / xResolution;
        yDiv = region_size_y / yResolution;
        
        for yCount = 1:yResolution
            for xCount = 1:xResolution
                densityMap_Array(xCount,yCount) = 0;
                temperatureMap_Array(xCount,yCount) = 0;
                
                for eCount = 1:numParticles
                    particleArray(eCount,1);
                    if(particleArray(eCount,1) >= (xCount*xDiv-xDiv) && particleArray(eCount,1) < (xCount*xDiv))
                        if(particleArray(eCount,2) >= (yCount*yDiv-yDiv) && particleArray(eCount,2) < (yCount*yDiv))
                            densityMap_Array(xCount,yCount) = densityMap_Array(xCount,yCount) + 1;
                            temptemp = sqrt(particleArray(eCount,4).^2+particleArray(eCount,5).^2)*mass_effective/(3*const_boltzman);
                            temperatureMap_Array(xCount,yCount) = temperatureMap_Array(xCount,yCount) + temptemp;
                        end                        
                    end
                end                
            end
        end   
        temperatureMap_Array = temperatureMap_Array ./ densityMap_Array;
    end
end

disp('----------');
disp(['Simulation ended successfully after ', num2str(toc), ' seconds.']);
disp('----------');
disp('Simulation Specs:');
disp(['Timestep Size: ', num2str(timestep), ' seconds']);
disp(['Number of Simulation Iterations: ', num2str(simLength)]);
disp(['Number of Particles: ', num2str(numParticles)]);
disp('----------');
disp('Simulation Results:');
disp(['Mean Free Path (MFP): ', num2str(scatterPathSum/numScatters), ' meters']);
disp(['Time between Collisions: ', num2str(scatterTimeSum/numScatters), ' seconds']);

%plot the trajectory of the particles
figure(10);
hold on;
for i = 1:7
    plot(particleXPos7Array(i,:),particleYPos7Array(i,:));
end
xlim([0,region_size_x]);
ylim([0,region_size_y]);
title(['Particle Trajectory, Simulation Count: ', num2str(simCount), ', Timestep: ', num2str(timestep)]);
xlabel('X (m)');
ylabel('Y (m)');

%Plot the X current over time
figure(11);
x = linspace(timestep,simLength*timestep,simLength);
plot(x,xCurrent_Array);
title('X Current over Time');
xlabel('Time (s)');
ylabel('Current (A)');

%plot the density map
figure(12);
[xMesh, yMesh] = meshgrid(linspace(0,region_size_x,xResolution),linspace(0,region_size_y,yResolution));
surf(xMesh,yMesh,transpose(densityMap_Array));
title('Density Map');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Count');

%plot the temperature map
figure(13);
[xMesh, yMesh] = meshgrid(linspace(0,region_size_x,xResolution),linspace(0,region_size_y,yResolution));
surf(xMesh,yMesh,transpose(temperatureMap_Array));
title('Temperature Map');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Temperature (K)');
