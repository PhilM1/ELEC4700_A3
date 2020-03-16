%ELEC 4700 Assignment 3
close all;
clear;
clc;
warning off;
tic;

% =====================================================
%                          #1
% =====================================================
disp('Question 1:');
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
pause(0.1);


% =====================================================
%                          #2 & 3
% =====================================================
disp('===========================');
disp('Question 2:');
ySize = 100; %size of the matrix
xSize = ySize*2; %size of the matrix
xDel = region_size_x/xSize; %delta X
yDel = region_size_y/ySize; %delta Y
maxIterations = 10000; %maximum number of iterations

vMatrix = zeros(ySize,xSize); %create our matrix of values
oldvMatrix = zeros(ySize,xSize); %another matrix to be used to clone the new one
vecMatrix = zeros(ySize,xSize); %matrix of the vectors
[xMesh, yMesh] = meshgrid(linspace(0,region_size_x,xSize),linspace(0,region_size_y,ySize));

%set boundary values
TBC_Left = 0.1;
TBC_Right = 0;
TBC_Top = 0; %this will be unbounded
TBC_Bottom = 0; %this will be unbounded

barrier_Height = 30;
barrier_Width = 40;

for i = 1:maxIterations %keep running until reached max iterations
    for m = 1:xSize %loop through columns
        for n = 1:ySize %loop through rows
            
            if(m == 1) %left boundary
                vMatrix(n,m) = TBC_Left;
            elseif(m == xSize) %right boundary
                vMatrix(n,m) = TBC_Right;
            elseif(n == 1) %top boundary
                %vMatrix(n,m) = TBC_Top;
                vMatrix(n,m) = (oldvMatrix(n,m-1) + oldvMatrix(n,m+1) + oldvMatrix(n+1,m))/3;
            elseif(n == ySize) %bottom boundary
                %vMatrix(n,m) = TBC_Bottom;
                vMatrix(n,m) = (oldvMatrix(n,m-1) + oldvMatrix(n,m+1) + oldvMatrix(n-1,m))/3;
            else
                vMatrix(n,m) = (oldvMatrix(n,m-1) + oldvMatrix(n,m+1) + oldvMatrix(n-1,m) + oldvMatrix(n+1,m))/4;
            end   
            
            if(m>=xSize/2-barrier_Width/2 && m<=xSize/2+barrier_Width/2 && (n<=barrier_Height || n>=ySize-barrier_Height)) %inside one of the boxes
                vMatrix(n,m) = vMatrix(n,m)*10^-2;
            end
        end
    end        
    [xVectors, yVectors] = gradient(oldvMatrix(:,:));   
    if(mod(i,1000) == 0) %Update the console with current simCount
        disp(['simcount: ', num2str(i), ' of ', num2str(maxIterations)]);
    end
    oldvMatrix = vMatrix;
end

figure(21)
surf(xMesh, yMesh, vMatrix);
shading interp
title('V(x,y) for Bottle-neck (0.1V left side, 0V right side)');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Voltage (V)');

figure(22)
quiver(xMesh, yMesh, -xVectors, -yVectors);
title('Electric Field for Bottle-neck (0.1V left side, 0V right side)');
xlabel('X (m)');
ylabel('Y (m)');
pause(0.1);
disp('FD Bottle-neck Simualation Completed');

%Now let's use the Efield to modify our velocities instead of same
%acceleration
disp('Starting MC simulation using new E-field (This might take a while)');
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


%initialize the particles
numParticles = 1000;
for i = 1:numParticles
    %init positions
    particleArray(i,1) = rand * region_size_x;
    particleArray(i,2) = rand * region_size_y;
    vth = randn*velocity_thermal + velocity_thermal;
    
    while(particleArray(i,1)>=(region_size_x/2-barrier_Width/2*1e-9) && particleArray(i,1)<=(region_size_x/2+barrier_Width/2*1e-9) && (particleArray(i,2)<=(barrier_Height*1e-9) || particleArray(i,2)>=(region_size_y-barrier_Height*1e-9)))
        %pick new location, you're in the boxes!
        particleArray(i,1) = rand * region_size_x;
        particleArray(i,2) = rand * region_size_y;
    end
    
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
simLength = 1000;
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
    for eCount = 1:numParticles
        for yCount = 1:ySize
            for xCount = 1:xSize            
                if(particleArray(eCount,1) >= (xCount*xDel-xDel) && particleArray(eCount,1) < (xCount*xDel))
                    if(particleArray(eCount,2) >= (yCount*yDel-yDel) && particleArray(eCount,2) < (yCount*yDel))
                        particleArray(eCount,4) = particleArray(eCount,4) + timestep * -xVectors(yCount,xCount)*electron_Charge/mass_effective/xDel;
                        particleArray(eCount,5) = particleArray(eCount,5) + timestep * -yVectors(yCount,xCount)*electron_Charge/mass_effective/xDel; 
                    end                        
                end
            end                
        end
    end  
    
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
       
        %You're in the boxes
        box_Left = region_size_x/2-barrier_Width/2*1e-9;
        box_Right = region_size_x/2+barrier_Width/2*1e-9;
        box_Bottom = barrier_Height*1e-9;
        box_Top = region_size_y-barrier_Height*1e-9;
        if(new_xPos(i) >= box_Left && new_xPos(i) <= box_Right && (new_yPos(i) <= box_Bottom || new_yPos(i) >= box_Top))
            if(particleArray(i,1) < box_Left)
                tempx = box_Left - abs(new_xPos(i) - particleArray(i,1));
                if(~(tempx >= box_Left && tempx <= box_Right && (new_yPos(i) <= box_Bottom || new_yPos(i) >= box_Top))) %not in box
                    new_xPos(i) = tempx;
                    particleArray(i,4) = particleArray(i,4) * -1; %swap direction
                end
            elseif(particleArray(i,1) > box_Right)
                tempx = box_Right + abs(new_xPos(i) - particleArray(i,1));
                if(~(tempx >= box_Left && tempx <= box_Right && (new_yPos(i) <= box_Bottom || new_yPos(i) >= box_Top))) %not in box
                    new_xPos(i) = tempx;
                    particleArray(i,4) = particleArray(i,4) * -1; %swap direction
                end
            elseif(particleArray(i,2) < box_Top)
                tempy = box_Top - abs(new_yPos(i) - particleArray(i,2));
                if(~(new_xPos(i) >= box_Left && new_xPos(i) <= box_Right && (tempy <= box_Bottom || tempy >= box_Top))) %not in box
                    new_yPos(i) = tempy;
                    particleArray(i,5) = particleArray(i,5) * -1; %swap direction
                end
            elseif(particleArray(i,2) > box_Bottom)
                tempy = box_Bottom + abs(new_yPos(i) - particleArray(i,2));
                if(~(new_xPos(i) >= box_Left && new_xPos(i) <= box_Right && (tempy <= box_Bottom || tempy >= box_Top))) %not in box
                    new_yPos(i) = tempy;
                    particleArray(i,5) = particleArray(i,5) * -1; %swap direction
                end
            end
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
        temperatureMap_Array(isnan(temperatureMap_Array)) = 0;
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
figure(23);
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
figure(31);
x = linspace(timestep,simLength*timestep,simLength);
plot(x,xCurrent_Array);
title('X Current over Time');
xlabel('Time (s)');
ylabel('Current (A)');

%plot the density map
figure(32);
[xMesh, yMesh] = meshgrid(linspace(0,region_size_x,xResolution),linspace(0,region_size_y,yResolution));
surf(xMesh,yMesh,transpose(densityMap_Array));
title('Density Map');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Count');

%plot the temperature map
figure(33);
[xMesh, yMesh] = meshgrid(linspace(0,region_size_x,xResolution),linspace(0,region_size_y,yResolution));
surf(xMesh,yMesh,transpose(temperatureMap_Array));
title('Temperature Map');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Temperature (K)');

avgCurrent = mean(xCurrent_Array(100:end));
BottleneckWidth = 100-barrier_Height*2;
disp(['Avg Curr: ', num2str(avgCurrent)]);
disp(['B-n width: ', num2str(BottleneckWidth)]);

%plot the avg current vs bottleneck width
avgCurArr = [0.053111, 0.069134, 0.088485, 0.096499, 0.10608, 0.11388, 0.1217, 0.12317, 0.12871];
bneckWArr = [10, 20, 30, 40, 50, 60, 70, 80, 90];
figure(34);
plot(bneckWArr,avgCurArr);
title('Average Current vs Bottle-neck Width');
xlabel('Bottle-neck Width (nm)');
ylabel('Current (A)');

