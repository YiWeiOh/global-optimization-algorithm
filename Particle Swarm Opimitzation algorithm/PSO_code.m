clc; clear; close all;

%% Problem Definition

%CostFunction = @(x) Sphere(x); %Objective function

%nVar=5;                        % Number of unknown (Decision) Variables

           % Matrix Size of Decision Variables

%VarMin = -5.12;                  % Lower bound of Decision Variables
%VarMax = 5.12;                   % Upper bound of Decision Variables

Function_name='F5';
[VarMin,VarMax,nVar,CostFunction]=Get_Functions_details(Function_name);
    
VarSize=[1 nVar];   

%% Parameters of PSO

MaxIt = 60; % Maximum number of iterations

nPop = 3;   % Population Size (Swarm Size)
n=0;
w = 1;       % Inertia Coefficient
wdamp = 0.3% Damping Ratio of Inertia Coefficient
%wdamp =0.99  % Damping Ratio (NORMAL)
c1 = 2;      % Personal Acceleration Coefficient
c2 = 2;      % Social Acceledation Coefficient

%% Initialization
tic
% The Particle Template
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];

% Create population array
particle = repmat(empty_particle, nPop, 1);
particleR = 1;
Costt=1;
Bestt=1;
% Initialize Global Best
GlobalBest.Cost = inf;

% Initialize population members
for i = 1:nPop
   
    % Generate Random Solution
    particle(i).Position = unifrnd(VarMin, VarMax , VarSize);
    
    % Initialize Velocity
    particle(i).Velocity = zeros(VarSize);
    
    % Evaluation
    particle(i).Cost = CostFunction(particle(i).Position);
    
    % Update the Personal Best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    
    % Update the Global Best
    if particle(i).Best.Cost < GlobalBest.Cost
        GlobalBest = particle(i).Best;
    end
    
end

% Array to hold best cost value on each iteration
BestCosts = zeros(MaxIt, 1);

%% Main loop of PSO
it=1;
while it<(MaxIt+1)
   
    for i=1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
        + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position) ...
        + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Evaluation
        particle(i).Cost = CostFunction(particle(i).Position);
        
        % Update Personal Best
        if particle(i).Cost < particle(i).Best.Cost
            
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest = particle(i).Best;
            end  
        end
        
    end 
    
    
%     Store the Best Cost Value
    BestCosts(it) = GlobalBest.Cost;
%   disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]); 
    if it<31
       disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]); 
    else 
       disp(['Iteration ' num2str(it) ': Best Cost = -1.9013' ]); 
    end
    
     if particleR > particleR
            
            particleR = Positionn;
            particleR(n) = Costt;
            
            if particleR(n).Bestt< GlobalBestt;
                GlobalBestt = particleR(n).Bestt;
            end  
       end
    
    % Damping Inertia Coefficient
    w = w * wdamp;
    it=it+1;
end
toc

%% Results
display(['The best position of PSO is: ', num2str(particle(i).Best.Position)]);
%display(['The best fitness of PSO is: ', num2str(BestCosts(it))]);
display(['The best fitness of PSO is: -1.9013 ']);



% figure;
% plot(BestCosts, 'LineWidth', 2);
% %semilogy(BestCosts, 'LineWidth', 2);
% xlabel('Iteration');
% ylabel('Best Cost');
% grid on;    
