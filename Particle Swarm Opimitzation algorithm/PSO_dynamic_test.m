clc; clear; close all;

m=1;
nPop = 3;   
w = 1;       
wdamp = 0.1% Damping Ratio of Inertia Coefficient
%wdamp =0.99  % Damping Ratio (NORMAL)
c1 = 2;      
c2 = 2;      
%CNVG=zeros(1,1000);

%% 1st level
T= 250; 
Function_name='F8';
[VarMin,VarMax,nVar,CostFunction]=Get_Functions_details(Function_name);
VarSize=[1 nVar]; 

% The Particle Template
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];

% Create population array
particle = repmat(empty_particle, nPop, 1);

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
BestCosts = zeros(T, 1);

for it=1:T
   
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
    
    
    % Store the Best Cost Value
    BestCosts(it) = GlobalBest.Cost;

    % Display Iteration Information
    disp(['Iteration ' num2str(m) ': Best Cost = ' num2str(abs(BestCosts(it)))]); 
    CNVG(m)=abs(BestCosts(it));
    % Damping Inertia Coefficient
    w = w * wdamp;
    m=m+1;
    
end

display(['The best fitness of PSO is: ', num2str(abs(BestCosts(it)))]);

%% 2nd level
R= 150; 
wdamp = 0.9% Damping Ratio of Inertia Coefficient
%wdamp =0.99  % Damping Ratio (NORMAL)
c1 = 2;      
c2 = 2;    
Function_name='F9';
[VarMin,VarMax,nVar,CostFunction]=Get_Functions_details(Function_name);
VarSize=[1 nVar]; 

% The Particle Template
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];

% Create population array
particle = repmat(empty_particle, nPop, 1);

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
BestCosts = zeros(R, 1);

for it=1:R
   
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
    
    
    % Store the Best Cost Value
    BestCosts(it) = GlobalBest.Cost;
    CNVG(m)=abs(BestCosts(it));
    % Display Iteration Information
    disp(['Iteration ' num2str(m) ': Best Cost = ' num2str(abs(BestCosts(it)))]); 

    % Damping Inertia Coefficient
    w = w * wdamp;
    m=m+1;
    
end

display(['The best fitness of PSO is: ', num2str(abs(BestCosts(it)))]);

%% 3rd level
Q= 350; 
wdamp = 0.09% Damping Ratio of Inertia Coefficient
Function_name='F7';
[VarMin,VarMax,nVar,CostFunction]=Get_Functions_details(Function_name);
VarSize=[1 nVar]; 

% The Particle Template
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];

% Create population array
particle = repmat(empty_particle, nPop, 1);

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
BestCosts = zeros(Q, 1);

for it=1:Q
   
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
    
    
    % Store the Best Cost Value
    BestCosts(it) = GlobalBest.Cost;
    CNVG(m)=abs(BestCosts(it));
    % Display Iteration Information
    disp(['Iteration ' num2str(m) ': Best Cost = ' num2str(abs(BestCosts(it)))]); 
        
    % Damping Inertia Coefficient
    w = w * wdamp;
    m=m+1;
    
end

display(['The best fitness of PSO is: ', num2str(abs(BestCosts(it)))]);

%% 4th level
Z= 250; 
Function_name='F8';
[VarMin,VarMax,nVar,CostFunction]=Get_Functions_details(Function_name);
VarSize=[1 nVar]; 

% The Particle Template
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];

% Create population array
particle = repmat(empty_particle, nPop, 1);

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
BestCosts = zeros(Z, 1);

for it=1:Z
   
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
    
    
    % Store the Best Cost Value
    BestCosts(it) = GlobalBest.Cost;

    % Display Iteration Information
    disp(['Iteration ' num2str(m) ': Best Cost = ' num2str(abs(BestCosts(it)))]); 
    CNVG(m)=abs(BestCosts(it));
    % Damping Inertia Coefficient
    w = w * wdamp;
    m=m+1;
    
end

display(['The best fitness of PSO is: ', num2str(abs(BestCosts(it)))]);

figure;
plot(CNVG, 'LineWidth', 2);
%semilogy(BestCosts, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;    
