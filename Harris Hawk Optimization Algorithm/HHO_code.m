clc; clear; close all;

%% Problem Definition

CostFunction = @(x) Sphere(x); %Objective function

nVar=5;                        % Number of unknown (Decision) Variables

VarSize=[1 nVar];              % Matrix Size of Decision Variables

VarMin = -10;                  % Lower bound of Decision Variables
VarMax = 10;                   % Upper bound of Decision Variables

%% Parameters of HHO

N=30; % Number of search agents
T=500; % Maximum number of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

[Rabbit_Energy,Rabbit_Location,CNVG]=HHO(N,T,lb,ub,dim,fobj);

%% Initialization

Rabbit_Location=zeros(1,dim);
Rabbit_Energy=inf;

%Initialize the locations of Harris' hawks
X=initialization(N,dim,ub,lb);

CNVG=zeros(1,T);

t=0; % Loop counter

%% Main loop of HHO
while t<T
    for i=1:size(X,1)
            E0=2*rand()-1; %-1<E0<1
            Escaping_Energy=E1*(E0);  % escaping energy of rabbit

        if abs(Escaping_Energy)>=1
            %% Exploration phase
            if q<0.5
                X(i,:)=X_rand-rand()*abs(X_rand-2*rand()*X(i,:));    
            elseif q>=0.5
                X(i,:)=(Rabbit_Location(1,:)-mean(X))-rand()*((ub-lb)*rand()+lb);
            end

        elseif abs(Escaping_Energy)<1
             %% Exploitation phase
             r=rand(); % probablity of each event

             if r>=0.5 && abs(Escaping_Energy)<0.5 % Hard besiege
                X(i,:)=(Rabbit_Location)-Escaping_Energy*abs(Rabbit_Location-X(i,:));
             end

             if r>=0.5 && abs(Escaping_Energy)>=0.5  % Soft besiege
                Jump_strength=2*(1-rand()); % random jump strength of the rabbit
                X(i,:)=(Rabbit_Location-X(i,:))-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
             end

             if r<0.5 && abs(Escaping_Energy)>=0.5, % Soft besiege % rabbit try to escape by many zigzag deceptive motions

                Jump_strength=2*(1-rand());
                X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));

                if fobj(X1)<fobj(X(i,:)) % improved move?
                    X(i,:)=X1;
                else % hawks perform levy-based short rapid dives around the rabbit
                    X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:))+rand(1,dim).*Levy(dim);
                    if (fobj(X2)<fobj(X(i,:))), % improved move?
                        X(i,:)=X2;
                      end
                  end
             end

              if r<0.5 && abs(Escaping_Energy)<0.5, % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                % hawks try to decrease their average location with the rabbit
                Jump_strength=2*(1-rand());
                X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X));

                if fobj(X1)<fobj(X(i,:)) % improved move?
                    X(i,:)=X1;
                else % Perform levy-based short rapid dives around the rabbit
                    X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X))+rand(1,dim).*Levy(dim);
                    if (fobj(X2)<fobj(X(i,:))), % improved move?
                        X(i,:)=X2;
                    end
                 end
              end
        end
    end 
    t=t+1;
    CNVG(t)=Rabbit_Energy;
end 

function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end
%% Results

%Draw objective space
figure,
hold on
semilogy(CNVG,'Color','b','LineWidth',4);
title('Convergence curve')
xlabel('Iteration');
ylabel('Best fitness obtained so far');
axis tight
grid off
box on
legend('HHO')

display(['The best location of HHO is: ', num2str(Rabbit_Location)]);
display(['The best fitness of HHO is: ', num2str(Rabbit_Energy)]);