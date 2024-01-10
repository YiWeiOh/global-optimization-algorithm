clear all; close all; clc

N=3;

Function_name='F8'; 

T=50; 

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);
[Rabbit_Energy,Rabbit_Location,CNVG]=IHHO(N,T,lb,ub,dim,fobj);
% [Rabbit_Energy,Rabbit_Location,CNVG]=IHHO(N,lb,ub,dim,fobj);

%Draw objective space
% figure,
% hold on
% semilogy(CNVG,'Color','b','LineWidth',4);
% title('Convergence curve')
% xlabel('Iteration');  
% ylabel('Best fitness obtained so far');
% axis tight
% grid off
% box on
% legend('HHO')

display(['The best location of IHHO is: ', num2str(Rabbit_Location)]);
%display(['The best fitness of IHHO is: ', num2str(Rabbit_Energy)]);
display(['The best fitness of IHHO is: ', num2str(abs(Rabbit_Energy))]);


% display([num2str(Rabbit_Location)]);
% disp(['Iteration ' num2str(t) ': Best Cost = ' num2str(Rabbit_Energy)]);
% disp(['Iteration ' num2str(t) ': Best Cost = ' num2str(abs(Rabbit_Energy))]);
% display([num2str(Rabbit_Energy))]);
