% This code computes the phase function and the isochrons of the Van der Pol model, using forward integration. It runs in about 2 minutes on a standard computer.
% The code can be easily adapted to other models with a limit cycle. The method also works with models of higher dimensions.
% Note that the use of the "contour" matlab function is delicate with the phase function, due to the singularity between -pi and pi. We rather use the "contourcs" function (Copyright (c)2010, Takeshi Ikuma).

% The method is based on the results presented in "A. Mauroy and I. Mezic, On the use of Fourier averages to compute the global isochrons of (quasi)periodic dynamics,
% Chaos, vol. 22, no. 3, p. 033112, 2012"

% For more information, please email me at "alexandre.mauroy 'at' unamur.be"

% Written by Alexandre Mauroy

close all;
clear all;

% time horizon for the averaging
t_end=200;

% 2D grid in state space
x_interv=linspace(-3,3,40);
y_interv=linspace(-3,3,40);
[x0,y0]=meshgrid(x_interv,y_interv);

% model parameters
options = odeset('RelTol',1e-6,'AbsTol',1e-50);
param.mu=1;
f_dyn     =   @(t,X) van_der_pol(t,X,param);
% frequency of the limit cycle (van der pol model with mu=1)
omega=0.942958; % vdp mu=1

%%

number_traject_per_loop=100;
number_traject=numel(x0);
number_loops=ceil(number_traject/number_traject_per_loop);

average=zeros(1,number_traject);

point=[x0(:) y0(:)]; %initial conditions of the trajectories

% 100 trajectories are integrated simultaneously
for test=1:number_loops;
    
    % setting the initial conditions
    index_low=1+(test-1)*number_traject_per_loop;
    index_high=min(test*number_traject_per_loop,number_traject);
    init_cond=[point(index_low:index_high,1);point(index_low:index_high,2)];
    number_traj_loop=index_high-index_low+1;
    
    % ode integration
    [t,x]=ode45(f_dyn,[0 t_end],init_cond,options);
    
    % we choose the observable f(x1,x2)=x1
    f=x(:,1:number_traj_loop);

    % averaging
    average(index_low:index_high)=1/t_end*trapz(t,f.*(exp(-1i*t*omega)*ones(1,number_traj_loop)));
     
    clear x t f
    
end

data=reshape(average,size(x0));

%% figure

figure(2)
h=colormap('hsv');
h=pcolor(x0,y0,angle(data));
set(h,'LineStyle','none')
[t,x]=ode45(f_dyn,[0 t_end],[1 1],options);
hold on
box on
plot(x(end-400:end,1),x(end-400:end,2),'k--','Linewidth',2)

% plot of the contours (remove the contours created by the artificial discontinuity between -pi and pi)

value_phases=[-180:10:180];
S=contourcs(x0(1,:),y0(:,1),angle(data)*180/pi,value_phases);

    for j=1:length(S)
 
        test_sin=interp2(x0,y0,sin(angle(data)),S(j).X,S(j).Y);
        test_cos=interp2(x0,y0,cos(angle(data)),S(j).X,S(j).Y);
        test=atan2(test_sin,test_cos)*180/pi;
  
        S(j).X(find(abs(test-S(j).Level)>2.5))=[];
        S(j).Y(find(abs(test-S(j).Level)>2.5))=[];

        plot(S(j).X,S(j).Y,'k','Linewidth',1)
    
    end
    
% particular case: 180
    
S=contourcs(x0(1,:),y0(:,1),sin(angle(data)),[0 0]);

for j=1:length(S)

    test_cos=interp2(x0,y0,cos(angle(data)),S(j).X,S(j).Y);   
  
    S(j).X(find(abs(test_cos-(-1))>1))=[];
    S(j).Y(find(abs(test_cos-(-1))>1))=[];

    plot(S(j).X,S(j).Y,'k','Linewidth',1)
    
end

colorbar
