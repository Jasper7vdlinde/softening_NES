clear all
close all
%% Add functions
 addpath('../');


%% SDOF system

m=1;
ep = 0.02;
xi_na =0.2;
xi = 0.2;

P = 9;

kappa= 5;


%% Nonlinear absorber
% Simulation parameters
% robustness
    Tl =2000;

    t=0:0.05:Tl;
    N = length(t);
    f = 0:10/N:(N-1)*10/N;
    
    w = 0.98; % Forcing frequency

    x0_real= 0;

    x0 = [0 0];
    x0_dot = [0 0];

 f1 = @(t,y)[y(3);y(4);...
    -(y(1) + ep*xi*y(3) - ep*kappa*atan(y(2)-y(1))-ep*xi_na*(y(4)-y(3))-ep*P*cos(w*t))/m;...
    -(ep*kappa*atan(y(2)-y(1))  + ep*xi_na*(y(4)-y(3)))/ep;...
    ];

    Prec = 1e-8;
    % Actual numerical simulation

    options = odeset('RelTol',Prec,'AbsTol',[Prec Prec Prec Prec]);
          [T1,Y1] = ode45(f1,t,[x0 x0_dot],options);
 

    Y1filt = bandpass(Y1,[0.90/(2*pi) 1.1/(2*pi)],10,'ImpulseResponse','iir','Steepness',0.999999);
envY  = envelope(Y1(:,1),200,'peaks'); 
envB  = envelope(abs(Y1(:,1)-Y1(:,2)),200,'peaks');

figure
box on

subplot(2,1,1);
plot(T1,Y1(:,1),'k','LineWidth',2)
hold on
plot(T1,envY,'LineWidth',2)
%xlim([1000, 1300])
   ax = gca; 
ax.FontSize = 15; 
ylabel('$A$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
 
subplot(2,1,2);
box on

plot(T1,Y1(:,2)-Y1(:,1),'k','LineWidth',2)
hold on
plot(T1,envB,'LineWidth',2)
%xlim([1000, 1300])
   ax = gca; 
ax.FontSize = 15; 
ylabel('$B$','FontSize',16.5,'interpreter','latex') 
xlabel('$\tau$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis

     Zb = 0:0.01:400;  
     Za = (Zb.*(w^2*xi_na^2+(w^2-2*kappa*(sqrt(Zb+1)-1)./(Zb)).^2))/(w^4);

     phase = asin(xi_na*w*sqrt(Zb)./(w^2*sqrt(Za)));
     
             a = w^4 + w^2*xi_na^2;
    b = -2*kappa*w^2 + 2*w^4 + 2*w^2*xi_na^2;
    c = -4*kappa*w^2 + w^4 + w^2*xi_na^2;
    d = 4*kappa^2 - 2*kappa*w^2;
         [r] = roots([a,b,c,d]);
    r_real=r(imag(r)==0);
    r_real=r_real(r_real>1);
         if(length(r_real)==2)
           b_max= max(r_real.^2-1);
           b_min= min(r_real.^2-1);
           a_min = ((xi_na*w)^2+ ( w^2-2*kappa*(sqrt(b_max+1)-1)/(b_max))^2)*b_max/w^4;
           a_max = ((xi_na*w)^2+ ( w^2-2*kappa*(sqrt(b_min+1)-1)/(b_min))^2)*b_min/w^4;
         end
        
    
     for i=2:length(Zb)
            m11= stability_m11( sqrt(Zb(i))*exp(1i*phase(i)),xi_na,kappa,w^2,xi);
            m22 = conj(m11);
            m12 = stability_m12( sqrt(Zb(i))*exp(1i*phase(i)),xi_na,kappa,w^2,xi);
            m21 = conj(m12);
            M = [m11, m12;
               m21, m22];
         if(any(real(eig(M))>0))
               B_unstable(i) =sqrt(Zb(i)); 
               A_unstable(i) =sqrt(Za(i)); 
               B_stable(i) = NaN; 
               A_stable(i) = NaN;
            else
           B_stable(i) =sqrt(Zb(i)); 
               A_stable(i) =sqrt(Za(i));  
               B_unstable(i) = NaN; 
               A_unstable(i) = NaN;
            end
     end
     
     
     Zb3 = -(-4*kappa^2 + 4*kappa*w^2 - w^4 - a_max*w^4 - w^2*xi_na^2)/(w^4+w^2*xi_na^2)/((min(r_real))^2);
     b3 = sqrt(Zb3-1);

%%
figure(501)
box on
hold on
plot(B_stable,A_stable,'k','LineWidth',2 )
plot(B_unstable,A_unstable,'k--','LineWidth',2 )
% plot(w,A_main/(ep*m*omega0^2*F)*sqrt(2),'-o','LineStyle', 'none' )
plot(envB,envY,'r' ,'LineWidth',2 )

   ax = gca; 
ax.FontSize = 15; 
xlabel('$b$','FontSize',16.5,'interpreter','latex') 
ylabel('$a$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
%legend('$\Omega=1$','$\Omega=0.965$','$\Omega=0.965$ $x(0)''=0.5$','Interpreter','latex')

%% SIMS STABILITY (DIFFERENT FROM FRF STABILITY)
function m11  = stability_m11(B,xi_na,kappa,X,xi)
        A = conj(B);
           m11  = -kappa/(sqrt(abs(B)^2+1));
          m11 = (-1i*X - xi_na*sqrt(X) -1i*m11);

       % m11 = (-1i*X - xi_na*sqrt(X) -1i*m11-1i*X*X/(sigma+1i*xi*sqrt(X)-X))/( sqrt(X)*2*(sigma+1i*xi*sqrt(X)) )*(sigma+1i*xi*sqrt(X)-X);
end

function m12  = stability_m12(B,xi_na,kappa,X,xi)
    A = conj(B);
       m12  = -kappa*(1/sqrt(abs(B)^2+1)-2*(sqrt(abs(B)^2+1)-1)/(abs(B)^2));
       m12 = -1i*m12
   % m12 = -1i*m12/( sqrt(X)*2*(sigma+1i*xi*sqrt(X)) )*(sigma+1i*xi*sqrt(X)-X);
end

