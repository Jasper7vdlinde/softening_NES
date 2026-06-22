%% Softening stiffness NES
% https://link.springer.com/article/10.1007/s11071-025-11164-6
% stiffness : kappa*atan(\bar{z})

clearvars;
%close all;
clc;
%% Define system

%Host system
w_i = 36;
xi = 0.2;
m = 1;
P =5; % Forcing

%NES 
xi_na = 0.2;
kappa = 2.5;
ep = 0.02;

Om = linspace(0.75,1.25,1000); % Frequency vector

for i=1:length(Om)
    
    sigma = (1-Om(i)^2)/(ep);
    
    %Polynomial relating NES amplitude and force
    a = Om(i)^4*(-sigma + Om(i)^2)^2 + Om(i)^6*xi^2- 2*(sigma - Om(i)^2)*Om(i)^4*xi*xi_na + (sigma - Om(i)^2)^2*Om(i)^2*xi_na^2 - 2*Om(i)^2*(-sigma + Om(i)^2)*(Om(i)^4 + Om(i)^2*xi*xi_na) + (Om(i)^4 + Om(i)^2*xi*xi_na)^2;
    b = -2*Om(i)^6*(-sigma+Om(i)^2)-4*kappa*Om(i)^2*(-sigma+Om(i)^2)^2+Om(i)^4*(-sigma+Om(i)^2)^2-4*kappa*Om(i)^4*xi^2+Om(i)^2^3*xi^2+4*kappa*(sigma-Om(i)^2)*Om(i)^2*xi*xi_na-2*(sigma-Om(i)^2)*Om(i)^2^2*xi*xi_na-2*Om(i)^4*(-sigma+Om(i)^2)*xi*xi_na+(sigma-Om(i)^2)^2*Om(i)^2*xi_na^2+2*Om(i)^4*(Om(i)^4+Om(i)^2*xi*xi_na)+4*kappa*(-sigma+Om(i)^2)*(Om(i)^4+Om(i)^2*xi*xi_na)+2*Om(i)^2*xi*xi_na*(Om(i)^4+Om(i)^2*xi*xi_na)-(Om(i)^4+Om(i)^2*xi*xi_na)^2;
    c= Om(i)^8+4*kappa*Om(i)^4*(-sigma+Om(i)^2)+4*kappa^2*(-sigma+Om(i)^2)^2-Om(i)^4*(-sigma+Om(i)^2)^2+4*kappa^2*Om(i)^2*xi^2-Om(i)^6*xi^2+2*(sigma-Om(i)^2)*Om(i)^4*xi*xi_na+2*Om(i)^6*xi*xi_na+4*kappa*Om(i)^2*(-sigma+Om(i)^2)*xi*xi_na-(sigma-Om(i)^2)^2*Om(i)^2*xi_na^2+Om(i)^4*xi^2*xi_na^2-2*Om(i)^4*(Om(i)^4+Om(i)^2*xi*xi_na)-4*kappa*(-sigma+Om(i)^2)*(Om(i)^4+Om(i)^2*xi*xi_na)+2*Om(i)^2*(-sigma+Om(i)^2)*(Om(i)^4+Om(i)^2*xi*xi_na)-2*Om(i)^2*xi*xi_na*(Om(i)^4+Om(i)^2*xi*xi_na)-Om(i)^4*P^2;
    d = -Om(i)^4*P^2-Om(i)^8-4*kappa*Om(i)^4*(-sigma+Om(i)^2)+2*Om(i)^6*(-sigma+Om(i)^2)-4*kappa^2*(-sigma+Om(i)^2)^2+4*kappa*Om(i)^2*(-sigma+Om(i)^2)^2-Om(i)^4*(-sigma+Om(i)^2)^2-4*kappa^2*Om(i)^2*xi^2+4*kappa*Om(i)^4*xi^2-Om(i)^6*xi^2-4*kappa*(sigma-Om(i)^2)*Om(i)^2*xi*xi_na+2*(sigma-Om(i)^2)*Om(i)^4*xi*xi_na-2*Om(i)^6*xi*xi_na-4*kappa*Om(i)^2*(-sigma+Om(i)^2)*xi*xi_na+2*Om(i)^4*(-sigma+Om(i)^2)*xi*xi_na-(sigma-Om(i)^2)^2*Om(i)^2*xi_na^2-Om(i)^4*xi^2*xi_na^2;
    
    [r] = roots([a,b,c,d]);
    r_real=r(imag(r)==0);
    r_real=r_real(r_real>1);
    
    % Single solution
    if(length(r_real) == 1)
          zna(i,1) = (r_real^2-1);
          zna(i,2)=NaN;
          zna(i,3)=NaN;
          zna(i,:) =sort(zna(i,:),'ascend') ;

          zo(i,1) = ((xi_na*Om(i))^2+  (Om(i)^2-2*kappa*(sqrt(zna(i,1)+1)-1)/(zna(i,1)))^2)*zna(i,1)/(Om(i))^4;
          zo(i,2)=NaN;
          zo(i,3)=NaN;
          %Phase
          alpha1(i,1) = asin( xi*sqrt(Om(i)^2)*sqrt(zo(i,1))/P); %- asin( (xi*sqrt(Om(i)^2)*sqrt(zo(i,1))+xi_na*Om(i)^2*sqrt(zna(i,1))/(sqrt(Om(i)^2)*sqrt(zo(i,1))))/P );
          beta1(i,1) = asin(xi_na* sqrt(zna(i,1)/ sqrt(zo(i,1))/sqrt(Om(i)^2)))-   alpha1(i,1);
           
           %stability
          M2= [-ep*sigma - 1i*xi*ep*sqrt(Om(i)^2), 0, 1i*ep*xi_na*sqrt(Om(i)^2)+ep*kappa/(sqrt(zna(i,1)+1)), ep*kappa*exp(1i*beta1(i,1))^2*(1/sqrt(zna(i,1)+1)-2*(sqrt(zna(i,1)+1)-1)/(zna(i,1)));
            0,  ep*sigma - 1i*xi*ep*sqrt(Om(i)^2), -conj(ep*kappa*exp(1i*beta1(i,1))^2*(1/sqrt(zna(i,1)+1)-2*(sqrt(zna(i,1)+1)-1)/(zna(i,1)))),  -conj(1i*ep*xi_na*sqrt(Om(i)^2)+ep*kappa/(sqrt(zna(i,1)+1)));
            ep*sigma + 1i*xi*ep*sqrt(Om(i)^2)+ Om(i)^2 , 0, -1i*xi_na*sqrt(Om(i)^2)*(1+ep)+Om(i)^2-(1+ep)*(kappa/(sqrt(zna(i,1)+1))),  -(1+ep)*kappa*exp(1i*beta1(i,1))^2*(1/sqrt(zna(i,1)+1)-2*(sqrt(zna(i,1)+1)-1)/(zna(i,1)));
            0,   -ep*sigma + 1i*xi*ep*sqrt(Om(i)^2) - Om(i)^2,  -conj(-(1+ep)*kappa*exp(1i*beta1(i,1))^2*(1/sqrt(zna(i,1)+1)-2*(sqrt(zna(i,1)+1)-1)/(zna(i,1)))), -conj(-1i*xi_na*sqrt(Om(i)^2)*(1+ep)+Om(i)^2-(1+ep)*(kappa/(sqrt(zna(i,1)+1))));]/(sqrt(Om(i)^2)*2*1i );
          M2eig(:,i)=eig(M2);
          if(any(real(eig(M2))>0))
              
           B_unstable(i,1) =sqrt(zna(i,1)); 
           A_unstable(i,1) =sqrt(zo(i,1)); 
           B_stable(i,1) = NaN; 
           A_stable(i,1) = NaN;
           
          else
              
           B_stable(i,1) =sqrt(zna(i,1)); 
           A_stable(i,1) =sqrt(zo(i,1));  
           B_unstable(i,1) = NaN; 
           A_unstable(i,1) = NaN;
           
           end
           B_stable(i,2) = NaN; 
           A_stable(i,2) = NaN;
           B_stable(i,3) = NaN; 
           A_stable(i,3) = NaN;
           B_unstable(i,2) = NaN; 
           A_unstable(i,2) = NaN;
           B_unstable(i,3) = NaN; 
           A_unstable(i,3) = NaN;
               
          
% three solution
    else
          zna(i,1) = (r_real(1)^2-1);
          zna(i,2)= (r_real(2)^2-1);
          zna(i,3)= (r_real(3)^2-1);
          zna(i,:) =sort(zna(i,:),'ascend') ;

          zo(i,1) = ((xi_na*Om(i))^2+  (Om(i)^2-2*kappa*(sqrt(zna(i,1)+1)-1)/(zna(i,1)))^2)*zna(i,1)/(Om(i))^4;
          zo(i,2) = ((xi_na*Om(i))^2+  (Om(i)^2-2*kappa*(sqrt(zna(i,2)+1)-1)/(zna(i,2)))^2)*zna(i,2)/(Om(i))^4;
          zo(i,3) = ((xi_na*Om(i))^2+  (Om(i)^2-2*kappa*(sqrt(zna(i,3)+1)-1)/(zna(i,3)))^2)*zna(i,3)/(Om(i))^4;

          for j=1:3
           alpha1(i,j) =  asin( xi*sqrt(Om(i)^2)*sqrt(zo(i,j))/P); %-asin( (xi*sqrt(Om(i)^2)*sqrt(zo(i,j))+xi_na*Om(i)^2*sqrt(zna(i,j))/(sqrt(Om(i)^2)*sqrt(zo(i,j))))/P );
           beta1(i,j) = asin(xi_na* sqrt(zna(i,j)/ sqrt(zo(i,j))/sqrt(Om(i)^2)))-   alpha1(i,j);
           
           M2= [-ep*sigma - 1i*xi*ep*sqrt(Om(i)^2), 0, 1i*ep*xi_na*sqrt(Om(i)^2)+ep*kappa/(sqrt(zna(i,j)+1)), ep*kappa*exp(1i*beta1(i,1))^2*(1/sqrt(zna(i,j)+1)-2*(sqrt(zna(i,j)+1)-1)/(zna(i,j)));
            0,  ep*sigma - 1i*xi*ep*sqrt(Om(i)^2), -conj(ep*kappa*exp(1i*beta1(i,1))^2*(1/sqrt(zna(i,j)+1)-2*(sqrt(zna(i,j)+1)-1)/(zna(i,j)))),  -conj(1i*ep*xi_na*sqrt(Om(i)^2)+ep*kappa/(sqrt(zna(i,j)+1)));
            ep*sigma + 1i*xi*ep*sqrt(Om(i)^2)+ Om(i)^2 , 0, -1i*xi_na*sqrt(Om(i)^2)*(1+ep)+Om(i)^2-(1+ep)*(kappa/(sqrt(zna(i,j)+1))),  -(1+ep)*kappa*exp(1i*beta1(i,1))^2*(1/sqrt(zna(i,j)+1)-2*(sqrt(zna(i,j)+1)-1)/(zna(i,j)));
            0,   -ep*sigma + 1i*xi*ep*sqrt(Om(i)^2) - Om(i)^2,  -conj(-(1+ep)*kappa*exp(1i*beta1(i,1))^2*(1/sqrt(zna(i,j)+1)-2*(sqrt(zna(i,j)+1)-1)/(zna(i,j)))), -conj(-1i*xi_na*sqrt(Om(i)^2)*(1+ep)+Om(i)^2-(1+ep)*(kappa/(sqrt(zna(i,j)+1))));]/(sqrt(Om(i)^2)*2*1i );
           
         eig(M2);
        
            if(any(real(eig(M2))>0))
               B_unstable(i,j) =sqrt(zna(i,j)); 
               A_unstable(i,j) =sqrt(zo(i,j)); 
               B_stable(i,j) = NaN; 
               A_stable(i,j) = NaN;
            else
                B_stable(i,j) =sqrt(zna(i,j)); 
               A_stable(i,j) =sqrt(zo(i,j));  
               B_unstable(i,j) = NaN; 
               A_unstable(i,j) = NaN;
            end
        end
    end
    
    %SIM computations (min/max and b++)
        a = Om(i)^4 + Om(i)^2*xi_na^2;
    b = -2*kappa*Om(i)^2 + 2*Om(i)^4 + 2*Om(i)^2*xi_na^2;
    c = -4*kappa*Om(i)^2 + Om(i)^4 + Om(i)^2*xi_na^2;
    d = 4*kappa^2 - 2*kappa*Om(i)^2;
     [r] = roots([a,b,c,d]);
    r_real=r(imag(r)==0);
    r_real=r_real(r_real>1);
         if(length(r_real)==2)
           b_max(i)= max(r_real.^2-1);
           b_min(i)= min(r_real.^2-1);
           a_min(i) = ((xi_na*Om(i))^2+ ( Om(i)^2-2*kappa*(sqrt(b_max(i)+1)-1)/(b_max(i)))^2)*b_max(i)/Om(i)^4;
           a_max(i) = ((xi_na*Om(i))^2+ ( Om(i)^2-2*kappa*(sqrt(b_min(i)+1)-1)/(b_min(i)))^2)*b_min(i)/Om(i)^4;
              Zb3 = -(-4*kappa^2 + 4*kappa*Om(i)^2 - Om(i)^4 - a_max(i)*Om(i)^4 - Om(i)^2*xi_na^2)/(Om(i)^4+Om(i)^2*xi_na^2)/((min(r_real))^2);      
     b3(i) = sqrt(Zb3^2-1); %b++
         end
    end

 zeta = ep*xi/2;

figure(1)
hold on
plot(Om,A_stable,'k','LineWidth',2)
hold on
plot(Om,A_unstable,'k--','LineWidth',2)
plot(Om,sqrt(a_max),'k:','LineWidth',1)
plot(Om,sqrt(a_min),'k:','LineWidth',1)
hold on
plot(Om,P*ep./sqrt( (1-Om.^2).^2+4*zeta^2*Om.^2),'k-.')
box on

   ax = gca; 
ax.FontSize = 15; 
xlabel('$\Omega$','FontSize',16.5,'interpreter','latex') 
ylabel('$a$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
xlim([min(Om) max(Om)])
figure(2)
hold on
hold on
plot(Om,B_stable,'k','LineWidth',2)
hold on
plot(Om,B_unstable,'k--','LineWidth',2)
plot(Om,sqrt(b_max),'k:','LineWidth',1)
plot(Om,sqrt(b_min),'k:','LineWidth',1)
plot(Om,b3,'r:','LineWidth',1)


   ax = gca; 
ax.FontSize = 15; 
xlabel('$\Omega$','FontSize',16.5,'interpreter','latex') 
ylabel('$b$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
xlim([min(Om) max(Om)])
box on

%Relative
figure(3)
hold on
plot(Om,A_stable/P,'k','LineWidth',2)
hold on
plot(Om,A_unstable/P,'k--','LineWidth',2)
plot(Om,sqrt(a_max)/P,'k:','LineWidth',1)
plot(Om,sqrt(a_min)/P,'k:','LineWidth',1)
hold on
plot(Om,ep./sqrt( (1-Om.^2).^2+4*zeta^2*Om.^2),'k-.')
box on

   ax = gca; 
ax.FontSize = 15; 
xlabel('$\Omega$','FontSize',16.5,'interpreter','latex') 
ylabel('$a$','FontSize',16.5,'interpreter','latex') 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
xlim([min(Om) max(Om)])



%% Time simulations, comment out as takes a while

% Time sim
%     Tl =5000;
% 
%     t=0:0.1:Tl;
%     N = length(t);
%     f = 0:10/N:(N-1)*10/N;
%     w = 0.9:0.002:1.1;
% 
% for i=1:length(w)
%     x0_real= 0;
% 
%     x0 = [0 0];
%     x0_dot = [0 0];
% 
%  f1 = @(t,y)[y(3);y(4);...
%     -(y(1) + ep*xi*y(3)- ep*kappa*atan(y(2)-y(1)) -ep*xi_na*(y(4)-y(3))-ep*P*cos(w(i)*t))/m;...
%     -(ep*kappa*atan(y(2)-y(1)) + ep*xi_na*(y(4)-y(3)))/ep;...
%     ];
% 
%     Prec = 1e-8;
%     % Actual numerical simulation
% 
%     options = odeset('RelTol',Prec,'AbsTol',[Prec Prec Prec Prec]);
%           [T1,Y1] = ode45(f1,t,[x0 x0_dot],options);
% 
%     % Plot simulation result      
% 
%     Y1filt = bandpass(Y1,[0.9*w(i)/(2*pi) 1.1*w(i)/(2*pi)],10,'ImpulseResponse','iir','Steepness',0.999999);
% %     A_main(i)=mean(envY(round(3*length(Y1)/4):end));
% % 
% %     B_main(i)=mean(envB(round(3*length(Y1)/4):end));
%     envY  = envelope(Y1(:,1),150,'peaks'); %sqrt(Y1(:,1).^2+Y1(:,3).^2);
% envB  = envelope(abs(Y1(:,1)-Y1(:,2)),150,'peaks'); %sqrt((Y1(:,1)-Y1(:,2)).^2+(Y1(:,3)-Y1(:,4)).^2);
% 
%      A_main(i)=rms(Y1(round(3*length(Y1)/4):end,1))*sqrt(2);
%     B_main(i)=rms(Y1(round(3*length(Y1)/4):end,2)-Y1(round(3*length(Y1)/4):end,1))*sqrt(2);
%     
%     TF = islocalmin(envY,'MinSeparation',100,'SamplePoints',T1);
%     if(length(TF)>0)
%             envY_min = envY(TF)
%      A_main_min(i)= min(envY_min(round(3*length(envY_min)/4):end))
%      if(A_main_min(i)/ A_main(i) > 0.9)
%          A_main_min(i)=NaN;
%      end
%                    % A_main_min(i)=min(envY(round(3*length(Y1)/4):end));
%                 A_main_max(i)=max(envY(round(3*length(Y1)/4):end));
%                     B_main_min(i)=min(envB(round(3*length(Y1)/4):end));
%                 B_main_max(i)=max(envB(round(3*length(Y1)/4):end));
%                 
%                  if(A_main_max(i)/ A_main(i) < 1.1)
%          A_main_max(i)=NaN;
%                  end
%      
%     else
%      A_main_min(i)= NaN;
%                    % A_main_min(i)=min(envY(round(3*length(Y1)/4):end));
%                 A_main_max(i)= NaN;
%                     B_main_min(i)= NaN;
%                 B_main_max(i)= NaN;
%     end
%     
% % 
% % 
% % 
% % % 
% figure(100)
% plot(T1,Y1(:,1))
% hold on
% plot(T1,envY,T1(TF),envY(TF),'r*')
% % 
% % xlim([500, 2000])
% %      pause(0.8)
% %  clf
% % figure(101)
% % plot(T1,Y1(:,2)-Y1(:,1))
% % hold on
% % plot(T1,envB)
% % xlim([500, 2000])
% % title(num2str(w(i)))
%      pause(0.2)
% clf
% 
%    i
% end
% 
% figure(1)
% % plot(w,A_main/(ep*m*omega0^2*F)*sqrt(2),'-o','LineStyle', 'none' )
%  hold on
% plot(w,A_main)%,'-o','LineStyle', 'none' )
%  hold on
% plot(w,A_main_min,'-o','LineStyle', 'none' )
%  hold on
% plot(w,A_main_max,'-o','LineStyle', 'none' )
% 
% figure(2)
% % hold on
% % plot(w,B_main*sqrt(2),'-o','LineStyle', 'none' )
% hold on
% plot(w,B_main,'-o','LineStyle', 'none' )
%  hold on
% plot(w,B_main_min,'-o','LineStyle', 'none' )
%  hold on
% plot(w,B_main_max,'-o','LineStyle', 'none' )
