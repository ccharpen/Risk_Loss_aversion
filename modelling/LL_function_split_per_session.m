% P is a "parameter" matrix, representing experimental design (gambles).  
% rows are trials, column 2 is the value of the sure option, 
% column 3 is the value of the potential loss,
% column 4 is the value of the potential gain, 
% column 8 is the binary choice variable: y value (1 if gamble
% chosen, 0 if sure option chosen).
      
function f = LL_function_split_per_session(x,P)

%transform parameters to make sure they are constrained between values that
%make sense, for ex
x(1) = exp(x(1)); % mu sess1 [0 Inf]
x(2) = exp(x(2)); % lambda sess1 [0 Inf]
x(3) = exp(x(3)); % rho sess1 [0 Inf] 
x(4) = exp(x(4)); % mu sess2 [0 Inf]
x(5) = exp(x(5)); % lambda sess2 [0 Inf]
x(6) = exp(x(6)); % rho sess2 [0 Inf] 
x(7) = exp(x(7)); % mu sess3 [0 Inf]
x(8) = exp(x(8)); % lambda sess3 [0 Inf]
x(9) = exp(x(9)); % rho sess3 [0 Inf]

Pg = zeros(length(P),1);

%sum log likelihood for each observation (trial)
for i=1:length(P)
    
    if P(i,10)==1 %baseline session
        Ugamble = 0.5*sign(P(i,4))*(abs(P(i,4)))^x(3) + 0.5*x(2)*sign(P(i,3))*(abs(P(i,3)))^x(3);
        Usure = sign(P(i,2))*(abs(P(i,2)))^x(3);
        xb = x(1)*(Ugamble - Usure);
    %x(1) is mu, x(2) is lambda, and x(3) is rho
    elseif P(i,10)==2 %anticipation session
        Ugamble = 0.5*sign(P(i,4))*(abs(P(i,4)))^x(6) + 0.5*x(5)*sign(P(i,3))*(abs(P(i,3)))^x(6);
        Usure = sign(P(i,2))*(abs(P(i,2)))^x(6);
        xb = x(4)*(Ugamble - Usure);
    elseif P(i,10)==3 %recovery session
        Ugamble = 0.5*sign(P(i,4))*(abs(P(i,4)))^x(9) + 0.5*x(8)*sign(P(i,3))*(abs(P(i,3)))^x(9);
        Usure = sign(P(i,2))*(abs(P(i,2)))^x(9);
        xb = x(7)*(Ugamble - Usure);
    end

    if xb<-709
      xb=-709;
    end
    if xb>36
      xb=36;
    end
    
    %if choice value is 1, use one part of likelihood contribution.
    if(P(i,8)==1)
      Pg(i) = (1+exp(-1*xb))^-1;                  
    %if choice value is 0, use other part of likelihood contribution    
    elseif(P(i,8)==0)
      Pg(i) = 1-(1+exp(-1*xb))^-1;          
    end
end

f=-sum(log(Pg)); %negative value of loglikelihood
          

