% P is a "parameter" matrix, representing experimental design (gambles).  
% rows are trials, column 2 is the value of the sure option, 
% column 3 is the value of the potential loss,
% column 4 is the value of the potential gain, 
% column 8 is the binary choice variable: y value (1 if gamble
% chosen, 0 if sure option chosen).
      
function Pgen = generate_choice_alldata(x,P)

Pg      = NaN(length(P),1);
Ugamble = NaN(length(P),1);
Usure   = NaN(length(P),1);
pred_ch = NaN(length(P),1);
ll      = NaN(length(P),1);
corr    = NaN(length(P),1);
ll_S    = NaN(length(P),1);

%sum log likelihood for each observation (trial)
for i=1:length(P)
    
    Ugamble(i) = 0.5*sign(P(i,4))*(abs(P(i,4)))^x(3) + 0.5*x(2)*sign(P(i,3))*(abs(P(i,3)))^x(3);
    Usure(i) = sign(P(i,2))*(abs(P(i,2)))^x(3);
    xb = x(1)*(Ugamble(i) - Usure(i));
    %x(1) is mu, x(2) is lambda, and x(3) is rho

    if xb<-709
      xb=-709;
    end
    if xb>36
      xb=36;
    end
    
    Pg(i) = (1+exp(-1*xb))^-1;
    
    n=rand();
    if n<Pg(i)
        pred_ch(i) = 1; %predicted choice is gamble   
        ll(i) = Pg(i); %loglikelihood
    else
        pred_ch(i) = 0; %predicted choice is right
        ll(i) = 1-Pg(i); %loglikelihood
    end 
    if pred_ch(i) == P(i,8) %if predicted choice is subject's choice
        corr(i) = 1;
    else
        corr(i) = 0;
    end
    if P(i,8)==1 %calculate likelihood of predicting subject choice
        ll_S(i) = Pg(i);
    elseif P(i,8)==0
        ll_S(i) = 1-Pg(i);
    end

end
Pgen = [Ugamble Usure Pg pred_ch ll corr ll_S];
          

