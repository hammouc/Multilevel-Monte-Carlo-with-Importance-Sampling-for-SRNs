function [t,X]=SSA(DefineReactions,tspan, IC, params)
%Implementation of Gillespie Exact Stochastic Algorithm.
% S = stoichiometry of C substrates in R reactions
% P = stoichiometry of products
% K = vector of reaction rates
% INITIALIZE by calling the passed function to define the reactions:
if ~isa(DefineReactions, 'function_handle') ...
        disp('pass a function handle with your reaction defs'); return; end
[S, P, K] = DefineReactions(params); % defines reaction stoichiometry
[R, C] = size(S);
if size(S) ~= size(P); disp('reaction defs inconsistent'); return;end
if size(IC,1) ~= 1; IC = IC'; end  % force row vector
if length(IC) ~= C; disp('ICs inconsistent with reactions'); return; end
if size(K,2) ~= 1; K = K'; end  % force row vector
if length(K) ~= R; disp ('parameters inconsistent with reactions');return;end
if length(tspan) ~= 2; disp('Time vector: [dt;maxtime] needed'); return;end;
%set the uniform random number generator 
rand('state',sum(100*clock));  
nRC  = 0;  %reaction counter
t(1) = tspan(1); 
maxT = tspan(2);
X = IC;  % initialize the matrix of chemicals 
%------------------------------------------------------------
while (t(end)<maxT )
    %step 1: Calculate a's (reaction rates given system state) 
    a = K;
    for r = 1:R
        for c = 1:C
            if S(r,c) == 0; a(r) = a(r);
            elseif S(r,c) == 1; a(r) = a(r)*X(end,c);
            elseif S(r,c) == 2; 
                a(r) = a(r)*X(end,c)*(X(end, c)-1)/2;
            elseif S(r,c) == 3; 
                a(r) = a(r)*X(end,c)*(X(end, c)-1)/2*(X(end,c)-2)/3;
            end; 
        end
    end 
    a0 = sum(a); % a0 is the total rate of change of system
    if a0 == 0 ;% system can't change; finish and exit;
        X(end+1,:) = X(end,:);
        t =[t; maxT];
        break; 
    end 

    %Step 2: calculate tau and r using random number generators
        % determine time of next reaction:
        
        p1  = rand;  tau = (1/a0)*log(1/p1);
        % determine which next reaction is:
        p2 = rand;
        for r=1:R
            if (sum(a(1:r)) >= p2*a0); break; end 
        end
    %Step 3: carry out the reaction
        t =[t; t(end) + tau]; % t is time array; add last entry to it.
        nRC = nRC + 1 ;       % nRC is number of reactions so far.
        X(end+1,:)=X(end,:)-S(r,:)+P(r,:);
end %end of while (t(end)<tspan(2) & any(X))
end % end of function

