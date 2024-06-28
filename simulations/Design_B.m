classdef Design_B
    % Simulate using essentially Kock & Callot [2015] Experimental Design B
    % and retrieve coefficient matrices (as one wide matrix).
    % Note: Since our dimensions (p) are all even, we use blocks of four
    % instead of the five in Kock & Callot.

    properties
    end

    methods (Static)
        function Yminus3ton = sim_data(n,p,seed,r)
            % INPUTS:
            %   n:          Effective sample size 
            %   p:          System dimension, divisible by four.  
            %   seed:       Seed for rng
            %   r:          MC iter
            % OUTPUT:
            %   Yminus3ton: (4+n) x p outcome matrix 
            if nargin < 2 || isempty(n) || isempty(p)
                error('Not enough input arguments.');
            end
            if nargin < 4 || isempty(seed) || isempty(r)
                rng('default');
            else
                rng(seed+r);
            end        
            sigma_eps=.1;                             % std(eps_0i) [common]
            nburn=10000;                              % burn-in period
            Yinit=zeros(p,4);
            Ylong=nan(p,nburn+4+n);
            Ylong(:,1:4)=Yinit;
            Eps=sigma_eps*randn(p,nburn+4+n);         % indep. eps_ti~N(0,sigma_eps^2)
            Theta=Design_B.coef_matrix(p);
            Theta1=Theta(1:p,1:p);
            Theta4=Theta(1:p,3*p+1:4*p);
            for t=5:nburn+4+n
                Ylong(:,t)=Theta1*Ylong(:,t-1)+Theta4*Ylong(:,t-4)+Eps(:,t);
            end
            Yminus3ton=Ylong(:,nburn+1:nburn+4+n)';   % return as (4+n) x p matrix
        end

        function Theta = coef_matrix(p)
            if mod(p,4)==0
                pover4=p/4;
            else
                warning('p must be divisible by four.');
            end
            Theta1=kron(eye(pover4),.15*ones(4,4));
            Theta4=kron(eye(pover4),-.1*ones(4,4));
            Theta=[Theta1,zeros(p,2*p),Theta4];
        end
    end
end