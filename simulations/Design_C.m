classdef Design_C
    % Simulate using essentially Kock & Callot [2015] Experimental Design D
    % (which we here renamed) and retrieve coefficient matrix.
    
    properties
    end

    methods (Static)
        function Yminus0ton = sim_data(n,p,seed,r)
            % INPUTS:
            %   n:          Effective sample size, positive integer
            %   p:          System dimension, positive integer
            %   seed:       Seed for rng
            %   r:          MC iter, integer
            % OUTPUT:
            %   Yminus0ton: (1+n) x p outcome matrix 
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
            Yinit=zeros(p,1);
            Ylong=nan(p,nburn+1+n);
            Ylong(:,1)=Yinit;
            Eps=sigma_eps*randn(p,nburn+1+n);         % indep. eps_ti~N(0,sigma_eps^2)
            Theta=Design_C.coef_matrix(p);
            for t=2:nburn+1+n
                Ylong(:,t)=Theta*Ylong(:,t-1)+Eps(:,t);
            end
            Yminus0ton=Ylong(:,nburn+1:nburn+1+n)';   % return as (1+n) x p matrix 
        end

        function Theta = coef_matrix(p)
            Theta=nan(p,p);
            for i=1:p
                for j=1:p
                    Theta(i,j)=((-1)^(abs(i-j)))*(.4^(abs(i-j)+1));
                end
            end
        end
    end
end