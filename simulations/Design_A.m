classdef Design_A
    
    % Simulate with Kock & Callot's [2015 J Economet] Experimental Design A
    % and retrieve coefficient matrix.

    properties
    end

    methods (Static)
    
        function Y0ton = sim_data(n,p,seed,r)
            % INPUTS:
            %   n:          Effective sample size, integer
            %   p:          System dimension, divisible by four
            %   seed:       Seed for rng
            %   r:          MC iter
            % OUTPUT:
            %   Y0ton:      (1+n) x p outcome matrix                                                        
            if nargin < 2 || isempty(n) || isempty(p)
                error('Not enough input arguments.');
            end
            if nargin < 3 || isempty(seed) || isempty(r)
                rng('default');
            else
                rng(seed+r);
            end
            sigma_eps=.1;                       % std(eps_0i) [common]
            theta=.5;                           % AR(1) coef [common]
            sigma_y=sigma_eps/sqrt(1-theta^2);  % => std(y_0i) [common]
            Yinit=sigma_y*randn(p,1);           % sample from stationary distribution
            Y0ton=nan(p,1+n);
            Y0ton(:,1)=Yinit;                   % initiate from stationary distribution
            Eps=sigma_eps*randn(p,1+n);         % independent eps_ti~N(0,sigma_eps^2)
            % Note: 1st epsilon (i.e. for t=0) not actually in use
            Theta=Design_A.coef_matrix(p);      % specify coef matrix (diagonal)
            for t=2:1+n
                Y0ton(:,t)=Theta*Y0ton(:,t-1)+Eps(:,t);
            end
            Y0ton=Y0ton'; % return as (1+n) x p matrix
        end

        function Theta = coef_matrix(p)
            theta=.5;
            Theta=theta*eye(p);
        end

    end
end