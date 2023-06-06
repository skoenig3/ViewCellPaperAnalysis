function [f, df, hessian] = lnPoissonGlm(param,data,modelType,mCCVParms,numParams)
%modified from glm_model in updated spline repo by SDK 7/24/22
%split so can use L2 and roughness
%updated code with forloop for J_g/J_h using Michael's code to streamline
%for roughness.


%---Process Inputs---%
X = data{1}; % subset of A
Y = data{2}; % number of spikes
if isrow(Y)
    Y = Y' ;
end


%---Compute Firing Rate---%
% compute the firing rate
u = X * param;
rate = exp(u);


%---Get Hessian---%
switch mCCVParms.regularizationType
    case {'L1','Lasso'} % L1/Lasso regularizer weight
        %seth thinks this is right but not 100% sure
        disp('please check math!')
        
        % L1 regularizer weight
        alpha = 1; % also known as lambda, shoulde be 1!
        f_L1 = 0.5*alpha*sum(abs(param)); %L1 magnitude cost
        
        f = sum(rate-Y.*u) + f_L1; %loss function F
        df = real(X' * (rate - Y)) + alpha; %derviate of L1 is alpha
        rX = bsxfun(@times,rate,X);
        hessian = rX'*X + alpha*eye(numel(param));
        
        
    case {'L2','Ridge'} % L2/Ridge regularization weight
        
        % L2 regularizer weight
        alpha = 1; % also known as lambda, shoulde be 1!
        f_L2 = 0.5*alpha*sum(param.^2); %L2 squared cost
        
        f = sum(rate-Y.*u) + f_L2; %loss function F
        df = real(X' * (rate - Y) + alpha*param);
        rX = bsxfun(@times,rate,X);
        hessian = rX'*X + alpha*eye(numel(param));
        
    case {'roughness'}
        % roughness regularizer weight - note: these are tuned using the sum of f,
        % and thus have decreasing influence with increasing amounts of data
        b = mCCVParms.regularizationWeights;
        
        % initialize parameter-relevant variables
        n_var = length(b); % total number of variables
        J = zeros(n_var,1);
        J_g = cell(n_var,1);
        J_h = cell(n_var,1);
        
        % parse the parameters
        [allparams] = find_param(param,modelType,numParams);% global variable from before
        
        % compute the contribution for f, df, and the hessian
        for i = 1:length(allparams)
            if ~isempty(allparams{i})
                switch mCCVParms.typeParams{i}
                    case '2d'
                        [J(i),J_g{i},J_h{i}] = rough_penalty_2d(allparams{i},b(i));
                    case '1d'
                        [J(i),J_g{i},J_h{i}] = rough_penalty_1d(allparams{i},b(i));
                    case '1dcirc'
                        [J(i),J_g{i},J_h{i}] = rough_penalty_1d_circ(allparams{i},b(i));
                end
            end
        end
        
        
        % compute f, the gradient, and the hessian
        f = sum(rate-Y.*u) + sum(J);
        df = real(X' * (rate - Y) + vertcat(J_g{:}));
        rX = bsxfun(@times,rate,X);
        hessian = rX'*X + blkdiag(J_h{:});
        
    otherwise
        error('regularization type not supported!')
end

end


function [J,J_g,J_h] = rough_penalty_2d(param,beta)
numParam = numel(param);
D1 = spdiags(ones(sqrt(numParam),1)*[-1 1],0:1,sqrt(numParam)-1,sqrt(numParam));
DD1 = D1'*D1;
M1 = kron(eye(sqrt(numParam)),DD1);
M2 = kron(DD1,eye(sqrt(numParam)));
M = (M1 + M2);

J = beta*0.5*param'*M*param;
J_g = beta*M*param;
J_h = beta*M;
end


function [J,J_g,J_h] = rough_penalty_1d_circ(param,beta)

numParam = numel(param);
D1 = spdiags(ones(numParam,1)*[-1 1],0:1,numParam-1,numParam);
DD1 = D1'*D1;

% to correct the smoothing across first and last bin
DD1(1,:) = circshift(DD1(2,:),[0 -1]);
DD1(end,:) = circshift(DD1(end-1,:),[0 1]);

J = beta*0.5*param'*DD1*param;
J_g = beta*DD1*param;
J_h = beta*DD1;
end


function [J,J_g,J_h] = rough_penalty_1d(param,beta)

numParam = numel(param);
D1 = spdiags(ones(numParam,1)*[-1 1],0:1,numParam-1,numParam);
DD1 = D1'*D1;
J = beta*0.5*param'*DD1*param;
J_g = beta*DD1*param;
J_h = beta*DD1;
end


function [allparams] = find_param(param,modelType,numParams)
% function to find the right parameters given the model type

assert(length(numParams)== size(modelType,2)); % make sure the length matches

numParams(modelType==0) = 0; % only include the parameters which are in the model
end_idx = cumsum(numParams);

allparams = cell(length(numParams),1);
allparams{1} = param(1:end_idx(1));
for jj = 2:length(numParams)
    allparams{jj} = param(end_idx(jj-1)+1:end_idx(jj));
end
end