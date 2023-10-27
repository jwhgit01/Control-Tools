function Mp = umat2aff(M,pNames)
%umat2aff 
%
% Copyright (c) 2022 Jeremy W. Hopwood. All rights reserved.
%
% This function ...
%
% Inputs:
%
%   in1     The first input...
%  
% Outputs:
%
%   out1    The first output...
%


% dimensions
n = size(M,1);
m = size(M,2);

% initialize variables
r = 1;
np = size(pNames,1);
MpCell = cell(np+1,1);
MpCell(:,:) = {zeros(n,m)};
Mp = cell2struct(MpCell,cat(1,'o',pNames),1);

% if matrix is not a umat, assign just the constant part
if isnumeric(M)
    Mp.o = M;
    return
end

% loop through elements of matrix
for ii = 1:n
    for jj = 1:m
        
        % is this element uncertain
        if size(fieldnames(M(ii,jj).Uncertainty),1) > 0
            
            % get the uncertain parameters in this block
            pijNames = fieldnames(M(ii,jj).Uncertainty);
            npij = size(pijNames);

            % we will need to set all parameters in a block to zero
            SijCell = cell(npij);
            SijCell(:,:) = {0};
            Sij = cell2struct(SijCell,pijNames,1);

            % loop through each ureal parameter in this block
            for k = 1:npij

                % increase the number of affine terms
                r = r + 1;
            
                % get the kth uncertain parameters in this block
                pijkName = pijNames{k};
                
                % store constant component
                Mp.o(ii,jj) = usubs(M(ii,jj),Sij);
                
                % store parameter-varying component
                Sijr1 = Sij;
                Sijr1.(pijkName) = 1; % set the rth parameter to 1
                Mp.(pijkName)(ii,jj) = usubs(M(ii,jj)-Mp.o(ii,jj),Sijr1);

            end
        
        else

            % store in M0
            Mp.o(ii,jj) = M(ii,jj).NominalValue;

        end
    end
end

end