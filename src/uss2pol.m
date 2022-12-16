function [Poly,legacy] = uss2pol(usys)
%
% Rules:
%   - system must be affine in the parameters
%   - more than one parameter may appear in a matrix element

% dimensions
nx = size(usys.A,1);
nu = size(usys.B,2);
ny = size(usys.C,1);
np = size(fieldnames(usys.Uncertainty),1);

% get the names and ranges of these parameters
unc = struct2cell(usys.Uncertainty);
pNames = cell(np,1);
range = cell(np,1);
for ii = 1:np
    pNames{ii} = unc{ii}.Name;
    range{ii} = unc{ii}.Range;
end

% convert ureal matrices to affine form
Ap = umat2aff(usys.A,pNames);
Bp = umat2aff(usys.B,pNames);
Cp = umat2aff(usys.C,pNames);
Dp = umat2aff(usys.D,pNames);

% loop though parameters and create ltisys "objects"
affSysCell = cell(np+1,1);
affSysCell(:,:) = {zeros(nx+ny+1,nx+nu+1)};
affSys = cell2struct(affSysCell,cat(1,'o',pNames));
affSys.o = ltisys(Ap.o,Bp.o,Cp.o,Dp.o);
for rr = 1:np
    affSys.(pNames{rr}) = ltisys(Ap.(pNames{rr}),...
                                Bp.(pNames{rr}),...
                                Cp.(pNames{rr}),...
                                Dp.(pNames{rr}));
end

% specify range of uncertain parameters for legacy functions
pv = pvec('box',cell2mat(range));

% create affine system using legacy functions
affSysCell = struct2cell(affSys);
affSysList = [affSysCell{:}];
affSys = psys(pv,affSysList);

% transform to polytopic system using legacy functions
legacy.PV = pv;
legacy.Vertex = polydec(pv);
legacy.Sys = aff2pol(affSys);

% put the legacy result into a readable format
Poly.Vertex = legacy.Vertex;
nv = size(Poly.Vertex,2);
Poly.Sys = cell(1,nv);
for rr = 1:nv
    sysMat = real(legacy.Sys(1:nx+ny,(2*rr)+(rr-1)*(nx+nu)+1:(2*rr)+(rr)*(nx+nu)));
    Ar = sysMat(1:nx,1:nx);
    Br = sysMat(1:nx,nx+1:nx+nu);
    Cr = sysMat(nx+1:nx+ny,1:nx);
    Dr = sysMat(nx+1:nx+ny,nx+1:nx+nu);
    Poly.Sys{1,rr} = ss(Ar,Br,Cr,Dr);
end

end