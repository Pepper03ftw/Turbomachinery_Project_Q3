function howell = read_howell_curves(paths)
%READ_HOWELL_CURVES Read and interpolate the digitized Howell charts.
%
% The supplied CSVs have no headers, so this function reads them robustly
% and exposes function handles:
%   howell.f_beta2(beta2_deg)
%   howell.Phi_Re(Re_over_1e5)
%   howell.Psi_sc(s_over_c)

fb = read_two_col(paths.f_beta2);
pr = read_two_col(paths.Phi_Re);
ps = read_two_col(paths.Psi_sc);

% Sort and build linear interpolants with extrapolation pinned to endpoints.
fb = sortrows(fb,1);
pr = sortrows(pr,1);
ps = sortrows(ps,1);

howell.f_beta2 = @(x) interp_clamped(fb(:,1), fb(:,2), x);
howell.Phi_Re  = @(x) interp_clamped(pr(:,1), pr(:,2), x);
howell.Psi_sc  = @(x) interp_clamped(ps(:,1), ps(:,2), x);

howell.raw.f_beta2 = fb;
howell.raw.Phi_Re  = pr;
howell.raw.Psi_sc  = ps;
end

function xy = read_two_col(file)
opts = detectImportOptions(file);
opts.VariableNamesLine = 0;
T = readtable(file, opts);
if width(T) ~= 2
    error('Expected 2 columns in %s', file);
end
xy = table2array(T);
if ~isnumeric(xy)
    xy = str2double(string(xy));
end
end

function yq = interp_clamped(x,y,xq)
x = x(:); y = y(:);
xq = double(xq);
yq = interp1(x,y,min(max(xq,min(x)),max(x)),'linear');
end
