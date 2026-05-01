function curves = read_lieblein_curves(paths)
%READ_LIEBLEIN_CURVES Read digitized Lieblein incidence/deviation curves.
% Interpolation is linear and extrapolation is disabled. Queries outside the
% digitized beta/solidity/t-c ranges return NaN.

if nargin < 1 || isempty(paths)
    paths = struct();
end
if ischar(paths) || isstring(paths)
    tmp.folder = char(paths);
    paths = tmp;
end
if ~isfield(paths,'folder') || isempty(paths.folder)
    paths.folder = 'lieblein_digitized';
end

files = struct();
files.i0_10       = get_path(paths,'i0_10','lieblein_i0_10.csv');
files.n           = get_path(paths,'n','lieblein_n.csv');
files.Ki_t        = get_path(paths,'Ki_t','lieblein_Ki_t.csv');
files.delta0_10   = get_path(paths,'delta0_10','lieblein_delta0_10.csv');
files.m_NACA65    = get_path(paths,'m_NACA65','lieblein_m_NACA65.csv');
files.Kdelta_t    = get_path(paths,'Kdelta_t','lieblein_Kdelta_t.csv');
files.b           = get_path(paths,'b','lieblein_b.csv');

T_i0 = readtable(files.i0_10);
T_n  = readtable(files.n);
T_Ki = readtable(files.Ki_t);
T_d0 = readtable(files.delta0_10);
T_m  = readtable(files.m_NACA65);
T_Kd = readtable(files.Kdelta_t);
T_b  = readtable(files.b);

F_i0 = scatteredInterpolant(T_i0.beta1_deg, T_i0.sigma, T_i0.i0_10_deg, 'linear', 'none');
F_n  = scatteredInterpolant(T_n.beta1_deg,  T_n.sigma,  T_n.n_coeff,    'linear', 'none');
F_d0 = scatteredInterpolant(T_d0.beta1_deg, T_d0.sigma, T_d0.delta0_10_deg, 'linear', 'none');

curves = struct();
curves.i0_10       = @(beta1_deg,sigma) interp2_no_extrap(F_i0,beta1_deg,sigma);
curves.n           = @(beta1_deg,sigma) interp2_no_extrap(F_n,beta1_deg,sigma);
curves.delta0_10   = @(beta1_deg,sigma) interp2_no_extrap(F_d0,beta1_deg,sigma);
curves.Ki_t        = @(t_over_c) interp1_no_extrap(T_Ki.t_over_c_max,T_Ki.Ki_t,t_over_c);
curves.m_NACA65    = @(beta1_deg) interp1_no_extrap(T_m.beta1_deg,T_m.m_coeff,beta1_deg);
curves.Kdelta_t    = @(t_over_c) interp1_no_extrap(T_Kd.t_over_c_max,T_Kd.Kdelta_t,t_over_c);
curves.b           = @(beta1_deg) interp1_no_extrap(T_b.beta1_deg,T_b.b_exp,beta1_deg);
curves.files = files;
curves.raw = struct('i0_10',T_i0,'n',T_n,'Ki_t',T_Ki,'delta0_10',T_d0, ...
                    'm_NACA65',T_m,'Kdelta_t',T_Kd,'b',T_b);
curves.ranges = struct( ...
    'i0_sigma',[min(T_i0.sigma) max(T_i0.sigma)], ...
    'i0_beta',[min(T_i0.beta1_deg) max(T_i0.beta1_deg)], ...
    'n_sigma',[min(T_n.sigma) max(T_n.sigma)], ...
    'n_beta',[min(T_n.beta1_deg) max(T_n.beta1_deg)], ...
    'delta0_sigma',[min(T_d0.sigma) max(T_d0.sigma)], ...
    'delta0_beta',[min(T_d0.beta1_deg) max(T_d0.beta1_deg)], ...
    'Ki_t_tc',[min(T_Ki.t_over_c_max) max(T_Ki.t_over_c_max)], ...
    'Kdelta_t_tc',[min(T_Kd.t_over_c_max) max(T_Kd.t_over_c_max)], ...
    'm_beta',[min(T_m.beta1_deg) max(T_m.beta1_deg)], ...
    'b_beta',[min(T_b.beta1_deg) max(T_b.beta1_deg)]);

    function p = get_path(s, field, default_name)
        if isfield(s,field) && ~isempty(s.(field))
            p = s.(field);
        else
            p = fullfile(s.folder, default_name);
        end
        if ~isfile(p)
            here = fileparts(mfilename('fullpath'));
            p2 = fullfile(here, p);
            if isfile(p2)
                p = p2;
            else
                error('Lieblein curve file not found: %s', p)
            end
        end
    end
end

function v = interp2_no_extrap(F, beta1_deg, sigma)
[bb,ss] = compatible_query(beta1_deg, sigma);
v = F(bb,ss);
end

function v = interp1_no_extrap(x, y, xq)
[xs,idx] = sort(x(:));
ys = y(idx);
[xu,ia] = unique(xs,'stable');
yu = ys(ia);
v = interp1(xu, yu, xq, 'linear', NaN);
end

function [a,b] = compatible_query(a,b)
if isscalar(a) && ~isscalar(b)
    a = a + zeros(size(b));
elseif isscalar(b) && ~isscalar(a)
    b = b + zeros(size(a));
elseif ~isequal(size(a),size(b))
    error('Lieblein interpolation queries must be scalar or equal-sized arrays.');
end
end
