function row = lieblein_row_estimate(curves, row_type, chi_in_deg, chi_out_deg, sigma, tmax_over_c, opts)
%LIEBLEIN_ROW_ESTIMATE Estimate optimum incidence/deviation using digitized Lieblein charts.
% chi_in_deg and chi_out_deg are positive row-frame flow-angle magnitudes.

if nargin < 7 || isempty(opts)
    opts = struct();
end
family = get_opt(opts,'profile_family','NACA65');
family_key = upper(strrep(family,'-',''));

switch family_key
    case {'NACA65','NACA65SERIES'}
        Ki_sh = get_opt(opts,'Ki_sh_NACA65',1.1);
        Kdelta_sh = get_opt(opts,'Kdelta_sh_NACA65',1.1);
        m_note = 'm coefficient uses digitized NACA-65 curve.';
    case {'DCA','DOUBLECIRCULARARC'}
        Ki_sh = get_opt(opts,'Ki_sh_DCA',0.7);
        Kdelta_sh = get_opt(opts,'Kdelta_sh_DCA',0.75);
        m_note = 'm coefficient still uses the digitized NACA-65 curve; no separate DCA m curve was supplied.';
    otherwise
        Ki_sh = get_opt(opts,'Ki_sh_NACA65',1.1);
        Kdelta_sh = get_opt(opts,'Kdelta_sh_NACA65',1.1);
        m_note = 'unknown family; using NACA-65 shape factors by default.';
end

chi_in_deg = abs(chi_in_deg);
chi_out_deg = abs(chi_out_deg);
flow_deflection_deg = abs(chi_in_deg - chi_out_deg);
tmax_over_c = max(tmax_over_c, 0);

vals.i0_10_deg     = curves.i0_10(chi_in_deg, sigma);
vals.n_coeff       = curves.n(chi_in_deg, sigma);
vals.Ki_t          = curves.Ki_t(tmax_over_c);
vals.delta0_10_deg = curves.delta0_10(chi_in_deg, sigma);
vals.m_coeff       = curves.m_NACA65(chi_in_deg);
vals.Kdelta_t      = curves.Kdelta_t(tmax_over_c);
vals.b_exp         = curves.b(chi_in_deg);

vec = [vals.i0_10_deg vals.n_coeff vals.Ki_t vals.delta0_10_deg vals.m_coeff vals.Kdelta_t vals.b_exp];
valid = all(isfinite(vec)) && sigma > 0;

if valid
    i0_deg = Ki_sh * vals.Ki_t * vals.i0_10_deg;
    delta0_deg = Kdelta_sh * vals.Kdelta_t * vals.delta0_10_deg;
    denom = 1 + vals.n_coeff - vals.m_coeff/(sigma^vals.b_exp);
    if abs(denom) < 1e-8
        valid = false;
        theta_deg = NaN;
        incidence_deg = NaN;
        deviation_deg = NaN;
    else
        theta_deg = (flow_deflection_deg - i0_deg + delta0_deg)/denom;
        incidence_deg = i0_deg + vals.n_coeff*theta_deg;
        deviation_deg = delta0_deg + vals.m_coeff/(sigma^vals.b_exp)*theta_deg;
    end
else
    i0_deg = NaN;
    delta0_deg = NaN;
    denom = NaN;
    theta_deg = NaN;
    incidence_deg = NaN;
    deviation_deg = NaN;
end

metal_in_deg = chi_in_deg - incidence_deg;
metal_out_deg = chi_out_deg - deviation_deg;
stagger_deg = 0.5*(metal_in_deg + metal_out_deg);

row = struct('row_type',row_type, 'valid',valid, 'profile_family',family, ...
    'chi_in_deg',chi_in_deg, 'chi_out_deg',chi_out_deg, ...
    'flow_deflection_deg',flow_deflection_deg, 'sigma',sigma, ...
    'tmax_over_c',tmax_over_c, 'Ki_sh',Ki_sh, 'Kdelta_sh',Kdelta_sh, ...
    'i0_10_deg',vals.i0_10_deg, 'n_coeff',vals.n_coeff, 'Ki_t',vals.Ki_t, ...
    'delta0_10_deg',vals.delta0_10_deg, 'm_coeff_NACA65',vals.m_coeff, ...
    'Kdelta_t',vals.Kdelta_t, 'b_exp',vals.b_exp, ...
    'i0_deg',i0_deg, 'delta0_deg',delta0_deg, 'denominator',denom, ...
    'camber_deg',theta_deg, 'incidence_deg',incidence_deg, ...
    'deviation_deg',deviation_deg, 'metal_in_deg',metal_in_deg, ...
    'metal_out_deg',metal_out_deg, 'stagger_deg',stagger_deg, 'note',m_note);
end

function v = get_opt(s, field, default_value)
if isfield(s,field) && ~isempty(s.(field))
    v = s.(field);
else
    v = default_value;
end
end
