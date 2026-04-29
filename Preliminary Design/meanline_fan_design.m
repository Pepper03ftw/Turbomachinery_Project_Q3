function out = meanline_fan_design(in)
%MEANLINE_FAN_DESIGN Preliminary meanline design for AE4206 fan project.
%
% OUT = MEANLINE_FAN_DESIGN(IN) computes a one-stage rotor-stator
% preliminary design using the project equations and the digitized Howell
% charts supplied by the user.
%
% IMPORTANT DESIGN LOGIC
% ----------------------
% The selected duty coefficients are interpreted as follows:
%   * psi and R are direct design choices.
%   * phi is either:
%       - an initial guess if solve_mode = 'solve_phi_from_geometry', or
%       - the enforced design phi if solve_mode = 'solve_area_from_phi'.
%
% The solver returns constraint residuals for Howell and Lieblein checks.
% Those residuals are suitable for later use in an optimizer. A residual
% <= 0 means the corresponding constraint is satisfied.

arguments
    in struct
end

in = fill_defaults(in, default_fan_inputs());

% ---------- 0) Ambient / inlet total conditions ----------
air = in.air;
amb = isa_atmosphere(in.altitude_m);
M0  = in.M_flight;

% Convert flight static conditions to total conditions at stage inlet.
T01 = amb.T .* (1 + 0.5*(air.gamma-1)*M0.^2);
p01 = amb.p .* (1 + 0.5*(air.gamma-1)*M0.^2).^(air.gamma/(air.gamma-1));

% Angular speed used throughout the meanline model.
omega = 2*pi*in.N_rpm/60;
phi_initial = in.phi;
psi = in.psi;
R   = in.R;

% ---------- 1) Continuity / geometry consistency ----------
% Two supported modes:
%   solve_phi_from_geometry : project sizing assumptions define the annulus
%                             first; continuity then returns phi_used.
%   solve_area_from_phi     : chosen phi is enforced; continuity then returns
%                             the required annulus area and radii.
a01 = sqrt(air.gamma*air.R*T01);

switch lower(in.design.solve_mode)
    case 'solve_phi_from_geometry'
        % Use project sizing assumptions to infer a first-pass annulus.
        Vx_target = in.Mx_design * a01;
        Utip_from_rel = sqrt(max((in.Mrel_tip_design*a01)^2 - Vx_target^2, 0));
        rt = Utip_from_rel / omega;
        rh = in.hub_to_tip * rt;
        rm = mean_radius_from_rule(rt, rh, in.design.mean_radius_rule);
        A1 = pi*(rt^2 - rh^2);

        % Iterate axial velocity so continuity is satisfied for the chosen
        % annulus area and the resulting inlet static density.
        Vx = Vx_target;
        for k = 1:100
            [T1,p1,rho1] = inlet_static_from_Vx(T01,p01,Vx,air);
            Vx_new = in.mdot/(rho1*A1);
            if abs(Vx_new - Vx) < 1e-8*max(1,Vx)
                Vx = Vx_new;
                break
            end
            Vx = 0.5*Vx + 0.5*Vx_new;
        end
        U = omega*rm;
        phi_used = Vx/U;
        phi_actual = phi_used;
        geometry_mode_note = 'Geometry fixed first; phi solved from continuity.';

    case 'solve_area_from_phi'
        % Enforce the chosen design phi. The annulus area then changes so
        % continuity is satisfied for that phi at the solved inlet density.
        phi_used = phi_initial;

        % First estimate mean radius from project tip-Mach sizing so that a
        % blade speed is available. This keeps the speed level tied to the
        % project sizing basis while allowing frontal area to resize.
        Vx_guess = in.Mx_design * a01;
        U_guess = Vx_guess / max(phi_used, 1e-8);
        rm_guess = U_guess / omega;

        % Use hub-to-tip ratio with chosen mean-radius rule to back out hub
        % and tip radii from the mean radius guess. Then iterate area.
        [rt_guess, rh_guess] = radii_from_mean_and_ratio(rm_guess, in.hub_to_tip, in.design.mean_radius_rule);
        rm = mean_radius_from_rule(rt_guess, rh_guess, in.design.mean_radius_rule);
        U = omega*rm;
        Vx = phi_used * U;

        % With phi fixed, continuity determines the required annulus area.
        % That area is then converted into hub and tip radii while holding
        % the hub-to-tip ratio fixed.
        for k = 1:100
            [T1,p1,rho1] = inlet_static_from_Vx(T01,p01,Vx,air);
            A1 = in.mdot/(rho1*Vx);
            [rt, rh] = radii_from_area_and_ratio(A1, in.hub_to_tip);
            rm_new = mean_radius_from_rule(rt, rh, in.design.mean_radius_rule);
            U_new = omega*rm_new;
            Vx_new = phi_used * U_new;
            if abs(Vx_new - Vx) < 1e-8*max(1,Vx)
                rm = rm_new;
                U = U_new;
                Vx = Vx_new;
                break
            end
            rm = rm_new;
            U = U_new;
            Vx = 0.5*Vx + 0.5*Vx_new;
        end
        [T1,p1,rho1] = inlet_static_from_Vx(T01,p01,Vx,air);
        A1 = pi*(rt^2 - rh^2);
        phi_actual = Vx/U;
        geometry_mode_note = 'phi fixed first; annulus area solved from continuity.';

    otherwise
        error('Unknown design.solve_mode: %s', in.design.solve_mode)
end

% Recompute inlet state one final time with converged values.
[T1,p1,rho1] = inlet_static_from_Vx(T01,p01,Vx,air);
U = omega*rm;
phi = phi_used;

% ---------- 2) Velocity triangles (project equations) ----------
% Once phi_used is known, all subsequent triangles and checks use it.
tan_a1 = (1 - R - psi/2)/phi;
tan_b2 = (psi - 1 + phi*tan_a1)/phi;
tan_b1 = tan_a1 - 1/phi;
tan_a2 = tan_b2 + 1/phi;

alpha1 = atan(tan_a1);
alpha2 = atan(tan_a2);
beta1  = atan(tan_b1);
beta2  = atan(tan_b2);

Vt1 = Vx * tan_a1;
Vt2 = Vx * tan_a2;
Wt1 = Vx * tan_b1;
Wt2 = Vx * tan_b2;

V1 = hypot(Vx, Vt1);
V2 = hypot(Vx, Vt2);
W1 = hypot(Vx, Wt1);
W2 = hypot(Vx, Wt2);

RR = (cos(beta2)^2 * tan(beta2))/phi ...
   - (cos(alpha1)^2 * cos(beta2)^2 * tan(alpha1)*tan(beta2))/phi^2 ...
   + (cos(alpha1)^2 * tan(alpha1))/phi;

% ---------- 3) Euler work and stage total-state rise ----------
Delta_h0 = psi * U^2;
T02 = T01 + Delta_h0/air.cp;
p02 = p01 * in.beta_tt_target;
T02s = T01 * in.beta_tt_target.^((air.gamma-1)/air.gamma);
w_is = air.cp*(T02s - T01);

% ---------- 4) Rotor / stator static states (idealized preliminary model) ----------
T2 = T02 - V2^2/(2*air.cp);
p2 = p02 / (T02/T2)^(air.gamma/(air.gamma-1));
rho2 = p2/(air.R*T2);
M2 = V2/sqrt(air.gamma*air.R*T2);
Mrel1 = W1/sqrt(air.gamma*air.R*T1);
Mrel2 = W2/sqrt(air.gamma*air.R*T2);

% The stator is taken to remove exit swirl to the user-specified alpha3.
alpha3 = deg2rad(in.alpha3_deg);
Vt3 = Vx*tan(alpha3);
V3  = hypot(Vx, Vt3);
T03 = T02;
p03 = p02;
T3 = T03 - V3^2/(2*air.cp);
p3 = p03 / (T03/T3)^(air.gamma/(air.gamma-1));
rho3 = p3/(air.R*T3);
M3 = V3/sqrt(air.gamma*air.R*T3);

% ---------- 5) Reynolds numbers ----------
mu1 = sutherland_mu(T1);
mu2 = sutherland_mu(T2);
Re_rotor  = rho1 * W1 * in.rotor_chord / mu1;
Re_stator = rho2 * V2 * in.stator_chord / mu2;

% ---------- 6) Howell turning-limit checks ----------
howell = read_howell_curves(in.howell_paths);

Dbeta_rotor_deg   = abs(rad2deg(beta1 - beta2));
Dalpha_stator_deg = abs(rad2deg(alpha2 - alpha3));

f_beta_rotor  = howell.f_beta2(abs(rad2deg(beta2)));
f_beta_stator = howell.f_beta2(abs(rad2deg(alpha3)));
phi_Re_rotor  = howell.Phi_Re(Re_rotor / 1e5);
phi_Re_stator = howell.Phi_Re(Re_stator / 1e5);

sc_rotor  = in.rotor_pitch / in.rotor_chord;
sc_stator = in.stator_pitch / in.stator_chord;
psi_sc_rotor  = howell.Psi_sc(sc_rotor);
psi_sc_stator = howell.Psi_sc(sc_stator);

D_rotor_star  = f_beta_rotor  * phi_Re_rotor  * psi_sc_rotor;
D_stator_star = f_beta_stator * phi_Re_stator * psi_sc_stator;

howell_rotor_residual_deg  = Dbeta_rotor_deg   - D_rotor_star;
howell_stator_residual_deg = Dalpha_stator_deg - D_stator_star;
howell_ok_rotor  = howell_rotor_residual_deg  <= 0;
howell_ok_stator = howell_stator_residual_deg <= 0;

% ---------- 7) Diffusion-factor checks ----------
DF_rotor_generic  = (1 - cos(beta1)/cos(beta2)) + cos(beta1)/(2*in.rotor_solidity) * (tan(beta1)-tan(beta2));
DF_stator_generic = (1 - cos(alpha2)/cos(alpha3)) + cos(alpha2)/(2*in.stator_solidity) * (tan(alpha2)-tan(alpha3));

DF_rotor_lieblein  = 1 - W2/W1 + (Wt1 - Wt2)/(2*in.rotor_solidity*W1);
DF_stator_lieblein = 1 - V3/V2 + (Vt2 - Vt3)/(2*in.stator_solidity*V2);

DF_rotor_residual  = DF_rotor_lieblein  - in.design.DF_limit;
DF_stator_residual = DF_stator_lieblein - in.design.DF_limit;
DF_rotor_ok  = DF_rotor_residual  <= 0;
DF_stator_ok = DF_stator_residual <= 0;

% ---------- 8) Literature-based modular loss build-up ----------
span = rt - rh;
loss = estimate_losses(in, struct( ...
    'air',air,'span',span,'rt',rt,'rh',rh,'rm',rm,'omega',omega,'Vx',Vx, ...
    'alpha1',alpha1,'alpha2',alpha2,'alpha3',alpha3, ...
    'beta1',beta1,'beta2',beta2,'Vt1',Vt1,'Vt2',Vt2,'Vt3',Vt3, ...
    'V1',V1,'V2',V2,'V3',V3,'W1',W1,'W2',W2, ...
    'M2',M2,'M3',M3,'Mrel1',Mrel1,'Mrel2',Mrel2,'T1',T1,'T2',T2,'T3',T3, ...
    'rotor_pitch',in.rotor_pitch,'stator_pitch',in.stator_pitch, ...
    'rotor_chord',in.rotor_chord,'stator_chord',in.stator_chord, ...
    'rotor_gap',in.rotor_tip_clearance,'stator_gap',in.stator_hub_clearance));

% ---------- 9) Stage performance ----------
% Losses are combined into a preliminary eta_tt estimate by weighting the
% row loss coefficients with the appropriate velocity scales and comparing
% them to the isentropic stage work.
zetaR = loss.zetaR;
zetaS = loss.zetaS;
eta_tt_est = 1 / (1 + (zetaS*V3^2 + zetaR*W2^2)/(2*w_is));

% Keep target total-pressure ratio as the selected design point. This code
% currently degrades eta_tt through row losses but does not yet close the
% loop to reduce achieved beta_tt.
p03_est = p01 * in.beta_tt_target;
beta_tt_est = p03_est/p01;
Power = in.mdot * Delta_h0;

% ---------- 10) Aerodynamic and centrifugal force estimates ----------
Torque = in.mdot * rm * (Vt2 - Vt1);
Ft_tan = Torque / max(rm,1e-8);
rotor_blades = max(1, round(2*pi*rm / in.rotor_pitch));
stator_vanes = max(1, round(2*pi*rm / in.stator_pitch));
blade_volume = in.rotor_chord * span * in.rotor_thickness * in.blade_height_factor;
blade_mass = in.material.rho_blade * blade_volume;
Fc_blade = blade_mass * omega^2 * rm;

% ---------- 11) Constraint bundle for later optimization ----------
residuals = struct();
residuals.howell_rotor_deg  = howell_rotor_residual_deg;
residuals.howell_stator_deg = howell_stator_residual_deg;
residuals.DF_rotor          = DF_rotor_residual;
residuals.DF_stator         = DF_stator_residual;
residuals.all_ineq = [howell_rotor_residual_deg; howell_stator_residual_deg; DF_rotor_residual; DF_stator_residual];

is_feasible = all(residuals.all_ineq <= 0);
if in.design.enforce_hard_limits && ~is_feasible
    error('Design violated Howell or diffusion-factor feasibility limits.');
end

% ---------- 12) Package outputs ----------
out.input = in;
out.ambient = amb;
out.station = [ ...
    make_station('1',T1,p1,rho1,Vx,Vt1,V1,alpha1,NaN,M0,T01,p01), ...
    make_station('2',T2,p2,rho2,Vx,Vt2,V2,alpha2,beta2,M2,T02,p02), ...
    make_station('3',T3,p3,rho3,Vx,Vt3,V3,alpha3,NaN,M3,T03,p03_est) ...
    ];
out.geometry = struct('rt',rt,'rh',rh,'rm',rm,'A1',A1, ...
    'rotor_pitch',in.rotor_pitch,'stator_pitch',in.stator_pitch, ...
    'rotor_chord',in.rotor_chord,'stator_chord',in.stator_chord, ...
    'rotor_tip_clearance',in.rotor_tip_clearance, ...
    'stator_hub_clearance',in.stator_hub_clearance, ...
    'rotor_blades',rotor_blades,'stator_vanes',stator_vanes, ...
    'span_sections',{{'hub','mean','tip'}},'r_sections',[rh rm rt], ...
    'solve_mode',in.design.solve_mode,'solve_mode_note',geometry_mode_note);

out.coeffs = struct('phi_initial',phi_initial,'phi_actual',phi_actual,'phi_used',phi, ...
    'psi',psi,'R',R,'RR',RR,'beta_tt_target',in.beta_tt_target,'beta_tt_est',beta_tt_est);
out.angles_deg = struct('alpha1',rad2deg(alpha1),'beta1',rad2deg(beta1), ...
    'beta2',rad2deg(beta2),'alpha2',rad2deg(alpha2),'alpha3',rad2deg(alpha3));
out.velocities = struct('Vx',Vx,'Vt1',Vt1,'Vt2',Vt2,'Vt3',Vt3, ...
    'Wt1',Wt1,'Wt2',Wt2,'V1',V1,'V2',V2,'V3',V3,'W1',W1,'W2',W2);
out.reynolds = struct('rotor',Re_rotor,'stator',Re_stator);
out.howell = struct('rotor_required_deg',Dbeta_rotor_deg,'rotor_limit_deg',D_rotor_star, ...
    'stator_required_deg',Dalpha_stator_deg,'stator_limit_deg',D_stator_star, ...
    'rotor_residual_deg',howell_rotor_residual_deg, ...
    'stator_residual_deg',howell_stator_residual_deg, ...
    'rotor_ok',howell_ok_rotor,'stator_ok',howell_ok_stator, ...
    'f_beta_rotor',f_beta_rotor,'f_beta_stator',f_beta_stator, ...
    'phi_re_rotor',phi_Re_rotor,'phi_re_stator',phi_Re_stator, ...
    'psi_sc_rotor',psi_sc_rotor,'psi_sc_stator',psi_sc_stator);
out.diffusion = struct('rotor_generic',DF_rotor_generic,'stator_generic',DF_stator_generic, ...
    'rotor_lieblein',DF_rotor_lieblein,'stator_lieblein',DF_stator_lieblein, ...
    'limit',in.design.DF_limit, ...
    'rotor_residual',DF_rotor_residual,'stator_residual',DF_stator_residual, ...
    'rotor_ok',DF_rotor_ok,'stator_ok',DF_stator_ok);
out.loss = loss;
out.performance = struct('Delta_h0',Delta_h0,'w_is',w_is,'eta_tt_est',eta_tt_est, ...
    'power_W',Power,'torque_Nm',Torque,'tangential_force_N',Ft_tan, ...
    'centrifugal_force_per_blade_N',Fc_blade);
out.constraints = struct('residuals',residuals,'is_feasible',is_feasible);
out.flags = struct('DF_rotor_high',~DF_rotor_ok, 'DF_stator_high',~DF_stator_ok, ...
    'Mrel1_high',Mrel1 > 1.4, 'Mrel2_high',Mrel2 > 1.4, 'M2_high',M2 > 1.0);
end

function [T,p,rho] = inlet_static_from_Vx(T0,p0,Vx,air)
T = T0 - Vx^2/(2*air.cp);
if T <= 1
    error('Inlet static temperature became non-physical.');
end
p = p0 / (T0/T)^(air.gamma/(air.gamma-1));
rho = p/(air.R*T);
end

function rm = mean_radius_from_rule(rt, rh, rule)
switch lower(rule)
    case 'arithmetic'
        rm = 0.5*(rt + rh);
    case 'area'
        rm = sqrt(0.5*(rt^2 + rh^2));
    otherwise
        error('Unknown design.mean_radius_rule: %s', rule)
end
end

function [rt, rh] = radii_from_area_and_ratio(A, lambda)
% Solve A = pi*(rt^2-rh^2) with rh = lambda*rt.
rt = sqrt(A/(pi*(1-lambda^2)));
rh = lambda*rt;
end

function [rt, rh] = radii_from_mean_and_ratio(rm, lambda, rule)
switch lower(rule)
    case 'arithmetic'
        rt = 2*rm/(1+lambda);
    case 'area'
        rt = rm / sqrt(0.5*(1+lambda^2));
    otherwise
        error('Unknown design.mean_radius_rule: %s', rule)
end
rh = lambda*rt;
end

function st = make_station(name,T,p,rho,Vx,Vt,V,alpha,beta,M,T0,p0)
st = struct('name',name,'T',T,'p',p,'rho',rho,'Vx',Vx,'Vt',Vt,'V',V, ...
    'alpha_deg',rad2deg(alpha),'beta_deg',rad2deg(beta), ...
    'Mach',M,'T0',T0,'p0',p0);
end

function out = fill_defaults(in, defaults)
out = defaults;
fn = fieldnames(in);
for i = 1:numel(fn)
    key = fn{i};
    if isstruct(in.(key)) && isfield(defaults,key) && isstruct(defaults.(key))
        out.(key) = fill_defaults(in.(key), defaults.(key));
    else
        out.(key) = in.(key);
    end
end
end

function loss = estimate_losses(in, s)
% loss_model.type is treated as a documentation label only. The active loss
% physics is selected solely through the enabled/disabled source terms in
% default_fan_inputs.m.
rotorRadial = make_radial_row_data('rotor', in.loss_model, s, ...
    s.rotor_pitch, s.rotor_chord, s.T1, s.T2);
statorRadial = make_radial_row_data('stator', in.loss_model, s, ...
    s.stator_pitch, s.stator_chord, s.T2, s.T3);

rotor = row_loss_published('rotor', in.loss_model, s.air, s.span, ...
    abs(s.beta1), abs(s.beta2), s.W1, s.W2, s.Mrel1, s.Mrel2, s.T1, s.T2, ...
    s.rotor_pitch, s.rotor_chord, s.rotor_gap, rotorRadial);
stator = row_loss_published('stator', in.loss_model, s.air, s.span, ...
    abs(s.alpha2), abs(s.alpha3), s.V2, s.V3, s.M2, s.M3, s.T2, s.T3, ...
    s.stator_pitch, s.stator_chord, s.stator_gap, statorRadial);

loss = make_loss_struct(rotor.total, stator.total, in.loss_model.type, ...
    ['Modular literature-based loss model. loss_model.type is a label only; ', ...
     'the included loss sources are controlled exclusively by the enabled ', ...
     'switches in default_fan_inputs.m. Howell/Lieblein are checks only.'], ...
     struct('rotor',rotor,'stator',stator));
loss.radial.rotor = rotorRadial;
loss.radial.stator = statorRadial;
end

function row = row_loss_published(rowType, lm, air, span, chi_in, chi_out, v_in, v_out, Min, Mout, Tin, Tout, pitch, chord, gap, radial)
% Literature-based modular row loss estimate with all switches controlled
% from default_fan_inputs.m.

zBL = 0;
if isfield(lm,'profile') && lm.profile.enabled
    % Denton idealized boundary-layer/profile loss model.
    r = max(lm.profile.dV_over_Vbar, 1e-6);
    Cd = lm.profile.Cd_bl;
    deltaTan = max(abs(tan(chi_in) - tan(chi_out)), 0);
    zBL = Cd * (2/r + 6*r) * deltaTan;
end

zTE = 0;
if isfield(lm,'trailing_edge') && lm.trailing_edge.enabled
    if lm.trailing_edge.use_direct_t_over_s
        ts = lm.trailing_edge.t_over_s;
    else
        ts = chord_to_pitch_times_thickness(chord, pitch, lm);
    end
    if lm.trailing_edge.use_bl_terms
        theta_over_s = lm.trailing_edge.theta_over_s;
        dstar_over_s = lm.trailing_edge.dstar_over_s;
        zTE = 2*theta_over_s + (dstar_over_s + ts)^2 - lm.trailing_edge.Cpb*ts;
    else
        zTE = ts^2 - lm.trailing_edge.Cpb*ts;
    end
    zTE = max(zTE, 0);
end

zSW = 0;
shockDetails = struct('model','none','evaluate_at','none', ...
    'section_names',{{}},'section_weights',[],'M1',[],'M2_cv',[], ...
    'beta1_deg',[],'chi1_deg',[],'incidence_deg',[], ...
    'thickness_over_gap',[],'p02_over_p01_rel',[], ...
    'Yrel',[],'ds_over_R',[],'h_loss',[],'zeta_sections',[], ...
    'velocity_scale',max(v_out,1e-8));
if isfield(lm,'shock') && lm.shock.enabled
    switch lower(field_or(lm.shock,'model','weak_shock_entropy'))
        case {'weak_shock_entropy','weak_normal_shock'}
            switch lower(field_or(lm.shock,'apply_to','outlet'))
                case 'outlet'
                    Mchar = Mout; Tchar = Tout; vref = max(v_out, 1e-6);
                case 'max'
                    if Mout >= Min
                        Mchar = Mout; Tchar = Tout; vref = max(v_out, 1e-6);
                    else
                        Mchar = Min; Tchar = Tin; vref = max(v_in, 1e-6);
                    end
                otherwise
                    error('Unknown shock.apply_to option: %s', lm.shock.apply_to)
            end
            if Mchar > lm.shock.M_crit
                ds_over_R = (2/3) * air.gamma / (air.gamma + 1)^2 * (Mchar^2 - 1)^3;
                Tds = Tchar * air.R * ds_over_R;
                zSW = Tds / (0.5*vref^2);
            else
                ds_over_R = 0;
            end
            shockDetails = struct('model','weak_shock_entropy','evaluate_at',field_or(lm.shock,'apply_to','outlet'), ...
                'section_names',{{'mean'}},'section_weights',1,'M1',Mchar,'M2_cv',NaN, ...
                'beta1_deg',NaN,'chi1_deg',NaN,'incidence_deg',NaN, ...
                'thickness_over_gap',NaN,'p02_over_p01_rel',NaN, ...
                'Yrel',NaN,'ds_over_R',ds_over_R,'h_loss',0.5*zSW*vref^2, ...
                'zeta_sections',zSW,'velocity_scale',vref);

        case {'freeman_cumpsty_inlet_cv','freeman_cumpsty_combined'}
            [zSW, shockDetails] = freeman_cumpsty_shock_loss(rowType, lm.shock, air, radial, v_out);

        otherwise
            error('Unknown shock.model option: %s', lm.shock.model)
    end
end

zTL = 0;
tipDetails = struct('model','none','surface_velocity_mode','none', ...
    'uSS_over_Vx',NaN,'uPS_over_Vx',NaN,'x_over_Cs',NaN, ...
    'stagger_deg',NaN,'solidity',NaN,'D_TL',NaN,'I_tip',NaN, ...
    'Vx',NaN,'velocity_scale',NaN);
if isfield(lm,'tip') && gap > 0
    enableTip = (strcmpi(rowType,'rotor')  && lm.tip.enabled_rotor) || ...
                (strcmpi(rowType,'stator') && lm.tip.enabled_stator);
    if enableTip
        switch lower(field_or(lm.tip,'model','legacy_area_proxy'))
            case 'legacy_area_proxy'
                % Previous one-line leakage-area proxy. Kept for back-to-back
                % comparison only.
                hb = max(span, 1e-6);
                throatFactor = max(cos(chi_in), 0.1);
                Cd_tip = field_or(lm.tip,'Cd_tip',0.002);
                zTL = 2 * Cd_tip * gap * chord / (hb * pitch * throatFactor);
                tipDetails.model = 'legacy_area_proxy';
                tipDetails.velocity_scale = max(v_out,1e-8);

            case 'denton_hall_unshrouded_optiona'
                % Denton/Hall unshrouded tip-clearance loss.
                % Hall Appendix B gives Phi/(mdot*Vx^2); convert it to the
                % code convention zeta such that loss = 0.5*zeta*v_out^2.
                hb = max(span, 1e-6);
                sigma = chord / max(pitch,1e-8);
                Cd_leak = field_or(lm.tip,'Cd_leak',0.8);
                Cs_over_c = field_or(lm.tip,'Cs_over_c',1.0);
                Cs = Cs_over_c * chord;
                Vx_row = 0.5*(abs(v_in*cos(chi_in)) + abs(v_out*cos(chi_out)));
                Vx_row = max(Vx_row, 1e-8);

                [uSS,uPS,x_over_Cs,stagger,modeUsed] = ...
                    tip_surface_velocity_ratios(rowType, lm.tip, chi_in, chi_out, sigma);

                % Guard against tiny numerical/optimizer excursions outside
                % the physical domain. If direct velocities are used, this
                % avoids complex values but does not make bad inputs valid.
                uSS = max(uSS, 0);
                uPS = max(uPS, 0);
                uSS = max(uSS, uPS);

                integrand = uSS .* max(uSS - uPS,0) .* sqrt(max(uSS.^2 - uPS.^2,0));
                if isscalar(integrand)
                    I_tip = integrand;
                else
                    I_tip = trapz(x_over_Cs, integrand);
                end

                D_TL = Cd_leak * (gap * Cs)/(hb * chord) * sigma * I_tip;
                zTL = 2 * D_TL * (Vx_row/max(v_out,1e-8))^2;

                tipDetails = struct('model','denton_hall_unshrouded_optionA', ...
                    'surface_velocity_mode',modeUsed, ...
                    'uSS_over_Vx',uSS,'uPS_over_Vx',uPS,'x_over_Cs',x_over_Cs, ...
                    'stagger_deg',rad2deg(stagger),'solidity',sigma, ...
                    'D_TL',D_TL,'I_tip',I_tip,'Vx',Vx_row, ...
                    'velocity_scale',max(v_out,1e-8));

            otherwise
                error('Unknown tip.model option: %s', lm.tip.model)
        end
    end
end

zEW = 0;
endwallDetails = struct('model','none','surface_velocity_mode','none', ...
    'evaluate_at','none','section_names',{{}},'section_weights',[], ...
    'uSS_over_Vx',[],'uPS_over_Vx',[],'x_over_Cs',[], ...
    'stagger_deg',[],'solidity',[],'Cd',NaN,'n_endwalls',NaN, ...
    'Cs_over_c',NaN,'I_endwall',[],'D_EW',[],'D_EW_total',NaN, ...
    'Vx',[],'velocity_scale',max(v_out,1e-8));
if isfield(lm,'endwall')
    enableEndwall = (strcmpi(rowType,'rotor')  && lm.endwall.enabled_rotor) || ...
                    (strcmpi(rowType,'stator') && lm.endwall.enabled_stator);
    if enableEndwall
        switch lower(field_or(lm.endwall,'model','user_constant'))
            case 'user_constant'
                if strcmpi(rowType,'rotor')
                    zEW = field_or(lm.endwall,'rotor_zeta',0.0);
                else
                    zEW = field_or(lm.endwall,'stator_zeta',0.0);
                end
                endwallDetails.model = 'user_constant';
                endwallDetails.evaluate_at = 'mean_combined';
                endwallDetails.velocity_scale = max(v_out,1e-8);

            case {'hall_denton_cd_baseline','hall_constant_cd','denton_cd'}
                % Hall/Denton constant-C_D endwall boundary-layer
                % dissipation only. This is a baseline endwall-wall
                % friction model; it deliberately excludes secondary-flow,
                % corner-separation, blockage and tip-vortex interaction
                % losses.
                [zEW, endwallDetails] = hall_denton_endwall_loss(rowType, lm, air, span, radial, v_out);

            otherwise
                error('Unknown endwall.model option: %s', lm.endwall.model)
        end
    end
end

row = struct('row_type',rowType, 'profile',zBL, 'trailing_edge',zTE, ...
    'shock',zSW, 'tip',zTL, 'endwall',zEW, ...
    'shock_details',shockDetails, 'tip_details',tipDetails, 'endwall_details',endwallDetails, ...
    'total',zBL + zTE + zSW + zTL + zEW);
end

function [uSS,uPS,x_over_Cs,stagger,modeUsed] = tip_surface_velocity_ratios(rowType, tip, chi_in, chi_out, sigma)
modeUsed = lower(field_or(tip,'surface_velocity_mode','approx'));
switch modeUsed
    case 'approx'
        % Hall Appendix-B option-A approximation:
        %   uSS - uPS ~= Delta(tan chi)/sigma
        %   uSS + uPS ~= 2/cos(stagger)
        % where u = V_surface/Vx.  For the rotor, Delta(tan chi)/sigma is
        % equivalent to psi/(sigma*phi) when chi are row-relative angles.
        stagger = tip_stagger_angle(rowType, tip, chi_in, chi_out);
        uDiff = abs(tan(chi_in) - tan(chi_out)) / max(sigma,1e-8);
        uSum  = 2 / max(cos(stagger),1e-8);
        uSS = 0.5*(uSum + uDiff);
        uPS = 0.5*(uSum - uDiff);
        x_over_Cs = NaN;

    case 'direct'
        % Direct values are expected to be normalized by Vx. Scalars give
        % the algebraic form; vectors give the original chordwise integral.
        [uSS,uPS,x_over_Cs] = tip_direct_surface_velocities(rowType, tip);
        if numel(uSS) ~= numel(uPS)
            error('Direct uSS_over_Vx and uPS_over_Vx must have the same length.')
        end
        if isscalar(uSS)
            x_over_Cs = NaN;
        elseif any(~isfinite(x_over_Cs))
            x_over_Cs = linspace(0,1,numel(uSS));
        end
        stagger = tip_stagger_angle(rowType, tip, chi_in, chi_out);

    otherwise
        error('Unknown tip.surface_velocity_mode option: %s', modeUsed)
end
end

function stagger = tip_stagger_angle(rowType, tip, chi_in, chi_out)
switch lower(field_or(tip,'stagger_mode','mean_flow_angles'))
    case 'mean_flow_angles'
        % Circular-arc / low-camber first pass: stagger is approximated by
        % the mean of the inlet and outlet row-frame flow angle magnitudes.
        stagger = 0.5*(chi_in + chi_out);
    case 'user'
        if strcmpi(rowType,'rotor')
            stagger_deg = field_or(tip,'rotor_stagger_deg',NaN);
        else
            stagger_deg = field_or(tip,'stator_stagger_deg',NaN);
        end
        if ~isfinite(stagger_deg)
            error('tip.stagger_mode = user requires a finite row-specific stagger angle.')
        end
        stagger = deg2rad(stagger_deg);
    otherwise
        error('Unknown tip.stagger_mode option: %s', tip.stagger_mode)
end
end

function [uSS,uPS,x_over_Cs] = tip_direct_surface_velocities(rowType, tip)
if strcmpi(rowType,'rotor')
    uSS = field_or(tip,'rotor_uSS_over_Vx',NaN);
    uPS = field_or(tip,'rotor_uPS_over_Vx',NaN);
    x_over_Cs = field_or(tip,'rotor_x_over_Cs',NaN);
else
    uSS = field_or(tip,'stator_uSS_over_Vx',NaN);
    uPS = field_or(tip,'stator_uPS_over_Vx',NaN);
    x_over_Cs = field_or(tip,'stator_x_over_Cs',NaN);
end
if any(~isfinite(uSS(:))) || any(~isfinite(uPS(:)))
    error('tip.surface_velocity_mode = direct requires finite uSS_over_Vx and uPS_over_Vx values.')
end
uSS = uSS(:).';
uPS = uPS(:).';
x_over_Cs = x_over_Cs(:).';
end


function [uSS,uPS,x_over_Cs,stagger,modeUsed] = endwall_surface_velocity_ratios(rowType, lm, chi_in, chi_out, sigma)
ew = lm.endwall;
modeUsed = lower(field_or(ew,'surface_velocity_mode','approx'));
switch modeUsed
    case 'approx'
        % Hall Appendix-D baseline: endwall velocity is linearly
        % interpolated between pressure-side and suction-side velocities.
        % Without blade-surface data, estimate the two edge velocities from
        % the same low-camber velocity split used in the Denton/Hall tip
        % model:
        %   uSS - uPS ~= Delta(tan chi)/sigma
        %   uSS + uPS ~= 2/cos(stagger)
        stagger = endwall_stagger_angle(rowType, ew, chi_in, chi_out);
        uDiff = abs(tan(chi_in) - tan(chi_out)) / max(sigma,1e-8);
        uSum  = 2 / max(cos(stagger),1e-8);
        uSS = 0.5*(uSum + uDiff);
        uPS = 0.5*(uSum - uDiff);
        x_over_Cs = NaN;

    case 'direct'
        % Direct values are expected to be normalized by Vx. Scalars give
        % the algebraic passage-average form; vectors trigger numerical
        % integration along x/Cs.
        [uSS,uPS,x_over_Cs] = endwall_direct_surface_velocities(rowType, ew);
        if numel(uSS) ~= numel(uPS)
            error('Direct endwall uSS_over_Vx and uPS_over_Vx must have the same length.')
        end
        if isscalar(uSS)
            x_over_Cs = NaN;
        elseif any(~isfinite(x_over_Cs))
            x_over_Cs = linspace(0,1,numel(uSS));
        end
        stagger = endwall_stagger_angle(rowType, ew, chi_in, chi_out);

    case 'tip_shared'
        % Convenience mode: reuse the tip model's surface-velocity settings.
        % This is useful while the same MULTALL/StageN extraction is used
        % for both tip leakage and endwall BL loss.
        if ~isfield(lm,'tip')
            error('endwall.surface_velocity_mode = tip_shared requires lm.tip to exist.')
        end
        [uSS,uPS,x_over_Cs,stagger,tipMode] = ...
            tip_surface_velocity_ratios(rowType, lm.tip, chi_in, chi_out, sigma);
        modeUsed = ['tip_shared_' tipMode];

    otherwise
        error('Unknown endwall.surface_velocity_mode option: %s', modeUsed)
end
end

function I_endwall = endwall_velocity_cubed_integral(uSS,uPS,x_over_Cs)
% Pitchwise mean of u_wall^3 for a linear pressure-side to suction-side
% interpolation:
%   u_wall(y) = uPS + y*(uSS-uPS), 0 <= y <= 1.
% The exact pitchwise average is:
%   integral_0^1 u_wall^3 dy = (uSS^4-uPS^4)/(4*(uSS-uPS)).
du = uSS - uPS;
mean_y = zeros(size(uSS));
idx = abs(du) > 1e-12;
mean_y(idx) = (uSS(idx).^4 - uPS(idx).^4) ./ (4*du(idx));
mean_y(~idx) = uSS(~idx).^3;

if isscalar(mean_y)
    I_endwall = mean_y;
else
    I_endwall = trapz(x_over_Cs, mean_y);
end
I_endwall = max(I_endwall, 0);
end

function stagger = endwall_stagger_angle(rowType, ew, chi_in, chi_out)
switch lower(field_or(ew,'stagger_mode','mean_flow_angles'))
    case 'mean_flow_angles'
        % Low-order first pass: use the mean of inlet and outlet row-frame
        % flow-angle magnitudes as the representative blade stagger.
        stagger = 0.5*(chi_in + chi_out);
    case 'user'
        if strcmpi(rowType,'rotor')
            stagger_deg = field_or(ew,'rotor_stagger_deg',NaN);
        else
            stagger_deg = field_or(ew,'stator_stagger_deg',NaN);
        end
        if ~isfinite(stagger_deg)
            error('endwall.stagger_mode = user requires a finite row-specific stagger angle.')
        end
        stagger = deg2rad(stagger_deg);
    otherwise
        error('Unknown endwall.stagger_mode option: %s', ew.stagger_mode)
end
end

function [uSS,uPS,x_over_Cs] = endwall_direct_surface_velocities(rowType, ew)
if strcmpi(rowType,'rotor')
    uSS = field_or(ew,'rotor_uSS_over_Vx',NaN);
    uPS = field_or(ew,'rotor_uPS_over_Vx',NaN);
    x_over_Cs = field_or(ew,'rotor_x_over_Cs',NaN);
else
    uSS = field_or(ew,'stator_uSS_over_Vx',NaN);
    uPS = field_or(ew,'stator_uPS_over_Vx',NaN);
    x_over_Cs = field_or(ew,'stator_x_over_Cs',NaN);
end
if any(~isfinite(uSS(:))) || any(~isfinite(uPS(:)))
    error('endwall.surface_velocity_mode = direct requires finite uSS_over_Vx and uPS_over_Vx values.')
end
uSS = uSS(:).';
uPS = uPS(:).';
x_over_Cs = x_over_Cs(:).';
end


function radial = make_radial_row_data(rowType, lm, s, pitch_mean, chord, Tin, Tout)
% Build hub/mean/tip row data from the current meanline triangles.
% This is intentionally low-order. Once MEANGEN/STAGEN/MULTALL data are
% available, these fields can be replaced by section-resolved profiles while
% keeping the loss-model interfaces unchanged.
r = [s.rh, s.rm, s.rt];
names = {'hub','mean','tip'};
Vx = s.Vx * ones(size(r));
VthetaMode = 'constant_vtheta';
if isfield(lm,'radial')
    VthetaMode = field_or(lm.radial,'vtheta_mode','constant_vtheta');
end
pitchMode = 'scale_with_radius';
if isfield(lm,'radial')
    pitchMode = field_or(lm.radial,'pitch_mode','scale_with_radius');
end

switch lower(rowType)
    case 'rotor'
        Vt_in  = spanwise_vtheta(s.Vt1, r, s.rm, VthetaMode);
        Vt_out = spanwise_vtheta(s.Vt2, r, s.rm, VthetaMode);
        U = s.omega .* r;
        Wt_in = Vt_in - U;
        Wt_out = Vt_out - U;
        v_in = hypot(Vx, Wt_in);
        v_out = hypot(Vx, Wt_out);
        chi_in_signed = atan2(Wt_in, Vx);
        chi_out_signed = atan2(Wt_out, Vx);

    case 'stator'
        Vt_in  = spanwise_vtheta(s.Vt2, r, s.rm, VthetaMode);
        Vt_out = spanwise_vtheta(s.Vt3, r, s.rm, VthetaMode);
        U = zeros(size(r));
        v_in = hypot(Vx, Vt_in);
        v_out = hypot(Vx, Vt_out);
        chi_in_signed = atan2(Vt_in, Vx);
        chi_out_signed = atan2(Vt_out, Vx);

    otherwise
        error('Unknown rowType: %s', rowType)
end

switch lower(pitchMode)
    case 'scale_with_radius'
        pitch = pitch_mean .* r ./ max(s.rm,1e-8);
    case 'constant'
        pitch = pitch_mean .* ones(size(r));
    otherwise
        error('Unknown loss_model.radial.pitch_mode: %s', pitchMode)
end

Min = v_in ./ sqrt(s.air.gamma*s.air.R*Tin);
Mout = v_out ./ sqrt(s.air.gamma*s.air.R*Tout);
radial = struct('row_type',rowType,'names',{names},'r',r,'Vx',Vx,'U',U, ...
    'pitch',pitch,'chord',chord.*ones(size(r)), ...
    'chi_in_signed',chi_in_signed,'chi_out_signed',chi_out_signed, ...
    'chi_in',abs(chi_in_signed),'chi_out',abs(chi_out_signed), ...
    'v_in',v_in,'v_out',v_out,'Min',Min,'Mout',Mout, ...
    'Tin',Tin.*ones(size(r)),'Tout',Tout.*ones(size(r)));
end

function Vt = spanwise_vtheta(Vt_mean, r, rm, mode)
switch lower(mode)
    case {'constant_vtheta','constant'}
        Vt = Vt_mean .* ones(size(r));
    case {'free_vortex','constant_rvtheta'}
        Vt = Vt_mean .* rm ./ max(r,1e-8);
    otherwise
        error('Unknown spanwise vtheta mode: %s', mode)
end
end

function [idx,weights,mode] = select_section_indices(settings, rowType, defaultMode)
if strcmpi(rowType,'rotor')
    mode = field_or(settings,'rotor_evaluate_at',field_or(settings,'evaluate_at',defaultMode));
else
    mode = field_or(settings,'stator_evaluate_at',field_or(settings,'evaluate_at',defaultMode));
end
mode = lower(mode);
switch mode
    case 'hub'
        idx = 1; weights = 1;
    case {'mean','meanline'}
        idx = 2; weights = 1;
    case {'tip','casing'}
        idx = 3; weights = 1;
    case {'hub_tip','endwalls'}
        idx = [1 3]; weights = [0.5 0.5];
    case {'hub_mean_tip','all','spanwise'}
        idx = [1 2 3]; weights = [0.25 0.5 0.25];
    case {'mean_combined','legacy_mean'}
        idx = 2; weights = 1;
    otherwise
        error('Unknown evaluate_at option: %s', mode)
end
end

function [zSW, details] = freeman_cumpsty_shock_loss(rowType, shock, air, radial, v_out_mean)
% Freeman & Cumpsty inlet-region control-volume shock model.
% Solves their combined mass/relative-stagnation-enthalpy/momentum equation
% for the subsonic post-inlet-region Mach number, then converts the relative
% stagnation-pressure loss into the code's zeta convention.
[idx,weights,mode] = select_section_indices(shock, rowType, 'tip');
if strcmpi(rowType,'stator') && ~field_or(shock,'enabled_stator',false)
    zSW = 0;
    details = empty_freeman_details('freeman_cumpsty_inlet_cv',mode,max(v_out_mean,1e-8));
    return
end
if strcmpi(rowType,'rotor') && ~field_or(shock,'enabled_rotor',true)
    zSW = 0;
    details = empty_freeman_details('freeman_cumpsty_inlet_cv',mode,max(v_out_mean,1e-8));
    return
end

n = numel(idx);
hLoss = zeros(1,n); zeta_i = zeros(1,n); M1v = zeros(1,n); M2v = NaN(1,n);
betaDeg = zeros(1,n); chiDeg = zeros(1,n); incDeg = zeros(1,n);
tauVec = zeros(1,n); p02p01 = ones(1,n); Yrel = zeros(1,n); dsR = zeros(1,n);
sectionNames = cell(1,n);

for k = 1:n
    j = idx(k);
    sectionNames{k} = radial.names{j};
    M1 = radial.Min(j);
    beta1 = radial.chi_in(j);
    chi1 = shock_blade_inlet_angle(rowType, shock, radial, j);
    tau = shock_thickness_over_gap(rowType, shock, j);
    tau = min(max(tau,0),0.95);

    M1v(k) = M1;
    betaDeg(k) = rad2deg(beta1);
    chiDeg(k) = rad2deg(chi1);
    incDeg(k) = rad2deg(beta1 - chi1);
    tauVec(k) = tau;

    if M1 <= field_or(shock,'M_crit',1.0)
        continue
    end

    [M2cv, ok] = solve_freeman_cumpsty_M2(M1, beta1, chi1, tau, air.gamma);
    if ~ok
        if strcmpi(field_or(shock,'fallback','normal_shock'),'normal_shock')
            M2cv = normal_shock_downstream_mach(M1, air.gamma);
        else
            warning('Freeman-Cumpsty shock solve failed at %s %s section; setting shock loss to zero.', rowType, sectionNames{k});
            continue
        end
    end
    M2v(k) = M2cv;

    A1 = 1 + 0.5*(air.gamma-1)*M1^2;
    A2 = 1 + 0.5*(air.gamma-1)*M2cv^2;
    p01_p1 = A1^(air.gamma/(air.gamma-1));

    % Static pressure ratio from the mass + stagnation-enthalpy equations.
    p2_p1 = (M1/M2cv) * sqrt(A1/A2) * cos(beta1)/(max(cos(chi1),1e-8)*max(1-tau,1e-8));
    p02_p2 = A2^(air.gamma/(air.gamma-1));
    p02_p01 = p2_p1 * p02_p2 / p01_p1;
    p02_p01 = min(max(p02_p01,eps),1.0);

    p02p01(k) = p02_p01;
    Yrel(k) = p01_p1*(1 - p02_p01) / max(p01_p1 - 1,1e-8);
    dsR(k) = -log(p02_p01);

    T2cv = radial.Tin(j) * A1/A2;
    hLoss(k) = T2cv * air.R * dsR(k);
    zeta_i(k) = 2*hLoss(k)/max(v_out_mean,1e-8)^2;
end

hLossWeighted = sum(weights .* hLoss);
zSW = 2*hLossWeighted/max(v_out_mean,1e-8)^2;
details = struct('model','freeman_cumpsty_inlet_cv','evaluate_at',mode, ...
    'section_names',{sectionNames},'section_weights',weights, ...
    'M1',M1v,'M2_cv',M2v,'beta1_deg',betaDeg,'chi1_deg',chiDeg, ...
    'incidence_deg',incDeg,'thickness_over_gap',tauVec, ...
    'p02_over_p01_rel',p02p01,'Yrel',Yrel,'ds_over_R',dsR, ...
    'h_loss',hLoss,'zeta_sections',zeta_i,'velocity_scale',max(v_out_mean,1e-8));
end

function details = empty_freeman_details(model,mode,vscale)
details = struct('model',model,'evaluate_at',mode, ...
    'section_names',{{}},'section_weights',[],'M1',[],'M2_cv',[], ...
    'beta1_deg',[],'chi1_deg',[],'incidence_deg',[], ...
    'thickness_over_gap',[],'p02_over_p01_rel',[], ...
    'Yrel',[],'ds_over_R',[],'h_loss',[],'zeta_sections',[], ...
    'velocity_scale',vscale);
end

function [M2, ok] = solve_freeman_cumpsty_M2(M1, beta1, chi1, tau, gamma)
ok = false;
M2 = NaN;
if M1 <= 1 || tau >= 1
    return
end
f = @(m) freeman_cumpsty_residual(m, M1, beta1, chi1, tau, gamma);
lo = 0.05; hi = 0.999;
N = 300;
grid = linspace(lo,hi,N);
vals = arrayfun(f, grid);
valid = isfinite(vals);
grid = grid(valid); vals = vals(valid);
if numel(grid) < 2
    return
end
sgn = vals(1:end-1).*vals(2:end);
br = find(sgn <= 0, 1, 'first');
if isempty(br)
    % Choose the minimum residual on the subsonic branch if no clean sign
    % change is found, but mark as failed unless the residual is very small.
    [minAbs,ii] = min(abs(vals));
    if minAbs < 1e-7
        M2 = grid(ii); ok = true;
    end
    return
end
try
    M2 = fzero(f, [grid(br) grid(br+1)]);
    ok = isfinite(M2) && M2 > 0 && M2 < 1;
catch
    ok = false;
    M2 = NaN;
end
end

function r = freeman_cumpsty_residual(M2, M1, beta1, chi1, tau, gamma)
A1 = 1 + 0.5*(gamma-1)*M1^2;
A2 = 1 + 0.5*(gamma-1)*M2.^2;
lhs = sqrt(A1) .* (1 + gamma*M2.^2*(1-tau)) ./ ...
    (max(M2,1e-12).*(1-tau).*sqrt(A2));
rhs = (cos(chi1)/max(cos(beta1),1e-8) + gamma*M1^2*cos(beta1-chi1)) / M1;
r = lhs - rhs;
end

function M2 = normal_shock_downstream_mach(M1, gamma)
M2 = sqrt((1 + 0.5*(gamma-1)*M1^2)/(gamma*M1^2 - 0.5*(gamma-1)));
end

function chi1 = shock_blade_inlet_angle(rowType, shock, radial, idx)
mode = lower(field_or(shock,'blade_angle_mode','zero_incidence'));
suffix = radial.names{idx};
switch mode
    case 'zero_incidence'
        if strcmpi(rowType,'rotor')
            inc = field_or(shock,'rotor_design_incidence_deg',field_or(shock,'design_incidence_deg',0));
        else
            inc = field_or(shock,'stator_design_incidence_deg',field_or(shock,'design_incidence_deg',0));
        end
        chi1 = radial.chi_in(idx) - deg2rad(inc);

    case 'user'
        fieldName = [lower(rowType) '_chi1_deg'];
        chiDeg = field_or(shock,fieldName,NaN);
        if ~isfinite(chiDeg)
            error('shock.blade_angle_mode = user requires %s.', fieldName)
        end
        chi1 = deg2rad(chiDeg);

    case {'user_sections','sectional'}
        fieldName = [lower(rowType) '_chi1_' suffix '_deg'];
        chiDeg = field_or(shock,fieldName,NaN);
        if ~isfinite(chiDeg)
            % Fallback to row-wise user value if available.
            chiDeg = field_or(shock,[lower(rowType) '_chi1_deg'],NaN);
        end
        if ~isfinite(chiDeg)
            error('shock.blade_angle_mode = user_sections requires %s or row-wise chi1.', fieldName)
        end
        chi1 = deg2rad(chiDeg);

    otherwise
        error('Unknown shock.blade_angle_mode option: %s', mode)
end
chi1 = min(max(abs(chi1),deg2rad(1e-6)),deg2rad(89.9));
end

function tau = shock_thickness_over_gap(rowType, shock, idx)
names = {'hub','mean','tip'};
suffix = names{idx};
if strcmpi(rowType,'rotor')
    tau = field_or(shock,['rotor_t_over_s_' suffix],NaN);
    if ~isfinite(tau)
        tau = field_or(shock,'rotor_t_over_s',field_or(shock,'t_over_s',0.075));
    end
else
    tau = field_or(shock,['stator_t_over_s_' suffix],NaN);
    if ~isfinite(tau)
        tau = field_or(shock,'stator_t_over_s',field_or(shock,'t_over_s',0.075));
    end
end
end

function [zEW, details] = hall_denton_endwall_loss(rowType, lm, air, span, radial, v_out_mean)
ew = lm.endwall;
[idx,weights,mode] = select_section_indices(ew, rowType, 'hub_tip');
if strcmpi(mode,'mean_combined') || strcmpi(mode,'legacy_mean')
    nEndwallsDefault = 2;
else
    nEndwallsDefault = 1;
end

n = numel(idx);
Cd_ew = field_or(ew,'Cd',0.002);
Cs_over_c = field_or(ew,'Cs_over_c',1.0);
sectionNames = cell(1,n);
D = zeros(1,n); I = zeros(1,n); VxVals = zeros(1,n); sigVals = zeros(1,n);
staggerDeg = zeros(1,n); nEnd = zeros(1,n);
uSS_store = cell(1,n); uPS_store = cell(1,n); x_store = cell(1,n);
for k = 1:n
    j = idx(k);
    sectionNames{k} = radial.names{j};
    pitch_j = radial.pitch(j);
    chord_j = radial.chord(j);
    sigma = chord_j/max(pitch_j,1e-8);
    Cs = Cs_over_c * chord_j;
    Vx_row = max(abs(radial.Vx(j)),1e-8);

    [uSS,uPS,x_over_Cs,stagger,modeEW] = ...
        endwall_surface_velocity_ratios(rowType, lm, radial.chi_in(j), radial.chi_out(j), sigma);
    uSS = max(uSS, 0);
    uPS = max(uPS, 0);
    I_endwall = endwall_velocity_cubed_integral(uSS, uPS, x_over_Cs);

    if strcmpi(mode,'mean_combined') || strcmpi(mode,'legacy_mean')
        if strcmpi(rowType,'rotor')
            nLocal = field_or(ew,'rotor_n_endwalls',field_or(ew,'n_endwalls',2));
        else
            nLocal = field_or(ew,'stator_n_endwalls',field_or(ew,'n_endwalls',2));
        end
    else
        nLocal = nEndwallsDefault;
    end

    D(k) = Cd_ew * nLocal * (Cs/max(span,1e-8)) * Vx_row^2 * I_endwall;
    I(k) = I_endwall;
    VxVals(k) = Vx_row;
    sigVals(k) = sigma;
    staggerDeg(k) = rad2deg(stagger);
    nEnd(k) = nLocal;
    uSS_store{k} = uSS;
    uPS_store{k} = uPS;
    x_store{k} = x_over_Cs;
end

% If both actual endwalls are selected, sum the hub and casing losses. If a
% spanwise diagnostic mode is selected, use the provided integration weights.
if strcmpi(mode,'hub_tip') || strcmpi(mode,'endwalls')
    Dtotal = sum(D);
else
    Dtotal = sum(weights .* D);
end
zEW = 2*Dtotal/max(v_out_mean,1e-8)^2;
details = struct('model','hall_denton_cd_baseline', ...
    'surface_velocity_mode',field_or(ew,'surface_velocity_mode','approx'), ...
    'evaluate_at',mode,'section_names',{sectionNames},'section_weights',weights, ...
    'uSS_over_Vx',{uSS_store},'uPS_over_Vx',{uPS_store},'x_over_Cs',{x_store}, ...
    'stagger_deg',staggerDeg,'solidity',sigVals, ...
    'Cd',Cd_ew,'n_endwalls',nEnd,'Cs_over_c',Cs_over_c, ...
    'I_endwall',I,'D_EW',D,'D_EW_total',Dtotal, ...
    'Vx',VxVals,'velocity_scale',max(v_out_mean,1e-8));
end

function val = field_or(s, name, defaultVal)
if isfield(s,name)
    val = s.(name);
else
    val = defaultVal;
end
end

function out = make_loss_struct(zR, zS, typeLabel, note, breakdown)
out = struct('type_label',typeLabel, 'note',note, 'zetaR',zR, 'zetaS',zS, ...
    'breakdown',breakdown);
end

function ts = chord_to_pitch_times_thickness(chord, pitch, lm)
% Fallback helper if the user later wants to switch from direct t/s input
% to a geometry-based estimate. No geometry-based TE thickness model is
% hard-coded yet, so keep this tied to the user-supplied first-pass value.
ts = lm.trailing_edge.t_over_s * chord/chord * pitch/pitch;
end

function mu = sutherland_mu(T)
mu = 1.716e-5 * (T/273.15).^(3/2) .* (273.15 + 110.4) ./ (T + 110.4);
end

function amb = isa_atmosphere(h)
% International Standard Atmosphere, troposphere approximation.
T0 = 288.15; p0 = 101325; L = 0.0065; g = 9.80665; R = 287.05;
T = T0 - L*h;
p = p0 * (T/T0)^(g/(R*L));
rho = p/(R*T);
amb = struct('T',T,'p',p,'rho',rho,'a',sqrt(1.4*R*T));
end
