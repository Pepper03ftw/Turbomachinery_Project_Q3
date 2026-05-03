function out = meanline_fan_design(in)
%MEANLINE_FAN_DESIGN Preliminary meanline design for AE4206 fan project.
%
% This file is intentionally split into two layers:
%   1) a fast meanline candidate evaluator; and
%   2) an optional duty-coefficient sweep that searches for good MEANGEN
%      starting values before the slower MEANGEN -> STAGEN -> MULTALL loop.
%
% The modular loss models below this meanline driver are not used as design
% variables. They are only called to score a candidate after the meanline
% geometry, velocity triangles, Howell check and diffusion-factor check have
% been established.

arguments
    in struct
end

in = fill_defaults(in, default_fan_inputs());

% Resolve local data files once at the top level so the code does not depend
% on the MATLAB current working directory.
in = resolve_local_data_paths(in);

status_msg(in,'START','meanline_fan_design entered');

% Recommended mode for this project chat: sweep phi, psi, and R inside the
% hard course limits while holding the project operating point as constraints.
if (strcmpi(in.design.solve_mode,'sweep_project_constraints') || ...
        strcmpi(in.design.solve_mode,'sweep_fixed_mean_radius_duty')) && ...
        ~field_or(in.design,'sweep_active',false)
    out = meanline_design_sweep_project(in);
    return
end

% Legacy searches remain available for back-comparison, but the new project
% sweep above is the default route.
if isfield(in,'design') && isfield(in.design,'beta_tt_mode') && ...
        strcmpi(in.design.beta_tt_mode,'auto_within_range') && ...
        ~strcmpi(in.design.solve_mode,'project_fixed_candidate') && ...
        ~strcmpi(in.design.solve_mode,'sweep_project_constraints') && ...
        ~strcmpi(in.design.solve_mode,'sweep_fixed_mean_radius_duty')

    if ~isfield(in.design,'beta_tt_search_active') || ~in.design.beta_tt_search_active
        out = meanline_design_auto_beta_tt(in);
        return
    end
end

if isfield(in,'design') && isfield(in.design,'R_mode') && ...
        strcmpi(in.design.R_mode,'auto_balance_DF_Howell') && ...
        ~strcmpi(in.design.solve_mode,'project_fixed_candidate') && ...
        ~strcmpi(in.design.solve_mode,'sweep_project_constraints') && ...
        ~strcmpi(in.design.solve_mode,'sweep_fixed_mean_radius_duty')

    if ~isfield(in.design,'R_search_active') || ~in.design.R_search_active
        out = meanline_design_auto_R(in);
        return
    end
end

% ---------- 0) Ambient / inlet total conditions ----------
status_msg(in,'0/12','ambient and inlet total conditions');

air = in.air;
amb = isa_atmosphere(in.altitude_m);
M0  = in.M_flight;

[T01,p01,inlet_total_note] = inlet_total_conditions(amb,M0,air, ...
    field_or(in.design,'inlet_total_mode','freestream_total'));

omega = 2*pi*in.N_rpm/60;
phi_initial = in.phi;
psi = in.psi;
R   = in.R;

% ---------- 1) Continuity / geometry consistency ----------
status_msg(in,'1/12','geometry and continuity');

a01 = sqrt(air.gamma*air.R*T01);
T02s_target = T01 * in.beta_tt_target.^((air.gamma-1)/air.gamma);
w_is_target = air.cp*(T02s_target - T01);
fixed_radius_duty_active = false;

switch lower(in.design.solve_mode)
    case {'project_fixed_candidate','sweep_project_constraints','sweep_fixed_mean_radius_duty'}
        % Candidate-evaluation mode for the project workflow. Two geometrical
        % closures are supported here:
        %   1) continuity-first: Mx, mdot and hub-to-tip define r_m and phi is computed;
        %   2) fixed-mean-radius duty sweep: r_m is supplied/chosen and phi is swept.
        % The second option mirrors the MEANGEN input philosophy where mean
        % radius/diameter and the duty coefficients are the main preliminary
        % quantities to pass downstream.
        fixed_radius_duty_active = field_or(in.design,'fixed_mean_radius_duty_active',false) || ...
            strcmpi(field_or(in.design,'solve_mode',''),'sweep_fixed_mean_radius_duty');

        % Continuity reference at the project axial Mach.  This is always
        % computed as a diagnostic, even when a separate fixed radius is used.
        Mx_target = in.Mx_design;
        [T1_ref,p1_ref,rho1_ref,Vx_ref] = inlet_static_from_Mx(T01,p01,Mx_target,air); %#ok<ASGLU>
        A1_continuity = in.mdot/(rho1_ref*Vx_ref);
        [rt_continuity,rh_continuity] = radii_from_area_and_ratio(A1_continuity, in.hub_to_tip);
        rm_continuity = mean_radius_from_rule(rt_continuity, rh_continuity, in.design.mean_radius_rule);

        tipmach_radius_match = struct('enabled',false, ...
            'status','disabled', ...
            'rm_continuity',rm_continuity, ...
            'rt_continuity',rt_continuity, ...
            'rh_continuity',rh_continuity, ...
            'A1_continuity',A1_continuity, ...
            'rm_selected',rm_continuity, ...
            'relative_shift',0, ...
            'Mrel_tip_target',in.Mrel_tip_design, ...
            'Mrel_tip_selected',NaN);

        if fixed_radius_duty_active
            rm_user = field_or(in.design,'fixed_mean_radius_m',NaN);
            if isfinite(rm_user) && rm_user > 0
                rm = rm_user;
                rm_source = 'user_fixed_mean_radius_m';
            else
                scale = field_or(in.design,'fixed_mean_radius_scale',1.0);
                rm = scale*rm_continuity;
                rm_source = sprintf('continuity_reference_scaled_by_%.4g',scale);
            end
            [rt,rh] = radii_from_mean_and_ratio(rm, in.hub_to_tip, in.design.mean_radius_rule);
            A1 = pi*(rt^2 - rh^2);

            phi_used = phi_initial;
            Vx = phi_used * omega * rm;
            [T1,p1,rho1] = inlet_static_from_Vx(T01,p01,Vx,air);
            phi_actual = phi_used;

            tipmach_radius_match.status = 'not_used_fixed_mean_radius_duty_sweep';
            tipmach_radius_match.rm_selected = rm;
            tipmach_radius_match.relative_shift = (rm-rm_continuity)/max(rm_continuity,1e-12);

            geometry_mode_note = ['Fixed-mean-radius duty sweep: r_m is fixed from ',rm_source, ...
                ', hub-to-tip sets rt/rh, and phi is swept by setting Vx = phi*Omega*r_m. ', ...
                'Continuity mdot is reported as a diagnostic rather than used to resize the annulus.'];
        else
            % Continuity-first project closure.  The project point fixes Mx,
            % mdot and hub-to-tip ratio, so continuity fixes the annulus area
            % and therefore the mean radius.  Phi is an output/check.
            T1 = T1_ref; p1 = p1_ref; rho1 = rho1_ref; Vx = Vx_ref;
            A1 = A1_continuity;
            rt = rt_continuity; rh = rh_continuity; rm = rm_continuity;

            if field_or(in.design,'match_Mrel_tip_by_mean_radius',false)
                alpha1_mode_for_match = field_or(in.design,'alpha1_mode_active',field_or(in.design,'alpha1_mode','free'));
                [rt,rh,rm,A1,tipmach_radius_match] = match_radius_to_tip_mach( ...
                    rm_continuity, in.hub_to_tip, in.design.mean_radius_rule, ...
                    Vx, T1, air, omega, psi, R, alpha1_mode_for_match, ...
                    field_or(in.design,'alpha1_deg',0), in.Mrel_tip_design, ...
                    in.loss_model.radial.vtheta_mode, in.design);
            end

            U = omega*rm;
            phi_used = Vx/max(U,1e-12);
            phi_actual = phi_used;

            geometry_mode_note = ['Continuity-fixed project candidate: Mx sets Vx and density, ', ...
                'mdot gives the reference annulus/r_m, and hub-to-tip ratio sets rt/rh. ', ...
                'Phi is computed, not swept.'];
            if field_or(in.design,'match_Mrel_tip_by_mean_radius',false)
                geometry_mode_note = [geometry_mode_note, ' Mean radius is then optionally rescaled to match Mrel_tip.'];
            end
        end

    case 'solve_phi_from_geometry'
        Vx_target = in.Mx_design * a01;
        Utip_from_rel = sqrt(max((in.Mrel_tip_design*a01)^2 - Vx_target^2, 0));
        rt = Utip_from_rel / omega;
        rh = in.hub_to_tip * rt;
        rm = mean_radius_from_rule(rt, rh, in.design.mean_radius_rule);
        A1 = pi*(rt^2 - rh^2);

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
        geometry_mode_note = 'Legacy mode: geometry fixed first; phi solved from continuity.';

    case 'solve_area_from_phi'
        phi_used = phi_initial;
        Vx_guess = in.Mx_design * a01;
        U_guess = Vx_guess / max(phi_used, 1e-8);
        rm_guess = U_guess / omega;

        [rt_guess, rh_guess] = radii_from_mean_and_ratio(rm_guess, in.hub_to_tip, in.design.mean_radius_rule);
        rm = mean_radius_from_rule(rt_guess, rh_guess, in.design.mean_radius_rule);
        U = omega*rm;
        Vx = phi_used * U;

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
        geometry_mode_note = 'Legacy mode: phi fixed first; annulus area solved from continuity.';

    case 'solve_phi_psi_from_mach_pr'
        % Retained for old runs. This mode deliberately computes phi and psi
        % from the Mach/pressure-ratio closure and should not be used for
        % the new duty-coefficient sweep.
        Mx_target = in.Mx_design;
        [T1,p1,rho1,Vx] = inlet_static_from_Mx(T01,p01,Mx_target,air);

        eta_work_guess = field_or(in.design,'eta_work_guess',0.90);
        eta_work_guess = min(max(eta_work_guess,0.50),0.999);
        Delta_h0_design = w_is_target/eta_work_guess;

        if field_or(in.design,'enforce_mdot_by_hub_radius',true)
            A_target = in.mdot/(rho1*Vx);
            [rt,rh,rm,phi_used,psi,Mrel_tip_solved,Vt1_tip_solved] = ...
                solve_tip_radius_from_mdot_tipmach(A_target, Vx, T1, air, omega, ...
                R, Delta_h0_design, in.Mrel_tip_design, in.design.mean_radius_rule, ...
                in.loss_model.radial.vtheta_mode); %#ok<ASGLU>

            A1 = pi*(rt^2-rh^2);
            U = omega*rm;
            phi_actual = phi_used;
            geometry_mode_note = ['Legacy Mach/PR mode: mdot fixes area and ', ...
                'tip radius is solved from Mrel_tip. Hub-to-tip ratio is an output.'];
        else
            lambda = in.hub_to_tip;
            f_rm_rt = mean_radius_fraction(lambda, in.design.mean_radius_rule);
            a1 = sqrt(air.gamma*air.R*T1);
            Wtip_target = in.Mrel_tip_design*a1;
            if Wtip_target <= abs(Vx)
                error(['Mrel_tip_design is too low for the requested axial Mach. ', ...
                       'Need Mrel_tip_design > Mx_design approximately.']);
            end
            Krel = sqrt(max(Wtip_target^2 - Vx^2,0));
            Aq = 1 - f_rm_rt*(1 - R);
            Bq = -Krel;
            Cq = Delta_h0_design/(2*max(f_rm_rt,1e-12));
            disc = Bq^2 - 4*Aq*Cq;
            if disc < 0
                error(['No real U_tip satisfies the requested Mx, Mrel_tip, PR, R and eta_work_guess. ', ...
                       'Try a higher eta_work_guess, lower PR, lower R, or relax the Mach constraints.']);
            end
            roots_Utip = [(-Bq + sqrt(disc))/(2*Aq), (-Bq - sqrt(disc))/(2*Aq)];
            roots_Utip = roots_Utip(isfinite(roots_Utip) & roots_Utip > 0);
            if isempty(roots_Utip)
                error('No positive U_tip root found for Mach/pressure-ratio constrained sizing.');
            end
            Utip = max(roots_Utip);

            rt = Utip/omega;
            rh = lambda*rt;
            rm = mean_radius_from_rule(rt, rh, in.design.mean_radius_rule);
            A1 = pi*(rt^2 - rh^2);
            U = omega*rm;

            phi_used = Vx/U;
            phi_actual = phi_used;
            psi = Delta_h0_design/U^2;
            geometry_mode_note = 'Legacy Mach/PR mode with fixed input hub-to-tip ratio.';
        end

    otherwise
        error('Unknown design.solve_mode: %s', in.design.solve_mode)
end

if ~exist('rm_continuity','var')
    rm_continuity = rm;
    rt_continuity = rt;
    rh_continuity = rh;
    A1_continuity = A1;
end
if ~exist('tipmach_radius_match','var')
    tipmach_radius_match = struct('enabled',false,'status','not_applicable', ...
        'rm_continuity',rm_continuity,'rt_continuity',rt_continuity, ...
        'rh_continuity',rh_continuity,'A1_continuity',A1_continuity, ...
        'rm_selected',rm,'relative_shift',(rm-rm_continuity)/max(rm_continuity,1e-12), ...
        'Mrel_tip_target',in.Mrel_tip_design,'Mrel_tip_selected',NaN);
end

[T1,p1,rho1] = inlet_static_from_Vx(T01,p01,Vx,air);
U = omega*rm;
phi = phi_used;

mdot_meanline = rho1*Vx*A1;
Mx_actual = Vx/sqrt(air.gamma*air.R*T1);
Utip_actual = omega*rt;
hub_to_tip_actual = rh/rt;

% Continuity/radius diagnostics. These make explicit which density and
% velocity were used to size the annulus. The active sizing uses the chosen
% inlet-total convention plus Mx_design. For comparison, ambient_direct uses
% the ISA static atmosphere directly with Vx=Mx*a_amb; it is not used in the
% design unless explicitly coded elsewhere.
a1 = sqrt(air.gamma*air.R*T1);
continuity_active = make_continuity_radius_diagnostic('active_inlet_model', ...
    in.mdot, rho1, Vx, in.hub_to_tip, in.design.mean_radius_rule);
continuity_ambient_direct = make_continuity_radius_diagnostic('ambient_static_direct', ...
    in.mdot, amb.rho, in.Mx_design*amb.a, in.hub_to_tip, in.design.mean_radius_rule);

% ---------- 2) Velocity triangles ----------
status_msg(in,'2/12','velocity triangles');

alpha1_mode = field_or(in.design,'alpha1_mode_active',field_or(in.design,'alpha1_mode','free'));
alpha1_forced = strcmpi(alpha1_mode,'zero') || strcmpi(alpha1_mode,'force_zero') || ...
                strcmpi(alpha1_mode,'fixed_zero');
if alpha1_forced
    alpha1_target_deg = field_or(in.design,'alpha1_deg',0);
    tan_a1 = tand(alpha1_target_deg);
    R = 1 - phi*tan_a1 - psi/2;
else
    tan_a1 = (1 - R - psi/2)/max(phi,1e-12);
end

tan_b2 = (psi - 1 + phi*tan_a1)/max(phi,1e-12);
tan_b1 = tan_a1 - 1/max(phi,1e-12);
tan_a2 = tan_b2 + 1/max(phi,1e-12);

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

Vt1_tip_for_check = spanwise_vtheta(Vt1, rt, rm, in.loss_model.radial.vtheta_mode);
Mrel_tip_actual = hypot(Vx, Vt1_tip_for_check - Utip_actual)/sqrt(air.gamma*air.R*T1);

% Recovery-ratio stability metric. The velocity triangles in this code are
% stored with signed angles in a fixed tangential coordinate system:
%   Vtheta = Vx*tan(alpha), Wtheta = Vx*tan(beta), Vtheta = U + Wtheta.
% The Cumpsty/course derivation writes the rotor-exit absolute swirl as
%   Vtheta2 = U - Vx2*tan(beta2_RR),
% so beta2_RR is the positive geometric rotor-exit relative angle measured
% opposite to the signed Wtheta direction. In contrast, alpha1 is used in the
% derivation as signed inlet absolute swirl, Vtheta1 = Vx1*tan(alpha1).
% Therefore the default stability metric keeps alpha1 signed and converts
% beta2 with beta2_RR = -beta2_code. RR_signed is retained only as a
% coordinate-convention diagnostic.
RR_angle_convention = field_or(in.design,'RR_angle_convention','lecture_signed_alpha_beta2_positive');
RR_signed = recovery_ratio_from_angles(alpha1, beta2, phi);
if strcmpi(RR_angle_convention,'signed')
    alpha1_RR = alpha1;
    beta2_RR  = beta2;
elseif strcmpi(RR_angle_convention,'positive_magnitudes_legacy')
    alpha1_RR = abs(alpha1);
    beta2_RR  = abs(beta2);
else
    alpha1_RR = alpha1;
    beta2_RR  = -beta2;
end
RR_beta2_convention_residual = -beta2_RR; % <=0 only when beta2_RR >= 0
RR = recovery_ratio_from_angles(alpha1_RR, beta2_RR, phi);

% ---------- 3) Euler work and preliminary total states ----------
status_msg(in,'3/12','Euler work and preliminary thermodynamic states');

Delta_h0 = psi * U^2;
T02 = T01 + Delta_h0/air.cp;

% Use the ideal total-pressure rise for the first pass through density- and
% Reynolds-dependent checks. The final achieved beta_tt is recomputed after
% the modular loss model has been evaluated.
beta_tt_ideal_from_work = max((T02/T01)^(air.gamma/(air.gamma-1)), 1.0);
p02 = p01 * beta_tt_ideal_from_work;
T02s = T02s_target; %#ok<NASGU>
w_is = w_is_target; %#ok<NASGU>

% ---------- 4) Rotor / stator static states ----------
status_msg(in,'4/12','rotor and stator static states');

T2 = T02 - V2^2/(2*air.cp);
if T2 <= 1
    error('Rotor-exit static temperature became non-physical.');
end
p2 = p02 / (T02/T2)^(air.gamma/(air.gamma-1));
rho2 = p2/(air.R*T2);
M2 = V2/sqrt(air.gamma*air.R*T2);
Mrel1 = W1/sqrt(air.gamma*air.R*T1);
Mrel2 = W2/sqrt(air.gamma*air.R*T2);

alpha3 = deg2rad(in.alpha3_deg);
Vt3 = Vx*tan(alpha3);
V3  = hypot(Vx, Vt3);
T03 = T02;
p03_prelim = p02;
T3 = T03 - V3^2/(2*air.cp);
if T3 <= 1
    error('Stator-exit static temperature became non-physical.');
end
p3 = p03_prelim / (T03/T3)^(air.gamma/(air.gamma-1));
rho3 = p3/(air.R*T3);
M3 = V3/sqrt(air.gamma*air.R*T3);

% ---------- 5) Reynolds numbers ----------
status_msg(in,'5/12','Reynolds numbers');

mu1 = sutherland_mu(T1);
mu2 = sutherland_mu(T2);

% ---------- 6) Howell and DF/Howell solidity sizing ----------
status_msg(in,'6/12','Howell turning, diffusion-factor checks, and chord closure');

howell = read_howell_curves(in.howell_paths);

Dbeta_rotor_deg   = abs(rad2deg(beta1 - beta2));
Dalpha_stator_deg = abs(rad2deg(alpha2 - alpha3));

DFR_base = 1 - W2/max(W1,1e-12);
DFR_turn = abs(Wt2 - Wt1)/(2*max(W1,1e-12));

DFS_base = 1 - V3/max(V2,1e-12);
DFS_turn = abs(Vt2 - Vt3)/(2*max(V2,1e-12));

rotor_solidity_design = struct('enabled',false,'sigma',in.rotor_solidity,'status','fixed');
stator_solidity_design = struct('enabled',false,'sigma',in.stator_solidity,'status','fixed');

chord_closure_mode = field_or(in.design,'chord_closure_mode','fixed_chord');
chord_closure = struct('mode',chord_closure_mode,'iterations',0, ...
    'rotor_blade_count_active',NaN,'stator_vane_count_active',NaN, ...
    'rotor_t_over_c',NaN,'stator_t_over_c',NaN);

if strcmpi(chord_closure_mode,'blade_count_from_solidity')
    % Reverse geometry closure requested for the MEANGEN-starting-value
    % workflow: DF/Howell gives sigma, the selected integer blade count gives
    % pitch, and chord follows from c = sigma*s.  Re depends on chord, while
    % Howell includes a Reynolds correction, so a small fixed-point loop is
    % used.  The loss-model equations themselves are untouched.
    Zr = field_or(in.design,'rotor_blade_count_active',field_or(in.design,'rotor_blade_count_default',30));
    Zs = field_or(in.design,'stator_vane_count_active',field_or(in.design,'stator_vane_count_default',40));
    Zr = max(1,round(Zr));
    Zs = max(1,round(Zs));

    rotor_t_over_c  = field_or(in.design,'rotor_t_over_c',field_or(in.design,'t_over_c',0.10));
    stator_t_over_c = field_or(in.design,'stator_t_over_c',field_or(in.design,'t_over_c',0.10));

    n_iter = max(1,round(field_or(in.design,'chord_solidity_iterations',4)));
    Cmean = 2*pi*rm;

    % Initial pitch from blade count. Initial chord uses the current/fallback
    % solidity, then the loop updates solidity and chord consistently.
    in.rotor_pitch = Cmean/Zr;
    in.stator_pitch = Cmean/Zs;
    in.rotor_chord = max(in.rotor_solidity*in.rotor_pitch,1e-6);
    in.stator_chord = max(in.stator_solidity*in.stator_pitch,1e-6);
    in.rotor_thickness = rotor_t_over_c*in.rotor_chord;
    in.stator_thickness = stator_t_over_c*in.stator_chord;

    for k_geom = 1:n_iter
        Re_rotor  = rho1 * W1 * in.rotor_chord / mu1;
        Re_stator = rho2 * V2 * in.stator_chord / mu2;

        if isfield(in.design,'solidity_mode') && strcmpi(in.design.solidity_mode,'from_DF_Howell')
            rotor_solidity_design = size_solidity_from_DF_Howell( ...
                'rotor', Dbeta_rotor_deg, abs(rad2deg(beta2)), Re_rotor, ...
                DFR_base, DFR_turn, howell, in);

            stator_solidity_design = size_solidity_from_DF_Howell( ...
                'stator', Dalpha_stator_deg, abs(rad2deg(alpha3)), Re_stator, ...
                DFS_base, DFS_turn, howell, in);

            in.rotor_solidity = rotor_solidity_design.sigma;
            in.stator_solidity = stator_solidity_design.sigma;
        end

        in.rotor_pitch = Cmean/Zr;
        in.stator_pitch = Cmean/Zs;
        in.rotor_chord = max(in.rotor_solidity*in.rotor_pitch,1e-6);
        in.stator_chord = max(in.stator_solidity*in.stator_pitch,1e-6);
        in.rotor_thickness = rotor_t_over_c*in.rotor_chord;
        in.stator_thickness = stator_t_over_c*in.stator_chord;
    end

    % Final Reynolds numbers and final solidity update at the converged chord.
    Re_rotor  = rho1 * W1 * in.rotor_chord / mu1;
    Re_stator = rho2 * V2 * in.stator_chord / mu2;
    if isfield(in.design,'solidity_mode') && strcmpi(in.design.solidity_mode,'from_DF_Howell')
        rotor_solidity_design = size_solidity_from_DF_Howell( ...
            'rotor', Dbeta_rotor_deg, abs(rad2deg(beta2)), Re_rotor, ...
            DFR_base, DFR_turn, howell, in);
        stator_solidity_design = size_solidity_from_DF_Howell( ...
            'stator', Dalpha_stator_deg, abs(rad2deg(alpha3)), Re_stator, ...
            DFS_base, DFS_turn, howell, in);
        in.rotor_solidity = rotor_solidity_design.sigma;
        in.stator_solidity = stator_solidity_design.sigma;
        in.rotor_chord = max(in.rotor_solidity*in.rotor_pitch,1e-6);
        in.stator_chord = max(in.stator_solidity*in.stator_pitch,1e-6);
        in.rotor_thickness = rotor_t_over_c*in.rotor_chord;
        in.stator_thickness = stator_t_over_c*in.stator_chord;
        Re_rotor  = rho1 * W1 * in.rotor_chord / mu1;
        Re_stator = rho2 * V2 * in.stator_chord / mu2;
    end

    chord_closure.mode = 'blade_count_from_solidity';
    chord_closure.iterations = n_iter;
    chord_closure.rotor_blade_count_active = Zr;
    chord_closure.stator_vane_count_active = Zs;
    chord_closure.rotor_t_over_c = rotor_t_over_c;
    chord_closure.stator_t_over_c = stator_t_over_c;

else
    % Original closure: chord is a user input, solidity gives pitch and blade
    % count follows from circumference/pitch.
    Re_rotor  = rho1 * W1 * in.rotor_chord / mu1;
    Re_stator = rho2 * V2 * in.stator_chord / mu2;

    if isfield(in.design,'solidity_mode') && strcmpi(in.design.solidity_mode,'from_DF_Howell')

        rotor_solidity_design = size_solidity_from_DF_Howell( ...
            'rotor', Dbeta_rotor_deg, abs(rad2deg(beta2)), Re_rotor, ...
            DFR_base, DFR_turn, howell, in);

        stator_solidity_design = size_solidity_from_DF_Howell( ...
            'stator', Dalpha_stator_deg, abs(rad2deg(alpha3)), Re_stator, ...
            DFS_base, DFS_turn, howell, in);

        in.rotor_solidity = rotor_solidity_design.sigma;
        in.stator_solidity = stator_solidity_design.sigma;

        in.rotor_pitch = in.rotor_chord / in.rotor_solidity;
        in.stator_pitch = in.stator_chord / in.stator_solidity;
    end
end

f_beta_rotor  = howell.f_beta2(abs(rad2deg(beta2)));
f_beta_stator = howell.f_beta2(abs(rad2deg(alpha3)));
phi_Re_rotor  = howell.Phi_Re(Re_rotor / 1e5);
phi_Re_stator = howell.Phi_Re(Re_stator / 1e5);

sc_rotor  = 1 / in.rotor_solidity;
sc_stator = 1 / in.stator_solidity;

psi_sc_rotor  = howell.Psi_sc(sc_rotor);
psi_sc_stator = howell.Psi_sc(sc_stator);

D_rotor_star  = f_beta_rotor  * phi_Re_rotor  * psi_sc_rotor;
D_stator_star = f_beta_stator * phi_Re_stator * psi_sc_stator;

howell_rotor_residual_deg  = Dbeta_rotor_deg   - D_rotor_star;
howell_stator_residual_deg = Dalpha_stator_deg - D_stator_star;

howell_ok_rotor  = howell_rotor_residual_deg  <= 0;
howell_ok_stator = howell_stator_residual_deg <= 0;

DF_rotor_generic  = DFR_base + DFR_turn/in.rotor_solidity;
DF_stator_generic = DFS_base + DFS_turn/in.stator_solidity;

DF_rotor_lieblein  = DF_rotor_generic;
DF_stator_lieblein = DF_stator_generic;

DF_rotor_residual  = DF_rotor_lieblein  - in.design.DF_limit;
DF_stator_residual = DF_stator_lieblein - in.design.DF_limit;

DF_rotor_low_residual  = in.design.DF_min - DF_rotor_lieblein;
DF_stator_low_residual = in.design.DF_min - DF_stator_lieblein;

DF_rotor_ok  = DF_rotor_residual <= 0 && DF_rotor_low_residual <= 0;
DF_stator_ok = DF_stator_residual <= 0 && DF_stator_low_residual <= 0;

% ---------- 7) Project and duty-coefficient constraints ----------
status_msg(in,'7/12','hard project constraints');

phi_min = field_or(in.design,'phi_min',-inf);
phi_max = field_or(in.design,'phi_max', inf);
psi_min = field_or(in.design,'psi_min',-inf);
psi_max = field_or(in.design,'psi_max', inf);
R_min = field_or(in.design,'R_min',0);
R_max = field_or(in.design,'R_max',1);

phi_low_residual  = phi_min - phi;
phi_high_residual = phi - phi_max;
psi_low_residual  = psi_min - psi;
psi_high_residual = psi - psi_max;
R_low_residual    = R_min - R;
R_high_residual   = R - R_max;

mdot_error = mdot_meanline - in.mdot;
mdot_rel_error = mdot_error/max(in.mdot,1e-12);

fixed_tol = field_or(in.design,'fixed_tolerance_rel',0.02);
Mx_residual = rel_dev_residual(Mx_actual, in.Mx_design, fixed_tol);
Mrel_tip_residual = rel_dev_residual(Mrel_tip_actual, in.Mrel_tip_design, fixed_tol);
htr_residual = rel_dev_residual(hub_to_tip_actual, in.hub_to_tip, fixed_tol);
mdot_residual = rel_dev_residual(mdot_meanline, in.mdot, fixed_tol);
altitude_residual = -fixed_tol;
Mflight_residual = -fixed_tol;
Nrpm_residual = -fixed_tol;

% ---------- 8) Lieblein optimum incidence/deviation estimates ----------
status_msg(in,'8/12','Lieblein incidence/deviation');

lieblein = struct('enabled',false);
if isfield(in,'lieblein') && isfield(in.lieblein,'enabled') && in.lieblein.enabled
    try
        lieblein_curves = read_lieblein_curves(in.lieblein_paths);
        rotor_tmax_over_c = row_tmax_over_c(in,'rotor');
        stator_tmax_over_c = row_tmax_over_c(in,'stator');
        lieblein = struct();
        lieblein.enabled = true;
        lieblein.rotor = lieblein_row_estimate(lieblein_curves, 'rotor', ...
            abs(rad2deg(beta1)), abs(rad2deg(beta2)), in.rotor_solidity, ...
            rotor_tmax_over_c, in.lieblein);
        lieblein.stator = lieblein_row_estimate(lieblein_curves, 'stator', ...
            abs(rad2deg(alpha2)), abs(rad2deg(alpha3)), in.stator_solidity, ...
            stator_tmax_over_c, in.lieblein);
        lieblein.curve_files = lieblein_curves.files;
        lieblein.curve_ranges = lieblein_curves.ranges;
        lieblein = add_signed_lieblein_geometry(lieblein, beta1, beta2, alpha2, alpha3);
    catch ME
        lieblein = struct('enabled',false,'error',ME.message);
        warning('Lieblein incidence/deviation evaluation failed: %s', ME.message);
    end
end

lieblein_rotor_residual = -1;
lieblein_stator_residual = -1;
if isfield(in,'lieblein') && isfield(in.lieblein,'enabled') && in.lieblein.enabled
    if ~(isfield(lieblein,'enabled') && lieblein.enabled && ...
            isfield(lieblein,'rotor') && isfield(lieblein.rotor,'valid') && lieblein.rotor.valid)
        lieblein_rotor_residual = 1;
    end
    if ~(isfield(lieblein,'enabled') && lieblein.enabled && ...
            isfield(lieblein,'stator') && isfield(lieblein.stator,'valid') && lieblein.stator.valid)
        lieblein_stator_residual = 1;
    end
end

in_loss = apply_lieblein_to_loss_model(in, lieblein);

% ---------- 9) Literature-based modular loss build-up ----------
status_msg(in,'9/12','loss-model call');

span = rt - rh;
loss = estimate_losses(in_loss, struct( ...
    'air',air,'span',span,'rt',rt,'rh',rh,'rm',rm,'omega',omega,'Vx',Vx, ...
    'alpha1',alpha1,'alpha2',alpha2,'alpha3',alpha3, ...
    'beta1',beta1,'beta2',beta2,'Vt1',Vt1,'Vt2',Vt2,'Vt3',Vt3, ...
    'V1',V1,'V2',V2,'V3',V3,'W1',W1,'W2',W2, ...
    'M2',M2,'M3',M3,'Mrel1',Mrel1,'Mrel2',Mrel2,'T1',T1,'T2',T2,'T3',T3, ...
    'rotor_pitch',in.rotor_pitch,'stator_pitch',in.stator_pitch, ...
    'rotor_chord',in.rotor_chord,'stator_chord',in.stator_chord, ...
    'rotor_gap',in.rotor_tip_clearance,'stator_gap',in.stator_hub_clearance));

% ---------- 10) Stage performance ----------
status_msg(in,'10/12','performance estimate');

zetaR = loss.zetaR;
zetaS = loss.zetaS;
loss_specific = (zetaS*V3^2 + zetaR*W2^2)/2;

% Row/source-level loss-to-efficiency bookkeeping.  The stage efficiency
% estimate below is eta ~= 1 - loss_specific/Delta_h0, so each zeta source
% has an explicit eta penalty.  This is useful for sensitivity checks, e.g.
% quantifying how much the selected endwall model changes eta_tt.
eta_penalty = loss_efficiency_penalty_breakdown(loss, W2, V3, Delta_h0);
eta_penalty_endwall = eta_penalty.rotor.endwall + eta_penalty.stator.endwall;

% Reader Eq. form: eta ~= 1 - (row loss terms)/(psi U^2).  This makes the
% achieved pressure ratio a consequence of the swept work coefficient rather
% than an imposed closure.
eta_tt_est = 1 - loss_specific/max(Delta_h0,eps);
eta_tt_no_endwall_est = min(1, eta_tt_est + eta_penalty_endwall);
w_is_est = max(eta_tt_est,0) * Delta_h0;
w_is_no_endwall_est = max(eta_tt_no_endwall_est,0) * Delta_h0;
beta_tt_est = (1 + w_is_est/(air.cp*T01))^(air.gamma/(air.gamma-1));
beta_tt_no_endwall_est = (1 + w_is_no_endwall_est/(air.cp*T01))^(air.gamma/(air.gamma-1));
p03_est = p01 * beta_tt_est;
Power = in.mdot * Delta_h0;

beta_min = field_or(in.design,'beta_tt_min',-inf);
beta_max = field_or(in.design,'beta_tt_max', inf);
beta_tt_min_residual = beta_min - beta_tt_est;
beta_tt_max_residual = beta_tt_est - beta_max;

RR_soft_min = field_or(in.design,'RR_soft_min',0.50);
RR_soft_residual = RR_soft_min - RR;

% ---------- 11) Aerodynamic and centrifugal force estimates ----------
status_msg(in,'11/12','forces, chords, and constraints');

Torque = in.mdot * rm * (Vt2 - Vt1);
Ft_tan = Torque / max(rm,1e-8);
rotor_blades = max(1, round(2*pi*rm / in.rotor_pitch));
stator_vanes = max(1, round(2*pi*rm / in.stator_pitch));
blade_volume = in.rotor_chord * span * in.rotor_thickness * in.blade_height_factor;
blade_mass = in.material.rho_blade * blade_volume;
Fc_blade = blade_mass * omega^2 * rm;

[rotor_axial_chord, rotor_stagger_for_Cx_deg] = axial_chord_projection( ...
    in.rotor_chord, beta1, beta2, lieblein, 'rotor');

[stator_axial_chord, stator_stagger_for_Cx_deg] = axial_chord_projection( ...
    in.stator_chord, alpha2, alpha3, lieblein, 'stator');

residuals = struct();
residuals.howell_rotor_deg  = howell_rotor_residual_deg;
residuals.howell_stator_deg = howell_stator_residual_deg;
residuals.DF_rotor          = DF_rotor_residual;
residuals.DF_stator         = DF_stator_residual;
residuals.DF_rotor_low      = DF_rotor_low_residual;
residuals.DF_stator_low     = DF_stator_low_residual;
residuals.phi_min           = phi_low_residual;
residuals.phi_max           = phi_high_residual;
residuals.psi_min           = psi_low_residual;
residuals.psi_max           = psi_high_residual;
residuals.R_min             = R_low_residual;
residuals.R_max             = R_high_residual;
residuals.beta_tt_min       = beta_tt_min_residual;
residuals.beta_tt_max       = beta_tt_max_residual;
residuals.beta_tt_excess    = max([0; beta_tt_min_residual; beta_tt_max_residual]);
residuals.lieblein_rotor    = lieblein_rotor_residual;
residuals.lieblein_stator   = lieblein_stator_residual;
residuals.Mx_design         = Mx_residual;
residuals.Mrel_tip_design   = Mrel_tip_residual;
residuals.Mrel_tip_excess_rel = max(0, abs((Mrel_tip_actual - in.Mrel_tip_design)/max(in.Mrel_tip_design,1e-12)) - fixed_tol);
residuals.Mrel_tip_radius_match = Mrel_tip_residual;
residuals.hub_to_tip        = htr_residual;
residuals.mdot              = mdot_residual;
residuals.mdot_abs          = mdot_error;       % signed kg/s error, for readable printout
residuals.mdot_rel          = mdot_rel_error;   % signed relative error, not tolerance-subtracted
residuals.mdot_excess_rel   = max(0, abs(mdot_rel_error) - fixed_tol);
residuals.altitude          = altitude_residual;
residuals.M_flight          = Mflight_residual;
residuals.N_rpm             = Nrpm_residual;

% Separate true aerodynamic hard constraints from project-matching diagnostics.
% The preliminary meanline sweep is meant to find useful MEANGEN starting
% values. It should not declare every candidate infeasible merely because a
% one-stage meanline estimate misses beta_tt, mdot or Mrel_tip by a few percent.
hard_ineq = [howell_rotor_residual_deg; howell_stator_residual_deg; ...
    DF_rotor_residual; DF_stator_residual; DF_rotor_low_residual; DF_stator_low_residual; ...
    phi_low_residual; phi_high_residual; psi_low_residual; psi_high_residual; ...
    R_low_residual; R_high_residual; lieblein_rotor_residual; lieblein_stator_residual];

beta_ineq = [beta_tt_min_residual; beta_tt_max_residual];
project_ineq = [Mx_residual; Mrel_tip_residual; htr_residual; mdot_residual; ...
    altitude_residual; Mflight_residual; Nrpm_residual];

if field_or(in.design,'enforce_beta_tt_range_as_hard',false)
    hard_ineq = [hard_ineq; beta_ineq];
end
if field_or(in.design,'enforce_project_matching_as_hard',false)
    hard_ineq = [hard_ineq; project_ineq];
end
if field_or(in.design,'enforce_Mrel_tip_as_hard',false)
    hard_ineq = [hard_ineq; Mrel_tip_residual];
end
if field_or(in.design,'enforce_RR_as_hard',false)
    hard_ineq = [hard_ineq; RR_soft_residual];
    if ~strcmpi(RR_angle_convention,'signed')
        hard_ineq = [hard_ineq; RR_beta2_convention_residual];
    end
end

residuals.hard_ineq = hard_ineq;
residuals.beta_ineq = beta_ineq;
residuals.project_ineq = project_ineq;
residuals.all_ineq = hard_ineq;  % backward-compatible: feasibility now means hard aero feasibility
residuals.all_diagnostics_ineq = [hard_ineq; beta_ineq; project_ineq];
residuals.soft_RR = RR_soft_residual;
residuals.RR_beta2_convention = RR_beta2_convention_residual;
residuals.RR_metric = RR;
residuals.RR_signed = RR_signed;

is_feasible = all(isfinite(residuals.hard_ineq)) && all(residuals.hard_ineq <= 0);
if in.design.enforce_hard_limits && ~is_feasible
    error('Design violated hard Howell, diffusion-factor, duty-coefficient, or Lieblein feasibility limits.');
end

% ---------- 12) Package outputs ----------
status_msg(in,'12/12','packaging output');

out.input = in;
out.ambient = amb;
out.station = [ ...
    make_station('1',T1,p1,rho1,Vx,Vt1,V1,alpha1,NaN,M0,T01,p01), ...
    make_station('2',T2,p2,rho2,Vx,Vt2,V2,alpha2,beta2,M2,T02,p02), ...
    make_station('3',T3,p3,rho3,Vx,Vt3,V3,alpha3,NaN,M3,T03,p03_est) ...
    ];
out.geometry = struct('rt',rt,'rh',rh,'rm',rm,'A1',A1, ...
    'rm_continuity',rm_continuity,'A1_continuity',A1_continuity, ...
    'rt_continuity',rt_continuity,'rh_continuity',rh_continuity, ...
    'tipmach_radius_match',tipmach_radius_match, ...
    'rotor_pitch',in.rotor_pitch,'stator_pitch',in.stator_pitch, ...
    'rotor_chord',in.rotor_chord,'stator_chord',in.stator_chord, ...
    'rotor_thickness',in.rotor_thickness,'stator_thickness',field_or(in,'stator_thickness',NaN), ...
    'rotor_t_over_c',in.rotor_thickness/max(in.rotor_chord,1e-12), ...
    'stator_t_over_c',field_or(in,'stator_thickness',NaN)/max(in.stator_chord,1e-12), ...
    'rotor_solidity',in.rotor_solidity,'stator_solidity',in.stator_solidity, ...
    'rotor_axial_chord',rotor_axial_chord, ...
    'stator_axial_chord',stator_axial_chord, ...
    'rotor_stagger_for_axial_chord_deg',rotor_stagger_for_Cx_deg, ...
    'stator_stagger_for_axial_chord_deg',stator_stagger_for_Cx_deg, ...
    'rotor_tip_clearance',in.rotor_tip_clearance, ...
    'stator_hub_clearance',in.stator_hub_clearance, ...
    'rotor_blades',rotor_blades,'stator_vanes',stator_vanes, ...
    'chord_closure',chord_closure, ...
    'span_sections',{{'hub','mean','tip'}},'r_sections',[rh rm rt], ...
    'solve_mode',in.design.solve_mode,'solve_mode_note',geometry_mode_note, ...
    'fixed_mean_radius_duty_active',fixed_radius_duty_active, ...
    'inlet_total_mode',field_or(in.design,'inlet_total_mode','freestream_total'), ...
    'inlet_total_note',inlet_total_note);

out.coeffs = struct('phi_initial',phi_initial,'phi_actual',phi_actual,'phi_used',phi, ...
    'psi',psi,'R',R,'RR',RR,'RR_signed',RR_signed, ...
    'RR_angle_convention',RR_angle_convention, ...
    'alpha1_RR_deg',rad2deg(alpha1_RR),'beta2_RR_deg',rad2deg(beta2_RR), ...
    'beta_tt_target',in.beta_tt_target, ...
    'beta_tt_ideal_from_work',beta_tt_ideal_from_work,'beta_tt_est',beta_tt_est, ...
    'alpha1_mode',alpha1_mode);
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
    'limit',in.design.DF_limit,'minimum',in.design.DF_min, ...
    'rotor_residual',DF_rotor_residual,'stator_residual',DF_stator_residual, ...
    'rotor_low_residual',DF_rotor_low_residual,'stator_low_residual',DF_stator_low_residual, ...
    'rotor_ok',DF_rotor_ok,'stator_ok',DF_stator_ok, ...
    'rotor_solidity_design',rotor_solidity_design, ...
    'stator_solidity_design',stator_solidity_design);
out.design_basis = struct('Mx_target',in.Mx_design,'Mx_actual',Mx_actual, ...
    'Mrel_tip_target',in.Mrel_tip_design,'Mrel_tip_actual',Mrel_tip_actual, ...
    'mdot_target',in.mdot,'mdot_meanline',mdot_meanline, ...
    'mdot_error',mdot_error,'mdot_rel_error',mdot_rel_error, ...
    'fixed_tolerance_rel',fixed_tol, ...
    'hub_to_tip_actual',hub_to_tip_actual,'hub_to_tip_input',in.hub_to_tip, ...
    'eta_work_guess',field_or(in.design,'eta_work_guess',NaN), ...
    'Delta_h0_from_phi_psi',Delta_h0, ...
    'phi_min',phi_min,'phi_max',phi_max,'psi_min',psi_min,'psi_max',psi_max, ...
    'R_min',R_min,'R_max',R_max, ...
    'Utip',Utip_actual,'Vt1_tip_for_check',Vt1_tip_for_check, ...
    'rho_radius_eval',rho1,'T_radius_eval',T1,'p_radius_eval',p1, ...
    'a_radius_eval',a1,'Vx_radius_eval',Vx,'A_radius_eval',A1, ...
    'rt_radius_eval',rt,'rh_radius_eval',rh,'rm_radius_eval',rm, ...
    'rt_continuity',rt_continuity,'rh_continuity',rh_continuity, ...
    'rm_continuity',rm_continuity,'A1_continuity',A1_continuity, ...
    'tipmach_radius_match',tipmach_radius_match, ...
    'mean_radius_rule',in.design.mean_radius_rule, ...
    'fixed_mean_radius_duty_active',fixed_radius_duty_active, ...
    'T01',T01,'p01',p01,'inlet_total_mode',field_or(in.design,'inlet_total_mode','freestream_total'), ...
    'ambient_T',amb.T,'ambient_p',amb.p,'ambient_rho',amb.rho,'ambient_a',amb.a, ...
    'continuity_active',continuity_active,'continuity_ambient_direct',continuity_ambient_direct, ...
    'altitude_m',in.altitude_m,'M_flight',in.M_flight,'N_rpm',in.N_rpm);
out.stability = struct('RR',RR,'RR_signed',RR_signed, ...
    'RR_angle_convention',RR_angle_convention, ...
    'alpha1_RR_deg',rad2deg(alpha1_RR),'beta2_RR_deg',rad2deg(beta2_RR), ...
    'alpha1_signed_deg',rad2deg(alpha1),'beta2_signed_deg',rad2deg(beta2), ...
    'beta2_RR_positive_residual',RR_beta2_convention_residual, ...
    'beta2_RR_positive_ok',RR_beta2_convention_residual <= 0, ...
    'RR_soft_min',RR_soft_min, ...
    'RR_soft_residual',RR_soft_residual,'RR_soft_ok',RR_soft_residual <= 0, ...
    'interpretation','RR closer to 1 is favourable; RR near 0 is adverse. Default RR keeps alpha1 signed and uses beta2_RR = -beta2_code, matching the course derivation Vtheta2 = U - Vx*tan(beta2). RR_signed is only a coordinate-convention diagnostic.');
out.lieblein = lieblein;
out.loss = loss;
out.loss.loss_model_used = in_loss.loss_model;
out.loss.eta_penalty = eta_penalty;
out.loss.eta_penalty_endwall = eta_penalty_endwall;
out.performance = struct('Delta_h0',Delta_h0,'loss_specific',loss_specific, ...
    'w_is_est',w_is_est,'w_is_target_at_beta_tt_target',w_is_target, ...
    'eta_tt_est',eta_tt_est,'beta_tt_est',beta_tt_est, ...
    'eta_tt_no_endwall_est',eta_tt_no_endwall_est, ...
    'beta_tt_no_endwall_est',beta_tt_no_endwall_est, ...
    'eta_penalty_endwall',eta_penalty_endwall, ...
    'power_W',Power,'torque_Nm',Torque,'tangential_force_N',Ft_tan, ...
    'centrifugal_force_per_blade_N',Fc_blade);
out.constraints = struct('residuals',residuals,'is_feasible',is_feasible);
out.flags = struct('DF_rotor_high',~DF_rotor_ok, 'DF_stator_high',~DF_stator_ok, ...
    'Mrel1_high',Mrel1 > 1.4, 'Mrel2_high',Mrel2 > 1.4, 'M2_high',M2 > 1.0);

status_msg(in,'DONE','meanline_fan_design completed');
end


function out_best = meanline_design_sweep_project(in)
% Sweep duty coefficients to find a good starting point for MEANGEN.  In the
% continuity-fixed project mode phi is computed from geometry.  In the
% Denton/MEANGEN-style fixed-radius mode phi and psi are swept at a fixed
% design mean radius, and alpha1=0 gives R=1-psi/2.

status_msg(in,'SWEEP','started project-constrained duty sweep');

use_continuity_fixed_radius = field_or(in.design,'use_continuity_fixed_radius',true);
use_fixed_mean_radius_duty = strcmpi(field_or(in.design,'solve_mode',''),'sweep_fixed_mean_radius_duty') || ...
    field_or(in.design,'fixed_mean_radius_duty_active',false);
if use_fixed_mean_radius_duty
    phi_grid = linspace(in.design.phi_min, in.design.phi_max, ...
                        field_or(in.design,'phi_search_points',21));
elseif use_continuity_fixed_radius
    phi_grid = NaN;  % placeholder; candidate evaluator computes phi from continuity-fixed r_m
else
    phi_grid = linspace(in.design.phi_min, in.design.phi_max, ...
                        field_or(in.design,'phi_search_points',21));
    phi_grid = augment_phi_grid_for_mdot_target(phi_grid, in);
end

psi_grid = linspace(in.design.psi_min, in.design.psi_max, ...
                    field_or(in.design,'psi_search_points',21));
R_grid   = linspace(in.design.R_min, in.design.R_max, ...
                    field_or(in.design,'R_search_points',41));

mode_request = field_or(in.design,'alpha1_mode','both');
switch lower(mode_request)
    case {'both','compare'}
        alpha_modes = {'free','zero'};
    case {'zero','force_zero','fixed_zero'}
        alpha_modes = {'zero'};
    otherwise
        alpha_modes = {'free'};
end

blade_count_selection_mode = field_or(in.design,'blade_count_selection_mode','postprocess_table_only');
use_blade_count_chord_closure = strcmpi(field_or(in.design,'chord_closure_mode','fixed_chord'),'blade_count_from_solidity');

if use_blade_count_chord_closure && strcmpi(blade_count_selection_mode,'sweep_in_objective') && ...
        field_or(in.design,'blade_count_sweep',true)
    % Legacy behavior: blade count participates in the optimizer objective.
    % This can introduce artificial Z preferences when the meanline loss model
    % has weak/neutral blade-count sensitivity.
    rotorZ_grid = field_or(in.design,'rotor_blade_count_range',15:40);
    statorZ_grid = field_or(in.design,'stator_vane_count_range',20:60);
elseif use_blade_count_chord_closure
    % Recommended behavior: do not let Z participate in objective selection.
    % Use a single reference count only to close the chord/Re-dependent
    % preliminary loss estimate; after the best duty/solidity candidate is
    % selected, generate a post-processing table for all Z_R/Z_S combinations.
    rotorZ_grid = field_or(in.design,'rotor_blade_count_default',30);
    statorZ_grid = field_or(in.design,'stator_vane_count_default',40);
else
    rotorZ_grid = NaN;
    statorZ_grid = NaN;
end
rotorZ_grid = rotorZ_grid(:).';
statorZ_grid = statorZ_grid(:).';

n_free = numel(phi_grid)*numel(psi_grid)*numel(R_grid)*numel(rotorZ_grid)*numel(statorZ_grid);
n_zero = numel(phi_grid)*numel(psi_grid)*numel(rotorZ_grid)*numel(statorZ_grid);
n_total = 0;
for im = 1:numel(alpha_modes)
    if strcmpi(alpha_modes{im},'zero')
        n_total = n_total + n_zero;
    else
        n_total = n_total + n_free;
    end
end

history = struct();
history.phi = nan(n_total,1);
history.psi = nan(n_total,1);
history.R = nan(n_total,1);
history.alpha1_mode = cell(n_total,1);
history.rotor_blades = nan(n_total,1);
history.stator_vanes = nan(n_total,1);
history.rotor_chord = nan(n_total,1);
history.stator_chord = nan(n_total,1);
history.rotor_solidity = nan(n_total,1);
history.stator_solidity = nan(n_total,1);
history.mdot = nan(n_total,1);
history.Mrel_tip = nan(n_total,1);
history.score = inf(n_total,1);
history.eta = nan(n_total,1);
history.beta_tt = nan(n_total,1);
history.RR = nan(n_total,1);
history.feasible = false(n_total,1);
history.max_residual = inf(n_total,1);
history.message = cell(n_total,1);

best_feasible_score = inf;
best_any_score = inf;
out_best_feasible = [];
out_best_any = [];

count = 0;
status_every = max(1,field_or(in.design,'status_every',100));

for im = 1:numel(alpha_modes)
    mode = alpha_modes{im};

    if strcmpi(mode,'zero')
        R_loop = NaN;
    else
        R_loop = R_grid;
    end

    for ip = 1:numel(phi_grid)
        for iw = 1:numel(psi_grid)
            for ir = 1:numel(R_loop)
                for izr = 1:numel(rotorZ_grid)
                    for izs = 1:numel(statorZ_grid)
                        count = count + 1;

                        trial = in;
                        trial.design.solve_mode = 'project_fixed_candidate';
                        trial.design.sweep_active = true;
                        trial.design.fixed_mean_radius_duty_active = use_fixed_mean_radius_duty;
                        trial.design.alpha1_mode_active = mode;
                        trial.design.status_print = false;
                        if use_fixed_mean_radius_duty || ~use_continuity_fixed_radius
                            trial.phi = phi_grid(ip);
                        end
                        trial.psi = psi_grid(iw);

                        if strcmpi(mode,'zero')
                            trial.R = 1 - trial.psi/2;
                        else
                            trial.R = R_loop(ir);
                        end

                        if ~isnan(rotorZ_grid(izr))
                            trial.design.rotor_blade_count_active = rotorZ_grid(izr);
                            trial.design.stator_vane_count_active = statorZ_grid(izs);
                        end

                        try
                            cand = meanline_fan_design(trial);
                            c = cand.constraints.residuals.all_ineq(:);
                            feasible = all(isfinite(c)) && all(c <= 0);
                            max_res = max(c);
                            J = meanline_preliminary_objective(cand, in);
                            score_any = J + hard_constraint_penalty(c);

                            history.phi(count) = cand.coeffs.phi_used;
                            history.psi(count) = cand.coeffs.psi;
                            history.R(count) = cand.coeffs.R;
                            history.alpha1_mode{count} = mode;
                            history.rotor_blades(count) = cand.geometry.rotor_blades;
                            history.stator_vanes(count) = cand.geometry.stator_vanes;
                            history.rotor_chord(count) = cand.geometry.rotor_chord;
                            history.stator_chord(count) = cand.geometry.stator_chord;
                            history.rotor_solidity(count) = cand.geometry.rotor_solidity;
                            history.stator_solidity(count) = cand.geometry.stator_solidity;
                            history.mdot(count) = cand.design_basis.mdot_meanline;
                            history.Mrel_tip(count) = cand.design_basis.Mrel_tip_actual;
                            history.score(count) = score_any;
                            history.eta(count) = cand.performance.eta_tt_est;
                            history.beta_tt(count) = cand.coeffs.beta_tt_est;
                            history.RR(count) = cand.coeffs.RR;
                            history.feasible(count) = feasible;
                            history.max_residual(count) = max_res;
                            history.message{count} = '';

                            if feasible && J < best_feasible_score
                                best_feasible_score = J;
                                out_best_feasible = cand;
                            end

                            if score_any < best_any_score
                                best_any_score = score_any;
                                out_best_any = cand;
                            end

                        catch ME
                            if use_fixed_mean_radius_duty || ~use_continuity_fixed_radius
                                history.phi(count) = trial.phi;
                            end
                            history.psi(count) = trial.psi;
                            history.R(count) = trial.R;
                            history.alpha1_mode{count} = mode;
                            if ~isnan(rotorZ_grid(izr))
                                history.rotor_blades(count) = rotorZ_grid(izr);
                                history.stator_vanes(count) = statorZ_grid(izs);
                            end
                            history.message{count} = ME.message;
                        end

                        if field_or(in.design,'status_print',true) && ...
                                (mod(count,status_every)==0 || count==1 || count==n_total)
                            if isempty(out_best_feasible)
                                best_note = sprintf('no feasible yet; best violation %.3g', best_any_score);
                            else
                                best_note = sprintf('best feasible J %.6g', best_feasible_score);
                            end
                            fprintf('[meanline][SWEEP] %d/%d candidates evaluated; %s\n', ...
                                count, n_total, best_note);
                        end
                    end
                end
            end
        end
    end
end

% Trim preallocated history in case a loop exited early in a future edit.
fn = fieldnames(history);
for i = 1:numel(fn)
    v = history.(fn{i});
    if size(v,1) >= count
        history.(fn{i}) = v(1:count,:);
    end
end

if ~isempty(out_best_feasible)
    out_best = out_best_feasible;
    selected_source = 'best_feasible';
else
    if isempty(out_best_any)
        error('Project duty/blade-count sweep failed: no candidate could be evaluated.');
    end
    out_best = out_best_any;
    selected_source = 'least_violating_no_feasible_candidate';
end

out_best.sweep = history;
out_best.sweep.selected_source = selected_source;
out_best.sweep.objective = meanline_preliminary_objective(out_best, in);
out_best.sweep.best_any_penalized_score = best_any_score;
out_best.sweep.best_feasible_objective = best_feasible_score;
out_best.sweep.feasible_count = nnz(history.feasible);
out_best.sweep.total_count = count;
out_best.sweep.use_continuity_fixed_radius = use_continuity_fixed_radius;
out_best.sweep.use_fixed_mean_radius_duty = use_fixed_mean_radius_duty;
out_best.sweep.rotor_blade_count_range = rotorZ_grid;
out_best.sweep.stator_vane_count_range = statorZ_grid;
out_best.sweep.blade_count_selection_mode = blade_count_selection_mode;

if use_blade_count_chord_closure && strcmpi(blade_count_selection_mode,'postprocess_table_only')
    out_best.blade_count_table = make_blade_count_geometry_table(out_best, in);
end

% A real beta_tt-mdot speedline is intentionally not generated by this
% preliminary on-design sweep. Such a speedline requires fixed-geometry
% off-design operating points at the same shaft speed (typically from
% MULTALL/back-pressure or an off-design meanline solver). The current sweep
% changes duty variables and may change geometry/radius, so it cannot locate
% choke reliably.
out_best.speedline = struct('computed',false, ...
    'status','not_computed_fixed_geometry_off_design_required', ...
    'minimum_points_required',5, ...
    'reason',['A true compressor speedline/choke map requires at least several fixed-geometry off-design ', ...
              'solutions at the same shaft speed. The preliminary meanline sweep changes design variables ', ...
              'and is not a throttle/back-pressure continuation.']);

status_msg(in,'SWEEP','completed project-constrained duty sweep');
end


function T = make_blade_count_geometry_table(out, in)
%MAKE_BLADE_COUNT_GEOMETRY_TABLE Generate geometry options without using
%blade count in the meanline objective.  The selected aerodynamic candidate
%sets rm, solidity and stagger; each integer blade-count pair is then mapped
%to pitch/chord/axial chord/Re/thickness for MEANGEN/STAGEN screening.

rotorZ = field_or(in.design,'rotor_blade_count_range',15:40);
statorZ = field_or(in.design,'stator_vane_count_range',20:60);
rotorZ = rotorZ(:).';
statorZ = statorZ(:).';

n = numel(rotorZ)*numel(statorZ);
ZR = nan(n,1); ZS = nan(n,1);
sR = nan(n,1); sS = nan(n,1);
cR = nan(n,1); cS = nan(n,1);
CxR = nan(n,1); CxS = nan(n,1);
tR = nan(n,1); tS = nan(n,1);
ReR = nan(n,1); ReS = nan(n,1);
AR_R = nan(n,1); AR_S = nan(n,1);

rm = out.geometry.rm;
span = out.geometry.rt - out.geometry.rh;
sigmaR = out.geometry.rotor_solidity;
sigmaS = out.geometry.stator_solidity;
stagR = out.geometry.rotor_stagger_for_axial_chord_deg;
stagS = out.geometry.stator_stagger_for_axial_chord_deg;
tcR = out.geometry.rotor_t_over_c;
tcS = out.geometry.stator_t_over_c;

rho1 = out.station(1).rho;
rho2 = out.station(2).rho;
T1 = out.station(1).T;
T2 = out.station(2).T;
mu1 = sutherland_mu(T1);
mu2 = sutherland_mu(T2);
W1 = out.velocities.W1;
V2 = out.velocities.V2;

k = 0;
for i = 1:numel(rotorZ)
    for j = 1:numel(statorZ)
        k = k + 1;
        ZR(k) = rotorZ(i);
        ZS(k) = statorZ(j);
        sR(k) = 2*pi*rm/ZR(k);
        sS(k) = 2*pi*rm/ZS(k);
        cR(k) = sigmaR*sR(k);
        cS(k) = sigmaS*sS(k);
        CxR(k) = cR(k)*abs(cosd(stagR));
        CxS(k) = cS(k)*abs(cosd(stagS));
        tR(k) = tcR*cR(k);
        tS(k) = tcS*cS(k);
        ReR(k) = rho1*W1*cR(k)/max(mu1,1e-12);
        ReS(k) = rho2*V2*cS(k)/max(mu2,1e-12);
        AR_R(k) = span/max(cR(k),1e-12);
        AR_S(k) = span/max(cS(k),1e-12);
    end
end

T = table(ZR,ZS,sR,sS,cR,cS,CxR,CxS,tR,tS,ReR,ReS,AR_R,AR_S, ...
    'VariableNames',{'ZR','ZS','pitch_R','pitch_S','chord_R','chord_S', ...
                     'Cx_R','Cx_S','t_R','t_S','Re_R','Re_S','AR_R','AR_S'});
end


function phi_grid = augment_phi_grid_for_mdot_target(phi_grid, in)
% Add local phi samples near the phi implied by mdot for the project-fixed
% Mx/htr/RPM closure. In this closure, Mx fixes Vx and rho, while phi fixes
% rm = (Vx/phi)/omega. With fixed hub-to-tip ratio, A scales with rm^2, so
% mdot scales approximately with 1/phi^2. This helper only improves sampling;
% it does not force mdot.
if ~isfield(in,'design') || ~field_or(in.design,'phi_grid_include_mdot_target',false)
    return
end

try
    air = in.air;
    amb = isa_atmosphere(in.altitude_m);
    [T01,p01,~] = inlet_total_conditions(amb, in.M_flight, air, ...
        field_or(in.design,'inlet_total_mode','freestream_total'));
    [~,~,rho1,Vx] = inlet_static_from_Mx(T01,p01,in.Mx_design,air);
    omega = 2*pi*in.N_rpm/60;

    % Compute mdot at phi = 1 using the same radius/area closure as the
    % candidate evaluator, then use mdot(phi)=mdot(phi=1)/phi^2.
    rm_phi1 = Vx/omega;
    [rt1,rh1] = radii_from_mean_and_ratio(rm_phi1, in.hub_to_tip, in.design.mean_radius_rule);
    A_phi1 = pi*(rt1^2-rh1^2);
    mdot_phi1 = rho1*Vx*A_phi1;
    phi_target = sqrt(mdot_phi1/max(in.mdot,1e-12));

    phi_min = field_or(in.design,'phi_min',min(phi_grid));
    phi_max = field_or(in.design,'phi_max',max(phi_grid));
    if isfinite(phi_target) && phi_target >= phi_min && phi_target <= phi_max
        half_width = field_or(in.design,'phi_grid_mdot_refine_half_width',0.08);
        n_ref = max(3,round(field_or(in.design,'phi_grid_mdot_refine_points',9)));
        lo = max(phi_min, phi_target*(1-half_width));
        hi = min(phi_max, phi_target*(1+half_width));
        local_phi = linspace(lo, hi, n_ref);
        phi_grid = unique(sort([phi_grid(:); phi_target; local_phi(:)]).');
    end
catch
    % Do not let a sampling helper stop the design run. The original coarse
    % grid remains valid if the target estimate cannot be formed.
end
end


function [T0,p0,note] = inlet_total_conditions(amb,M0,air,mode)
switch lower(mode)
    case {'freestream_total','ram_total','flight_total'}
        T0 = amb.T .* (1 + 0.5*(air.gamma-1)*M0.^2);
        p0 = amb.p .* (1 + 0.5*(air.gamma-1)*M0.^2).^(air.gamma/(air.gamma-1));
        note = 'Inlet total conditions computed from cruise Mach and ISA static ambient conditions.';
    case {'ambient_static','static'}
        T0 = amb.T;
        p0 = amb.p;
        note = 'Inlet total conditions set equal to ISA static ambient conditions; cruise Mach is reported but no ram rise is applied.';
    otherwise
        error('Unknown design.inlet_total_mode: %s', mode)
end
end


function r = rel_dev_residual(actual,target,tol_rel)
r = abs(actual-target)/max(abs(target),1e-12) - tol_rel;
end


function p = hard_constraint_penalty(c)
c = c(:);
bad = c(isfinite(c) & c > 0);
if isempty(bad)
    p = 0;
else
    p = 1e3 + 1e4*sum(bad.^2) + 1e2*max(bad);
end
if any(~isfinite(c))
    p = p + 1e6;
end
end


function status_msg(in,tag,msg)
if isfield(in,'design') && field_or(in.design,'status_print',true)
    fprintf('[meanline][%s] %s\n', tag, msg);
end
end


function in = resolve_local_data_paths(in)
base_dir = fileparts(mfilename('fullpath'));

if isfield(in,'howell_paths')
    names = fieldnames(in.howell_paths);
    for i = 1:numel(names)
        f = names{i};
        in.howell_paths.(f) = resolve_one_path(in.howell_paths.(f), base_dir);
    end
end

if isfield(in,'lieblein_paths')
    names = fieldnames(in.lieblein_paths);
    for i = 1:numel(names)
        f = names{i};
        if strcmpi(f,'folder')
            continue
        end
        in.lieblein_paths.(f) = resolve_one_path(in.lieblein_paths.(f), base_dir);
    end
end
end


function p = resolve_one_path(p, base_dir)
if isempty(p) || ~(ischar(p) || isstring(p))
    return
end
p = char(p);
if isfile(p)
    return
end
[~,name,ext] = fileparts(p);
p2 = fullfile(base_dir,[name ext]);
if isfile(p2)
    p = p2;
end
end


function out_best = meanline_design_auto_beta_tt(in)
% Searches beta_tt within the allowed range and selects the design minimizing:
% J = 1 - eta_tt + w_beta*((beta_tt - beta_tt_preferred)/beta_tt_preferred)^2.
% Infeasible candidates receive a large penalty, but the least-bad candidate
% is still returned so diagnostics remain available.

beta_grid = linspace(in.design.beta_tt_min, ...
                     in.design.beta_tt_max, ...
                     in.design.beta_tt_search_points);

best_score = inf;
out_best = [];

history = struct();
history.beta_tt = [];
history.score = [];
history.eta = [];
history.feasible = [];
history.max_residual = [];
history.message = {};

found_feasible = false;

for ii = 1:numel(beta_grid)
    trial = in;
    trial.beta_tt_target = beta_grid(ii);
    trial.design.beta_tt_mode = 'fixed';
    trial.design.beta_tt_search_active = true;

    try
        cand = meanline_fan_design(trial);

        c = cand.constraints.residuals.all_ineq(:);
        max_res = max(c);
        feasible = all(isfinite(c)) && all(c <= 0);

        J = meanline_preliminary_objective(cand, in);

        if ~feasible
            J = 1e3 + J + 10*sum(max(c,0).^2);
        end

        history.beta_tt(end+1,1) = beta_grid(ii);
        history.score(end+1,1) = J;
        history.eta(end+1,1) = cand.performance.eta_tt_est;
        history.feasible(end+1,1) = feasible;
        history.max_residual(end+1,1) = max_res;
        history.message{end+1,1} = '';

        if feasible
            found_feasible = true;
        end

        if J < best_score
            best_score = J;
            out_best = cand;
        end

    catch ME
        history.beta_tt(end+1,1) = beta_grid(ii);
        history.score(end+1,1) = inf;
        history.eta(end+1,1) = NaN;
        history.feasible(end+1,1) = false;
        history.max_residual(end+1,1) = inf;
        history.message{end+1,1} = ME.message;
    end
end

if isempty(out_best)
    error('beta_tt search failed: no candidate could be evaluated.')
end

out_best.beta_tt_search = history;
out_best.beta_tt_search.selected_beta_tt = out_best.coeffs.beta_tt_est;
out_best.beta_tt_search.mode = in.design.beta_tt_mode;
out_best.beta_tt_search.found_feasible = found_feasible;
out_best.beta_tt_search.objective = meanline_preliminary_objective(out_best, in);
end


function J = meanline_preliminary_objective(out, in)
% Preliminary scalar objective used for beta_tt and R searches.

eta = out.performance.eta_tt_est;

if isfield(out,'coeffs') && isfield(out.coeffs,'beta_tt_est')
    beta_eval = out.coeffs.beta_tt_est;
elseif isfield(out,'performance') && isfield(out.performance,'beta_tt_est')
    beta_eval = out.performance.beta_tt_est;
else
    beta_eval = in.beta_tt_target;
end

if isfield(in.design,'beta_tt_preferred')
    beta_pref = in.design.beta_tt_preferred;
else
    beta_pref = in.beta_tt_target;
end

if isfield(in.design,'beta_tt_penalty_weight')
    w_beta = in.design.beta_tt_penalty_weight;
else
    w_beta = 0;
end

J = 1 - eta + w_beta*((beta_eval - beta_pref)/beta_pref)^2;

relax_mdot_penalty = (field_or(in.design,'relax_mdot_penalty_when_matching_tip_radius',false) && ...
    field_or(in.design,'match_Mrel_tip_by_mean_radius',false)) || ...
    (field_or(in.design,'relax_mdot_penalty_in_fixed_radius_duty',false) && ...
    (strcmpi(field_or(in.design,'solve_mode',''),'sweep_fixed_mean_radius_duty') || ...
     field_or(in.design,'fixed_mean_radius_duty_active',false)));

% Soft project-matching terms. These influence ranking but do not turn the
% candidate into a hard failure unless the corresponding hard switches are
% enabled in default_fan_inputs.
if ~relax_mdot_penalty && isfield(in.design,'mdot_penalty_weight') && in.design.mdot_penalty_weight > 0 && ...
        isfield(out,'design_basis') && isfield(out.design_basis,'mdot_rel_error')
    J = J + in.design.mdot_penalty_weight * out.design_basis.mdot_rel_error^2;
end

if isfield(in.design,'Mrel_tip_penalty_weight') && in.design.Mrel_tip_penalty_weight > 0 && ...
        isfield(out,'design_basis') && isfield(out.design_basis,'Mrel_tip_actual')
    Mrel_target = max(out.design_basis.Mrel_tip_target,1e-12);
    Mrel_err = (out.design_basis.Mrel_tip_actual - out.design_basis.Mrel_tip_target)/Mrel_target;
    J = J + in.design.Mrel_tip_penalty_weight * Mrel_err^2;
end

% Bound-excess penalties: no penalty inside the project tolerance band, then
% quadratic growth outside the band. This is useful when mdot/Mrel/beta are
% soft diagnostics but should still dominate ranking once they are clearly
% outside the allowed band.
if ~relax_mdot_penalty && isfield(in.design,'mdot_excess_penalty_weight') && in.design.mdot_excess_penalty_weight > 0 && ...
        isfield(out,'constraints') && isfield(out.constraints.residuals,'mdot_excess_rel')
    J = J + in.design.mdot_excess_penalty_weight * out.constraints.residuals.mdot_excess_rel^2;
end

if isfield(in.design,'Mrel_tip_excess_penalty_weight') && in.design.Mrel_tip_excess_penalty_weight > 0 && ...
        isfield(out,'constraints') && isfield(out.constraints.residuals,'Mrel_tip_excess_rel')
    J = J + in.design.Mrel_tip_excess_penalty_weight * out.constraints.residuals.Mrel_tip_excess_rel^2;
end

if isfield(in.design,'beta_tt_excess_penalty_weight') && in.design.beta_tt_excess_penalty_weight > 0 && ...
        isfield(out,'constraints') && isfield(out.constraints.residuals,'beta_tt_excess')
    J = J + in.design.beta_tt_excess_penalty_weight * out.constraints.residuals.beta_tt_excess^2;
end

if isfield(in.design,'R_penalty_weight') && in.design.R_penalty_weight > 0 && ...
        isfield(out,'coeffs') && isfield(out.coeffs,'R')
    R_pref = field_or(in.design,'R_preferred',0.5);
    R_span = max(field_or(in.design,'R_max',1) - field_or(in.design,'R_min',0),1e-12);
    J = J + in.design.R_penalty_weight * ((out.coeffs.R - R_pref)/R_span)^2;
end

if isfield(in.design,'RR_penalty_weight') && in.design.RR_penalty_weight > 0 && ...
        isfield(out,'stability') && isfield(out.stability,'RR_soft_min')
    J = J + in.design.RR_penalty_weight * max(0,out.stability.RR_soft_min - out.stability.RR)^2;
end
end


function out_best = meanline_design_auto_R(in)
R_grid = linspace(in.design.R_min, in.design.R_max, in.design.R_search_points);

best_score = inf;
out_best = [];

history = struct();
history.R = [];
history.score = [];
history.eta = [];
history.feasible = [];
history.max_residual = [];
history.message = {};

found_feasible = false;
for ii = 1:numel(R_grid)
    trial = in;
    trial.R = R_grid(ii);
    trial.design.R_mode = 'fixed';
    trial.design.R_search_active = true;

    try
        cand = meanline_fan_design(trial);

        c = cand.constraints.residuals.all_ineq(:);
        max_res = max(c);

        feasible = all(isfinite(c)) && all(c <= 0) && ...
                   cand.coeffs.phi_used >= in.design.phi_min && ...
                   cand.coeffs.phi_used <= in.design.phi_max && ...
                   cand.coeffs.psi >= in.design.psi_min && ...
                   cand.coeffs.psi <= in.design.psi_max;

        J = meanline_preliminary_objective(cand, in);

        if feasible
            score = J;
        else
            score = 1e3 + J + 10*sum(max(c,0).^2) + ...
                    10*max(0,in.design.phi_min-cand.coeffs.phi_used)^2 + ...
                    10*max(0,cand.coeffs.phi_used-in.design.phi_max)^2 + ...
                    10*max(0,in.design.psi_min-cand.coeffs.psi)^2 + ...
                    10*max(0,cand.coeffs.psi-in.design.psi_max)^2;
        end

        history.R(end+1,1) = R_grid(ii);
        history.score(end+1,1) = score;
        history.eta(end+1,1) = cand.performance.eta_tt_est;
        history.feasible(end+1,1) = feasible;
        history.max_residual(end+1,1) = max_res;
        history.message{end+1,1} = '';

        if feasible
            found_feasible = true;
        end

        if score < best_score
            best_score = score;
            out_best = cand;
        end

    catch ME
        history.R(end+1,1) = R_grid(ii);
        history.score(end+1,1) = inf;
        history.eta(end+1,1) = NaN;
        history.feasible(end+1,1) = false;
        history.max_residual(end+1,1) = inf;
        history.message{end+1,1} = ME.message;
    end
end

if isempty(out_best)
    error('R search failed: no candidate could be evaluated.')
end

out_best.R_search = history;
out_best.R_search.selected_R = out_best.coeffs.R;
out_best.R_search.mode = in.design.R_mode;
out_best.R_search.found_feasible = found_feasible;
end

function sol = size_solidity_from_DF_Howell(rowName, required_turn_deg, exit_angle_deg, Re, DF_base, DF_turn, howell, in)
target = in.design.DF_sizing_target;
DF_min = in.design.DF_min;
DF_max = in.design.DF_limit;

sigma_min = in.design.solidity_min;
sigma_max = in.design.solidity_max;
N = in.design.solidity_search_points;

den = target - DF_base;
if den > 0 && DF_turn > 0
    sigma_from_DF = DF_turn / den;
else
    sigma_from_DF = NaN;
end

sigma_grid = linspace(sigma_min, sigma_max, N);
if isfinite(sigma_from_DF)
    sigma_grid = unique(sort([sigma_grid, min(max(sigma_from_DF,sigma_min),sigma_max)]));
end
DF_grid = DF_base + DF_turn ./ sigma_grid;

f_beta = howell.f_beta2(exit_angle_deg);
Phi_Re = howell.Phi_Re(Re/1e5);
sc_grid = 1 ./ sigma_grid;

Howell_limit_grid = f_beta .* Phi_Re .* howell.Psi_sc(sc_grid);
Howell_residual_grid = required_turn_deg - Howell_limit_grid;

ok_DF = DF_grid >= DF_min & DF_grid <= DF_max;
ok_Howell = Howell_residual_grid <= 0;
ok = ok_DF & ok_Howell & isfinite(DF_grid) & isfinite(Howell_residual_grid);

selected_reason = '';
if any(ok)
    idx_ok = find(ok);

    switch lower(in.design.solidity_policy)
        case 'closest_to_df_target'
            [~,j] = min(abs(DF_grid(idx_ok) - target));
            idx = idx_ok(j);
            selected_reason = 'closest feasible DF to target';

        case 'minimum_solidity_feasible'
            [~,j] = min(sigma_grid(idx_ok));
            idx = idx_ok(j);
            selected_reason = 'minimum feasible solidity';

        case {'df_target_then_howell','df_target_then_howell_increment','df_target_howell'}
            if isfinite(sigma_from_DF)
                sigma_start = min(max(sigma_from_DF,sigma_min),sigma_max);
            else
                sigma_start = sigma_min;
            end

            % Start from the DF=target solidity and only increase solidity
            % if Howell turning requires it. This avoids the previous bias
            % toward the minimum feasible solidity.
            idx_after = idx_ok(sigma_grid(idx_ok) >= sigma_start - 1e-12);
            if ~isempty(idx_after)
                [~,j] = min(sigma_grid(idx_after));
                idx = idx_after(j);
                if abs(sigma_grid(idx) - sigma_start) <= 1e-6*max(1,sigma_start)
                    selected_reason = 'DF target satisfied Howell';
                else
                    selected_reason = 'DF target increased until Howell passed';
                end
            else
                % If increasing from the DF target cannot satisfy all checks
                % inside the search range, fall back to the feasible point
                % closest to the target DF instead of jumping to sigma_min.
                [~,j] = min(abs(DF_grid(idx_ok) - target));
                idx = idx_ok(j);
                selected_reason = 'fallback closest feasible DF to target';
            end

        otherwise
            error('Unknown design.solidity_policy: %s', in.design.solidity_policy)
    end

    status = 'feasible';
else
    violation = max(0, DF_min - DF_grid).^2 + ...
                max(0, DF_grid - DF_max).^2 + ...
                max(0, Howell_residual_grid/max(required_turn_deg,1)).^2;

    if isfinite(sigma_from_DF)
        % In infeasible cases, prefer the least-violating point near/above
        % the DF-target solidity rather than blindly selecting sigma_min.
        sigma_start = min(max(sigma_from_DF,sigma_min),sigma_max);
        violation = violation + 1e-3*((sigma_grid - sigma_start)./max(sigma_start,1e-12)).^2;
    end

    [~,idx] = min(violation);
    status = 'no_feasible_sigma_in_search_range';
    selected_reason = 'least violation in search range';
end

sol = struct();
sol.enabled = true;
sol.row = rowName;
sol.sigma_initial_from_DF_target = sigma_from_DF;
sol.sigma = sigma_grid(idx);
sol.s_over_c = 1/sol.sigma;
sol.DF = DF_grid(idx);
sol.DF_target = target;
sol.DF_min = DF_min;
sol.DF_limit = DF_max;
sol.required_turn_deg = required_turn_deg;
sol.howell_limit_deg = Howell_limit_grid(idx);
sol.howell_residual_deg = Howell_residual_grid(idx);
sol.howell_ok = Howell_residual_grid(idx) <= 0;
sol.DF_ok = DF_grid(idx) >= DF_min && DF_grid(idx) <= DF_max;
sol.status = status;
sol.selected_reason = selected_reason;
sol.policy = in.design.solidity_policy;
sol.Re = Re;
sol.exit_angle_deg = exit_angle_deg;
end

function [Cx, stagger_deg] = axial_chord_projection(chord, chi_in_signed, chi_out_signed, lieblein, rowType)
% Returns axial projection of true chord:
%   Cx = c |cos(stagger)|.
% Uses Lieblein signed stagger if available, otherwise signed mean flow angle.

stagger_deg = rad2deg(0.5*(chi_in_signed + chi_out_signed));

if isstruct(lieblein) && isfield(lieblein,'enabled') && lieblein.enabled && ...
        isfield(lieblein,rowType)
    row = lieblein.(rowType);
    if isfield(row,'stagger_deg') && isfinite(row.stagger_deg)
        stagger_deg = row.stagger_deg;
    elseif isfield(row,'stagger_signed_deg') && isfinite(row.stagger_signed_deg)
        stagger_deg = row.stagger_signed_deg;
    end
end

Cx = chord * abs(cosd(stagger_deg));
end


function [rt,rh,rm,phi,psi,Mrel_tip,Vt1_tip] = solve_tip_radius_from_mdot_tipmach(A, Vx, T1, air, omega, R, Delta_h0, Mrel_tip_target, mean_radius_rule, vtheta_mode)
%SOLVE_TIP_RADIUS_FROM_MDOT_TIPMACH Enforce mdot area and tip relative Mach.
% A and Vx are fixed by mass flow and axial Mach. For a trial tip radius,
% the hub radius follows from A = pi*(rt^2-rh^2), the mean radius gives
% phi and psi, and the tip relative inlet Mach can be evaluated. The root
% in rt is then solved so Mrel_tip equals the requested value.

a1 = sqrt(air.gamma*air.R*T1);
rt_min = sqrt(A/pi)*(1 + 1e-7);
rt_max = max(1.25*rt_min, 0.25);

    function [fval, vals] = eval_rt(rt_trial)
        vals = struct('valid',false);
        if rt_trial <= sqrt(A/pi)
            fval = NaN;
            return
        end
        rh_trial = sqrt(rt_trial^2 - A/pi);
        if ~isreal(rh_trial) || rh_trial <= 0 || rh_trial >= rt_trial
            fval = NaN;
            return
        end
        rm_trial = mean_radius_from_rule(rt_trial, rh_trial, mean_radius_rule);
        Umean = omega*rm_trial;
        Utip = omega*rt_trial;
        phi_trial = Vx/Umean;
        psi_trial = Delta_h0/Umean^2;
        tan_a1_trial = (1 - R - psi_trial/2)/phi_trial;
        Vt1_mean = Vx*tan_a1_trial;
        Vt1_tip_trial = spanwise_vtheta(Vt1_mean, rt_trial, rm_trial, vtheta_mode);
        Mrel = hypot(Vx, Vt1_tip_trial - Utip)/a1;
        fval = Mrel - Mrel_tip_target;
        vals = struct('valid',true,'rt',rt_trial,'rh',rh_trial,'rm',rm_trial, ...
            'phi',phi_trial,'psi',psi_trial,'Mrel_tip',Mrel,'Vt1_tip',Vt1_tip_trial);
    end

% Expand upper radius until the target is bracketed or the search becomes unreasonable.
[flo,~] = eval_rt(rt_min);
[fhi,~] = eval_rt(rt_max);
expand_count = 0;
while (isnan(flo) || isnan(fhi) || sign(flo)==sign(fhi)) && expand_count < 80
    rt_max = 1.20*rt_max + 0.02;
    [fhi,~] = eval_rt(rt_max);
    expand_count = expand_count + 1;
end

if ~(isfinite(flo) && isfinite(fhi) && sign(flo) ~= sign(fhi))
    % Fall back to a dense search to provide a useful error diagnostic.
    rt_grid = linspace(rt_min, rt_max, 1000);
    f_grid = nan(size(rt_grid));
    for ii = 1:numel(rt_grid)
        f_grid(ii) = eval_rt(rt_grid(ii));
    end
    [~,ii_best] = min(abs(f_grid),[],'omitnan');
    if isempty(ii_best) || ~isfinite(f_grid(ii_best))
        error('Could not evaluate any valid tip radius for the mdot/Mrel_tip sizing problem.');
    end
    error(['No tip radius root found for the requested mdot, Mx, Mrel_tip, PR and R. ', ...
           'Closest sampled Mrel_tip residual was %.4g at rt=%.4f m. ', ...
           'Try changing R bounds, eta_work_guess, Mrel_tip, Mx or mdot.'], f_grid(ii_best), rt_grid(ii_best));
end

rt = fzero(@(x) eval_rt(x), [rt_min rt_max]);
[~, vals] = eval_rt(rt);
rh = vals.rh;
rm = vals.rm;
phi = vals.phi;
psi = vals.psi;
Mrel_tip = vals.Mrel_tip;
Vt1_tip = vals.Vt1_tip;
end


function [rt,rh,rm,A,match] = match_radius_to_tip_mach(rm0, lambda, mean_radius_rule, Vx, T1, air, omega, psi, R_in, alpha1_mode, alpha1_deg, Mrel_target, vtheta_mode, design)
%MATCH_RADIUS_TO_TIP_MACH Rescale mean radius about the continuity value to
% match the target rotor-inlet relative tip Mach number. The hub-to-tip ratio
% remains fixed, so area and meanline mdot become diagnostics rather than the
% quantity used to set r_m. This is intended for MEANGEN/MULTALL starting
% values where the downstream simulation will prescribe mdot and adjust blade
% shape/throttle.
max_rel = field_or(design,'Mrel_tip_radius_match_max_rel',0.20);
tol_rel = field_or(design,'Mrel_tip_radius_match_tol_rel',field_or(design,'fixed_tolerance_rel',0.02));
rm_lo = max(rm0*(1-max_rel), 1e-6);
rm_hi = max(rm0*(1+max_rel), rm_lo*(1+1e-6));

a1 = sqrt(air.gamma*air.R*T1);
forced_alpha1 = any(strcmpi(alpha1_mode,{'zero','force_zero','fixed_zero'}));

    function [fval,Mrel,vals] = residual_for_rm(rm_trial)
        [rt_trial,rh_trial] = radii_from_mean_and_ratio(rm_trial, lambda, mean_radius_rule);
        U_trial = omega*rm_trial;
        phi_trial = Vx/max(U_trial,1e-12);
        if forced_alpha1
            tan_a1_trial = tand(alpha1_deg);
        else
            tan_a1_trial = (1 - R_in - psi/2)/max(phi_trial,1e-12);
        end
        Vt1_trial = Vx*tan_a1_trial;
        Vt1_tip_trial = spanwise_vtheta(Vt1_trial, rt_trial, rm_trial, vtheta_mode);
        Utip_trial = omega*rt_trial;
        Mrel = hypot(Vx, Vt1_tip_trial - Utip_trial)/max(a1,1e-12);
        fval = Mrel - Mrel_target;
        vals = struct('rt',rt_trial,'rh',rh_trial,'rm',rm_trial, ...
            'A',pi*(rt_trial^2-rh_trial^2), ...
            'phi',phi_trial,'Vt1_tip',Vt1_tip_trial,'Utip',Utip_trial, ...
            'Mrel_tip',Mrel);
    end

[f_lo,M_lo] = residual_for_rm(rm_lo);
[f_hi,M_hi] = residual_for_rm(rm_hi);
status = 'matched_by_fzero';

if isfinite(f_lo) && isfinite(f_hi) && f_lo*f_hi <= 0
    try
        rm = fzero(@(x) residual_for_rm(x), [rm_lo rm_hi]);
        [f_sel,M_sel,vals] = residual_for_rm(rm);
    catch
        rm_grid = linspace(rm_lo,rm_hi,81);
        f_grid = nan(size(rm_grid)); M_grid = nan(size(rm_grid)); vals_grid = cell(size(rm_grid));
        for ii = 1:numel(rm_grid)
            [f_grid(ii),M_grid(ii),vals_grid{ii}] = residual_for_rm(rm_grid(ii));
        end
        [~,ii] = min(abs(f_grid));
        vals = vals_grid{ii}; rm = vals.rm; f_sel = f_grid(ii); M_sel = M_grid(ii);
        status = 'closest_grid_after_fzero_failure';
    end
else
    rm_grid = linspace(rm_lo,rm_hi,81);
    f_grid = nan(size(rm_grid)); M_grid = nan(size(rm_grid)); vals_grid = cell(size(rm_grid));
    for ii = 1:numel(rm_grid)
        [f_grid(ii),M_grid(ii),vals_grid{ii}] = residual_for_rm(rm_grid(ii));
    end
    [~,ii] = min(abs(f_grid));
    vals = vals_grid{ii}; rm = vals.rm; f_sel = f_grid(ii); M_sel = M_grid(ii);
    if abs(f_sel)/max(Mrel_target,1e-12) <= tol_rel
        status = 'matched_by_grid_no_bracket';
    else
        status = 'closest_within_radius_bound_not_within_tolerance';
    end
end

rt = vals.rt;
rh = vals.rh;
A = vals.A;
match = struct();
match.enabled = true;
match.status = status;
match.rm_continuity = rm0;
[rt0,rh0] = radii_from_mean_and_ratio(rm0, lambda, mean_radius_rule);
match.rt_continuity = rt0;
match.rh_continuity = rh0;
match.A1_continuity = pi*(rt0^2-rh0^2);
match.rm_selected = rm;
match.rt_selected = rt;
match.rh_selected = rh;
match.A1_selected = A;
match.relative_shift = (rm-rm0)/max(rm0,1e-12);
match.max_relative_shift_allowed = max_rel;
match.Mrel_tip_target = Mrel_target;
match.Mrel_tip_selected = M_sel;
match.Mrel_tip_rel_error = f_sel/max(Mrel_target,1e-12);
match.within_tolerance = abs(match.Mrel_tip_rel_error) <= tol_rel;
match.tolerance_rel = tol_rel;
match.Mrel_tip_at_lower_bound = M_lo;
match.Mrel_tip_at_upper_bound = M_hi;
end

function d = make_continuity_radius_diagnostic(label, mdot, rho, Vx, hub_to_tip, mean_radius_rule)
%MAKE_CONTINUITY_RADIUS_DIAGNOSTIC Return area/radius implied by continuity.
A = mdot/(rho*Vx);
[rt,rh] = radii_from_area_and_ratio(A, hub_to_tip);
rm = mean_radius_from_rule(rt, rh, mean_radius_rule);
d = struct('label',label,'mdot',mdot,'rho',rho,'Vx',Vx,'A',A, ...
    'rt',rt,'rh',rh,'rm',rm,'hub_to_tip',hub_to_tip, ...
    'mean_radius_rule',mean_radius_rule);
end

function [T,p,rho] = inlet_static_from_Vx(T0,p0,Vx,air)
T = T0 - Vx^2/(2*air.cp);
if T <= 1
    error('Inlet static temperature became non-physical.');
end
p = p0 / (T0/T)^(air.gamma/(air.gamma-1));
rho = p/(air.R*T);
end

function [T,p,rho,Vx] = inlet_static_from_Mx(T0,p0,Mx,air)
% Convert a prescribed axial Mach number Mx = Vx/a_static to static inlet
% conditions and axial velocity at the fan face.
T = T0/(1 + 0.5*(air.gamma-1)*Mx^2);
p = p0/(1 + 0.5*(air.gamma-1)*Mx^2)^(air.gamma/(air.gamma-1));
rho = p/(air.R*T);
Vx = Mx*sqrt(air.gamma*air.R*T);
end

function f = mean_radius_fraction(lambda, rule)
% Return r_m/r_t for rh=lambda*rt.
switch lower(rule)
    case 'arithmetic'
        f = 0.5*(1 + lambda);
    case 'area'
        f = sqrt(0.5*(1 + lambda^2));
    otherwise
        error('Unknown design.mean_radius_rule: %s', rule)
end
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


function tmax_over_c = row_tmax_over_c(in, rowType)
switch lower(rowType)
    case 'rotor'
        if isfield(in,'lieblein') && isfield(in.lieblein,'rotor_tmax_over_c') && ...
                isfinite(in.lieblein.rotor_tmax_over_c)
            tmax_over_c = in.lieblein.rotor_tmax_over_c;
        elseif isfield(in,'rotor_thickness') && isfield(in,'rotor_chord') && in.rotor_chord > 0
            tmax_over_c = in.rotor_thickness/in.rotor_chord;
        elseif isfield(in,'lieblein') && isfield(in.lieblein,'default_tmax_over_c')
            tmax_over_c = in.lieblein.default_tmax_over_c;
        else
            tmax_over_c = 0.07;
        end
    case 'stator'
        if isfield(in,'lieblein') && isfield(in.lieblein,'stator_tmax_over_c') && ...
                isfinite(in.lieblein.stator_tmax_over_c)
            tmax_over_c = in.lieblein.stator_tmax_over_c;
        elseif isfield(in,'stator_thickness') && isfield(in,'stator_chord') && in.stator_chord > 0
            tmax_over_c = in.stator_thickness/in.stator_chord;
        elseif isfield(in,'lieblein') && isfield(in.lieblein,'default_tmax_over_c')
            tmax_over_c = in.lieblein.default_tmax_over_c;
        else
            tmax_over_c = 0.07;
        end
    otherwise
        error('Unknown rowType for row_tmax_over_c: %s', rowType)
end
end

function lieblein = add_signed_lieblein_geometry(lieblein, beta1, beta2, alpha2, alpha3)
if ~isfield(lieblein,'enabled') || ~lieblein.enabled
    return
end
if isfield(lieblein,'rotor') && isfield(lieblein.rotor,'valid') && lieblein.rotor.valid
    lieblein.rotor = add_signed_row_geometry(lieblein.rotor, beta1, beta2);
end
if isfield(lieblein,'stator') && isfield(lieblein.stator,'valid') && lieblein.stator.valid
    lieblein.stator = add_signed_row_geometry(lieblein.stator, alpha2, alpha3);
end
end

function row = add_signed_row_geometry(row, chi_in_signed, chi_out_signed)
% Preserve the positive-angle chart values for traceability.
if ~isfield(row,'metal_in_magnitude_deg')
    row.metal_in_magnitude_deg = row.metal_in_deg;
end
if ~isfield(row,'metal_out_magnitude_deg')
    row.metal_out_magnitude_deg = row.metal_out_deg;
end
if ~isfield(row,'stagger_magnitude_deg')
    row.stagger_magnitude_deg = row.stagger_deg;
end

inc_signed = signed_chart_angle(chi_in_signed, chi_in_signed, deg2rad(row.incidence_deg));
dev_signed = signed_chart_angle(chi_out_signed, chi_in_signed, deg2rad(row.deviation_deg));
metal_in_signed = chi_in_signed - inc_signed;
metal_out_signed = chi_out_signed - dev_signed;
stagger_signed = 0.5*(metal_in_signed + metal_out_signed);

row.incidence_signed_deg = rad2deg(inc_signed);
row.deviation_signed_deg = rad2deg(dev_signed);
row.metal_in_signed_deg = rad2deg(metal_in_signed);
row.metal_out_signed_deg = rad2deg(metal_out_signed);
row.stagger_signed_deg = rad2deg(stagger_signed);

% Make the primary geometry fields signed, because these are the quantities
% that should be reported and used as blade-geometry estimates.
row.metal_in_deg = row.metal_in_signed_deg;
row.metal_out_deg = row.metal_out_signed_deg;
row.stagger_deg = row.stagger_signed_deg;
end

function in_loss = apply_lieblein_to_loss_model(in, lieblein)
in_loss = in;
if ~isfield(in_loss,'loss_model') || ...
        ~field_or(in_loss.loss_model,'auto_use_lieblein_stagger',false)
    return
end
if ~isfield(lieblein,'enabled') || ~lieblein.enabled
    return
end
if ~isfield(lieblein,'rotor') || ~isfield(lieblein,'stator')
    return
end
if ~isfield(lieblein.rotor,'valid') || ~lieblein.rotor.valid || ...
   ~isfield(lieblein.stator,'valid') || ~lieblein.stator.valid
    return
end

contexts = {'tip','endwall','shock'};
for i = 1:numel(contexts)
    ctx = contexts{i};
    if isfield(in_loss.loss_model,ctx)
        autoCtx = field_or(in_loss.loss_model.(ctx),'auto_use_lieblein_stagger',true);
        if autoCtx
            in_loss.loss_model.(ctx).stagger_mode = 'lieblein_metal_angles';
            in_loss.loss_model.(ctx).rotor_incidence_deg = lieblein.rotor.incidence_deg;
            in_loss.loss_model.(ctx).rotor_deviation_deg = lieblein.rotor.deviation_deg;
            in_loss.loss_model.(ctx).stator_incidence_deg = lieblein.stator.incidence_deg;
            in_loss.loss_model.(ctx).stator_deviation_deg = lieblein.stator.deviation_deg;
        end
    end
end
in_loss.loss_model.lieblein_stagger_auto_applied = true;
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
    'section_names',{{}},'section_weights',[],'M1',[],'Mout',[], ...
    'Mss_max',[],'Mrep',[],'Mtarget',[],'beta1_deg',[], ...
    'tLE_over_pitch',[],'ds_le_over_R',[],'ds_passage_over_R',[], ...
    'ds_over_R',[],'h_loss',[],'zeta_sections',[], ...
    'passage_p02_over_p01',[],'passage_Mn1',[],'passage_M2_actual',[], ...
    'status',{{}},'Mss_mode','none','velocity_scale',max(v_out,1e-8));
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

        case {'koch_smith','koch_smith_le_and_passage','koch_smith_shock'}
            [zSW, shockDetails] = koch_smith_shock_loss(rowType, lm.shock, air, radial, v_out);

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
                    tip_surface_velocity_ratios(rowType, lm.tip, chi_in, chi_out, ...
                    radial.chi_in_signed(2), radial.chi_out_signed(2), sigma);

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

function [uSS,uPS,x_over_Cs,stagger,modeUsed] = tip_surface_velocity_ratios(rowType, tip, chi_in, chi_out, chi_in_signed, chi_out_signed, sigma)
modeUsed = lower(field_or(tip,'surface_velocity_mode','approx'));
switch modeUsed
    case 'approx'
        % Hall Appendix-B option-A approximation:
        %   uSS - uPS ~= Delta(tan chi)/sigma
        %   uSS + uPS ~= 2/cos(stagger)
        % where u = V_surface/Vx.  For the rotor, Delta(tan chi)/sigma is
        % equivalent to psi/(sigma*phi) when chi are row-relative angles.
        stagger = tip_stagger_angle(rowType, tip, chi_in_signed, chi_out_signed, 'mean');
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
        stagger = tip_stagger_angle(rowType, tip, chi_in_signed, chi_out_signed, 'mean');

    otherwise
        error('Unknown tip.surface_velocity_mode option: %s', modeUsed)
end
end

function stagger = tip_stagger_angle(rowType, tip, chi_in_signed, chi_out_signed, sectionName)
% Return signed blade stagger for the tip-leakage surface-velocity model.
% Stagger is a blade-geometry angle, not the average of absolute flow-angle
% magnitudes. The preferred preliminary estimate uses signed metal angles:
%   stagger = 0.5*((chi_in - i) + (chi_out - delta)).
stagger = generic_signed_stagger_angle(rowType, tip, chi_in_signed, chi_out_signed, sectionName, 'tip');
end

function stagger = generic_signed_stagger_angle(rowType, settings, chi_in_signed, chi_out_signed, sectionName, contextName)
%GENERIC_SIGNED_STAGGER_ANGLE Geometry-consistent signed stagger estimate.
%
% Supported settings.stagger_mode values:
%   'signed_flow_angles'       : 0.5*(chi_in + chi_out), signed fallback.
%   'lieblein_metal_angles'    : uses incidence/deviation from Lieblein-style
%                                positive-angle convention and converts them
%                                to signed row angles.
%   'signed_metal_angles'      : uses already signed incidence/deviation.
%   'user' / 'user_sections'   : reads prescribed signed stagger angle(s).
%   'magnitude_mean_flow_angles' : old absolute-angle fallback.
%
% For the Lieblein option, incidence and deviation are interpreted as chart
% values based on positive inlet/outlet angle magnitudes. They are converted
% to signed corrections before applying:
%   metal_in  = chi_in  - i_signed
%   metal_out = chi_out - delta_signed
%   stagger   = 0.5*(metal_in + metal_out)
mode = lower(field_or(settings,'stagger_mode','signed_flow_angles'));
row = lower(rowType);
sectionName = lower(sectionName);

switch mode
    case {'signed_flow_angles','signed_mean_flow_angles','mean_flow_angles'}
        % New default: signed average of inlet/outlet flow angles. The name
        % 'mean_flow_angles' is retained for backward compatibility, but it
        % no longer uses absolute values.
        stagger = 0.5*(chi_in_signed + chi_out_signed);

    case {'magnitude_mean_flow_angles','abs_mean_flow_angles','absolute_mean_flow_angles'}
        % Legacy fallback only. This is not a true geometry angle.
        stagger = 0.5*(abs(chi_in_signed) + abs(chi_out_signed));

    case {'lieblein_metal_angles','metal_angles_from_lieblein'}
        [inc_deg, dev_deg] = incidence_deviation_from_settings(settings, row, sectionName, contextName, false);
        inc_signed = signed_chart_angle(chi_in_signed, chi_in_signed, deg2rad(inc_deg));
        dev_signed = signed_chart_angle(chi_out_signed, chi_in_signed, deg2rad(dev_deg));
        metal_in = chi_in_signed - inc_signed;
        metal_out = chi_out_signed - dev_signed;
        stagger = 0.5*(metal_in + metal_out);

    case {'signed_metal_angles','metal_angles_signed'}
        [inc_deg, dev_deg] = incidence_deviation_from_settings(settings, row, sectionName, contextName, true);
        metal_in = chi_in_signed - deg2rad(inc_deg);
        metal_out = chi_out_signed - deg2rad(dev_deg);
        stagger = 0.5*(metal_in + metal_out);

    case {'user','user_sections','sectional'}
        stagger_deg = signed_stagger_from_settings(settings, row, sectionName, contextName);
        stagger = deg2rad(stagger_deg);

    otherwise
        error('Unknown %s.stagger_mode option: %s', contextName, mode)
end

stagger = clamp_signed_stagger(stagger);
end

function [inc_deg, dev_deg] = incidence_deviation_from_settings(settings, row, sectionName, contextName, allowDefault)
inc_deg = row_section_field(settings, row, 'incidence', sectionName, NaN);
dev_deg = row_section_field(settings, row, 'deviation', sectionName, NaN);

if ~isfinite(inc_deg)
    inc_deg = field_or(settings,[row '_incidence_deg'],NaN);
end
if ~isfinite(dev_deg)
    dev_deg = field_or(settings,[row '_deviation_deg'],NaN);
end
if ~isfinite(inc_deg)
    inc_deg = field_or(settings,'incidence_deg',NaN);
end
if ~isfinite(dev_deg)
    dev_deg = field_or(settings,'deviation_deg',NaN);
end

if (~isfinite(inc_deg) || ~isfinite(dev_deg)) && ~allowDefault
    error('%s.stagger_mode = lieblein_metal_angles requires %s_incidence_deg and %s_deviation_deg, or section-specific fields.', ...
        contextName, row, row)
end
if ~isfinite(inc_deg); inc_deg = 0; end
if ~isfinite(dev_deg); dev_deg = 0; end
end

function stagger_deg = signed_stagger_from_settings(settings, row, sectionName, contextName)
stagger_deg = row_section_field(settings, row, 'stagger', sectionName, NaN);
if ~isfinite(stagger_deg)
    stagger_deg = field_or(settings,[row '_stagger_deg'],NaN);
end
if ~isfinite(stagger_deg)
    stagger_deg = field_or(settings,'stagger_deg',NaN);
end
if ~isfinite(stagger_deg)
    error('%s.stagger_mode = user requires row-specific or section-specific signed stagger angle.', contextName)
end
end

function val = row_section_field(settings, row, quantity, sectionName, defaultVal)
% Accept both rotor_stagger_tip_deg and rotor_tip_stagger_deg naming styles.
val = field_or(settings,[row '_' quantity '_' sectionName '_deg'],NaN);
if ~isfinite(val)
    val = field_or(settings,[row '_' sectionName '_' quantity '_deg'],NaN);
end
if ~isfinite(val)
    val = defaultVal;
end
end

function a_signed = signed_chart_angle(primaryAngle, fallbackAngle, a_chart)
% Convert a positive-angle chart correction to the signed row convention.
sgn = sign(primaryAngle);
if sgn == 0
    sgn = sign(fallbackAngle);
end
if sgn == 0
    sgn = 1;
end
a_signed = sgn*a_chart;
end

function stagger = clamp_signed_stagger(stagger)
% Preserve sign but avoid exactly +/-90 deg in cos(stagger) denominators.
lim = deg2rad(89.9);
if abs(stagger) > lim
    stagger = sign(stagger)*lim;
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


function [uSS,uPS,x_over_Cs,stagger,modeUsed] = endwall_surface_velocity_ratios(rowType, lm, chi_in, chi_out, chi_in_signed, chi_out_signed, sigma)
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
        stagger = endwall_stagger_angle(rowType, ew, chi_in_signed, chi_out_signed, 'mean');
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
        stagger = endwall_stagger_angle(rowType, ew, chi_in_signed, chi_out_signed, 'mean');

    case 'tip_shared'
        % Convenience mode: reuse the tip model's surface-velocity settings.
        % This is useful while the same MULTALL/StageN extraction is used
        % for both tip leakage and endwall BL loss.
        if ~isfield(lm,'tip')
            error('endwall.surface_velocity_mode = tip_shared requires lm.tip to exist.')
        end
        [uSS,uPS,x_over_Cs,stagger,tipMode] = ...
            tip_surface_velocity_ratios(rowType, lm.tip, chi_in, chi_out, chi_in_signed, chi_out_signed, sigma);
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

function stagger = endwall_stagger_angle(rowType, ew, chi_in_signed, chi_out_signed, sectionName)
% Return signed blade stagger for the endwall surface-velocity model.
stagger = generic_signed_stagger_angle(rowType, ew, chi_in_signed, chi_out_signed, sectionName, 'endwall');
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

function [zSW, details] = koch_smith_shock_loss(rowType, shock, air, radial, v_out_mean)
% Koch & Smith shock loss model for axial-flow compressors.
% Separates shock loss into leading-edge bluntness entropy rise and passage
% shock entropy rise. The passage shock is represented as one oblique shock
% that reduces a representative passage inlet Mach number to unity, or to
% the exit Mach number if the exit Mach number is supersonic.
[idx,weights,mode] = select_section_indices(shock, rowType, 'tip');
if strcmpi(rowType,'stator') && ~field_or(shock,'enabled_stator',false)
    zSW = 0;
    details = empty_koch_smith_details('koch_smith_le_and_passage',mode,max(v_out_mean,1e-8));
    return
end
if strcmpi(rowType,'rotor') && ~field_or(shock,'enabled_rotor',true)
    zSW = 0;
    details = empty_koch_smith_details('koch_smith_le_and_passage',mode,max(v_out_mean,1e-8));
    return
end

n = numel(idx);
hLoss = zeros(1,n); zeta_i = zeros(1,n);
M1v = zeros(1,n); Moutv = zeros(1,n); Mssv = NaN(1,n); Mrepv = NaN(1,n); Mtargetv = NaN(1,n);
betaDeg = zeros(1,n); tlePitch = zeros(1,n);
dsLE = zeros(1,n); dsPass = zeros(1,n); dsTot = zeros(1,n);
p02p01Pass = ones(1,n); Mn1Pass = NaN(1,n); M2Pass = NaN(1,n);
status = cell(1,n); sectionNames = cell(1,n);
MssModeUsed = cell(1,n);

useLE = field_or(shock,'use_le_bluntness',true);
usePassage = field_or(shock,'use_passage_shock',true);
Mmin = field_or(shock,'M_crit',field_or(shock,'M1_rel_min',1.0));
ksWeight = field_or(shock,'Mss_weight',6.0);

for k = 1:n
    j = idx(k);
    sectionNames{k} = radial.names{j};
    M1 = radial.Min(j);
    Mout = radial.Mout(j);
    beta1 = radial.chi_in(j);
    tLE_over_pitch = koch_smith_tLE_over_pitch(rowType, shock, radial, j);

    M1v(k) = M1;
    Moutv(k) = Mout;
    betaDeg(k) = rad2deg(beta1);
    tlePitch(k) = tLE_over_pitch;
    status{k} = 'no_shock';

    if useLE && M1 > Mmin
        machTerm = 1.28*(M1 - 1) + 0.96*(M1 - 1)^2;
        arg = 1 - (tLE_over_pitch/max(cos(beta1),1e-8))*machTerm;
        if arg <= 0
            arg = eps;
            status{k} = 'le_log_argument_clipped';
        else
            status{k} = 'le_bluntness';
        end
        dsLE(k) = max(-log(arg),0);
    end

    if usePassage
        [Mss, modeUsed] = koch_smith_surface_mach(rowType, shock, air, radial, j);
        Mssv(k) = Mss;
        MssModeUsed{k} = modeUsed;
        if isfinite(Mss)
            Mrep = (ksWeight*Mss + M1)/max(ksWeight + 1,1e-8);
        else
            Mrep = M1;
        end
        Mrepv(k) = Mrep;

        if Mout > 1
            Mtarget = Mout;
        else
            Mtarget = 1.0;
        end
        Mtargetv(k) = Mtarget;

        if Mrep > max(Mtarget,Mmin)
            [dsR,p02p01,Mn1,M2actual,ok] = oblique_shock_entropy_to_target(Mrep, Mtarget, air.gamma);
            dsPass(k) = max(dsR,0);
            p02p01Pass(k) = p02p01;
            Mn1Pass(k) = Mn1;
            M2Pass(k) = M2actual;
            if ok
                if strcmp(status{k},'no_shock')
                    status{k} = 'passage_shock';
                else
                    status{k} = [status{k} '+passage_shock'];
                end
            else
                if strcmp(status{k},'no_shock')
                    status{k} = 'passage_normal_shock_fallback';
                else
                    status{k} = [status{k} '+passage_normal_shock_fallback'];
                end
            end
        end
    else
        Mssv(k) = NaN;
        Mrepv(k) = M1;
        Mtargetv(k) = max(1.0,Mout);
        MssModeUsed{k} = 'disabled';
    end

    dsTot(k) = dsLE(k) + dsPass(k);
    if dsTot(k) > 0
        Tref = koch_smith_temperature_reference(shock, radial, j, Mssv(k), air);
        hLoss(k) = Tref * air.R * dsTot(k);
        zeta_i(k) = 2*hLoss(k)/max(v_out_mean,1e-8)^2;
    end
end

hLossWeighted = sum(weights .* hLoss);
zSW = 2*hLossWeighted/max(v_out_mean,1e-8)^2;

if n == 1
    MssMode = MssModeUsed{1};
else
    MssMode = strjoin(MssModeUsed, ',');
end

details = struct('model','koch_smith_le_and_passage','evaluate_at',mode, ...
    'section_names',{sectionNames},'section_weights',weights, ...
    'M1',M1v,'Mout',Moutv,'Mss_max',Mssv,'Mrep',Mrepv,'Mtarget',Mtargetv, ...
    'beta1_deg',betaDeg,'tLE_over_pitch',tlePitch, ...
    'ds_le_over_R',dsLE,'ds_passage_over_R',dsPass,'ds_over_R',dsTot, ...
    'h_loss',hLoss,'zeta_sections',zeta_i, ...
    'passage_p02_over_p01',p02p01Pass,'passage_Mn1',Mn1Pass, ...
    'passage_M2_actual',M2Pass,'status',{status}, ...
    'Mss_mode',MssMode,'velocity_scale',max(v_out_mean,1e-8));
end

function details = empty_koch_smith_details(model,mode,vscale)
details = struct('model',model,'evaluate_at',mode, ...
    'section_names',{{}},'section_weights',[],'M1',[],'Mout',[], ...
    'Mss_max',[],'Mrep',[],'Mtarget',[],'beta1_deg',[], ...
    'tLE_over_pitch',[],'ds_le_over_R',[],'ds_passage_over_R',[], ...
    'ds_over_R',[],'h_loss',[],'zeta_sections',[], ...
    'passage_p02_over_p01',[],'passage_Mn1',[],'passage_M2_actual',[], ...
    'status',{{}},'Mss_mode','none','velocity_scale',vscale);
end

function tLE_over_pitch = koch_smith_tLE_over_pitch(rowType, shock, radial, idx)
names = {'hub','mean','tip'};
suffix = names{idx};
row = lower(rowType);
tLE_over_pitch = field_or(shock,[row '_tLE_over_pitch_' suffix],NaN);
if ~isfinite(tLE_over_pitch)
    tLE_over_pitch = field_or(shock,[row '_tLE_over_pitch'],NaN);
end
if ~isfinite(tLE_over_pitch)
    tLE_over_pitch = field_or(shock,'tLE_over_pitch',NaN);
end
if ~isfinite(tLE_over_pitch)
    tLE_over_chord = field_or(shock,[row '_tLE_over_chord_' suffix],NaN);
    if ~isfinite(tLE_over_chord)
        tLE_over_chord = field_or(shock,[row '_tLE_over_chord'],field_or(shock,'tLE_over_chord',0.008));
    end
    tLE_over_pitch = tLE_over_chord * radial.chord(idx)/max(radial.pitch(idx),1e-8);
end
tLE_over_pitch = max(tLE_over_pitch,0);
end

function [Mss, modeUsed] = koch_smith_surface_mach(rowType, shock, air, radial, idx)
mode = lower(field_or(shock,'Mss_mode',field_or(shock,'surface_mach_mode','estimated')));
row = lower(rowType);
suffix = radial.names{idx};
M1 = radial.Min(idx);
modeUsed = mode;

switch mode
    case {'direct','user','direct_sections','user_sections'}
        Mss = field_or(shock,[row '_Mss_max_' suffix],NaN);
        if ~isfinite(Mss)
            Mss = field_or(shock,[row '_Mss_max'],field_or(shock,'Mss_max',NaN));
        end
        if ~isfinite(Mss)
            [Mss, modeUsed] = koch_smith_surface_mach_estimated(rowType, shock, air, radial, idx);
            modeUsed = ['direct_missing_' modeUsed];
        end

    case {'factor','multiplier'}
        fac = field_or(shock,[row '_Mss_factor'],field_or(shock,'Mss_factor',1.0));
        Mss = fac*M1;

    case {'estimated','estimate','surface_velocity_estimate','velocity_split'}
        [Mss, modeUsed] = koch_smith_surface_mach_estimated(rowType, shock, air, radial, idx);

    otherwise
        error('Unknown shock.Mss_mode option: %s', mode)
end
Mss = max(Mss,0);
end

function [Mss, modeUsed] = koch_smith_surface_mach_estimated(rowType, shock, air, radial, idx)
chi_in = radial.chi_in(idx);
chi_out = radial.chi_out(idx);
sigma = radial.chord(idx)/max(radial.pitch(idx),1e-8);
stagger = shock_stagger_angle(rowType, shock, radial, idx);
uDiff = abs(tan(chi_in) - tan(chi_out)) / max(sigma,1e-8);
uSum  = 2 / max(cos(stagger),1e-8);
uSS = max(0.5*(uSum + uDiff),0);
Wss = uSS * max(abs(radial.Vx(idx)),1e-8);
T0rel = radial.Tin(idx) + radial.v_in(idx)^2/(2*air.cp);
Tss = T0rel - Wss^2/(2*air.cp);
if Tss <= 1
    Tss = 1;
end
Mss = Wss/sqrt(air.gamma*air.R*Tss);
modeUsed = 'estimated_velocity_split';
end

function stagger = shock_stagger_angle(rowType, shock, radial, idx)
% Return signed blade stagger for the Koch--Smith Mss estimate.
% Section-resolved hub/mean/tip signed flow angles are used.
sectionName = radial.names{idx};
stagger = generic_signed_stagger_angle(rowType, shock, ...
    radial.chi_in_signed(idx), radial.chi_out_signed(idx), sectionName, 'shock');
end

function Tref = koch_smith_temperature_reference(shock, radial, idx, Mss, air)
mode = lower(field_or(shock,'temperature_reference','inlet_static'));
switch mode
    case 'inlet_static'
        Tref = radial.Tin(idx);
    case 'mean_static'
        Tref = 0.5*(radial.Tin(idx) + radial.Tout(idx));
    case 'surface_static'
        if isfinite(Mss) && Mss > 0
            T0rel = radial.Tin(idx) + radial.v_in(idx)^2/(2*air.cp);
            Tref = T0rel/(1 + 0.5*(air.gamma-1)*Mss^2);
        else
            Tref = radial.Tin(idx);
        end
    otherwise
        error('Unknown shock.temperature_reference option: %s', mode)
end
Tref = max(Tref,1);
end

function [dsR,p02_p01,Mn1,M2actual,ok] = oblique_shock_entropy_to_target(M1, Mtarget, gamma)
ok = true;
dsR = 0; p02_p01 = 1; Mn1 = NaN; M2actual = M1;
if M1 <= 1 || Mtarget >= M1
    return
end

[M2normal,p02normal] = normal_shock_downstream_total_mach_and_p02(M1, gamma);
if Mtarget <= M2normal
    Mn1 = M1;
    M2actual = M2normal;
    p02_p01 = p02normal;
    dsR = -log(max(p02_p01,eps));
    ok = false;
    return
end

f = @(mn) oblique_downstream_mach(M1, mn, gamma) - Mtarget;
lo = 1 + 1e-8;
hi = M1;
try
    Mn1 = fzero(f, [lo hi]);
catch
    a = lo; b = hi;
    fa = f(a); fb = f(b);
    if fa*fb > 0
        Mn1 = M1;
        ok = false;
    else
        for it = 1:80
            c = 0.5*(a+b);
            fc = f(c);
            if abs(fc) < 1e-10
                a = c; b = c; break
            end
            if fa*fc <= 0
                b = c; fb = fc;
            else
                a = c; fa = fc;
            end
        end
        Mn1 = 0.5*(a+b);
    end
end
M2actual = oblique_downstream_mach(M1, Mn1, gamma);
p02_p01 = normal_shock_total_pressure_ratio(Mn1, gamma);
dsR = -log(max(p02_p01,eps));
end

function M2 = oblique_downstream_mach(M1, Mn1, gamma)
Mn1 = min(max(Mn1,1+1e-10),M1);
Vt1_over_a1 = sqrt(max(M1^2 - Mn1^2,0));
Mn2 = normal_shock_downstream_mach_from_Mn(Mn1, gamma);
[~,~,T2_T1] = normal_shock_static_ratios(Mn1, gamma);
M2 = sqrt(Mn2^2 + Vt1_over_a1^2/max(T2_T1,1e-12));
end

function [M2,p02_p01] = normal_shock_downstream_total_mach_and_p02(M1, gamma)
M2 = normal_shock_downstream_mach_from_Mn(M1, gamma);
p02_p01 = normal_shock_total_pressure_ratio(M1, gamma);
end

function M2 = normal_shock_downstream_mach_from_Mn(Mn1, gamma)
M2 = sqrt((1 + 0.5*(gamma-1)*Mn1^2)/(gamma*Mn1^2 - 0.5*(gamma-1)));
end

function p02_p01 = normal_shock_total_pressure_ratio(Mn1, gamma)
[p2_p1,~,~] = normal_shock_static_ratios(Mn1, gamma);
Mn2 = normal_shock_downstream_mach_from_Mn(Mn1, gamma);
A1 = 1 + 0.5*(gamma-1)*Mn1^2;
A2 = 1 + 0.5*(gamma-1)*Mn2^2;
p02_p01 = p2_p1 * (A2/A1)^(gamma/(gamma-1));
p02_p01 = min(max(p02_p01,eps),1.0);
end

function [p2_p1,rho2_rho1,T2_T1] = normal_shock_static_ratios(Mn1, gamma)
p2_p1 = 1 + 2*gamma/(gamma+1)*(Mn1^2 - 1);
rho2_rho1 = ((gamma+1)*Mn1^2)/(2 + (gamma-1)*Mn1^2);
T2_T1 = p2_p1/rho2_rho1;
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
        endwall_surface_velocity_ratios(rowType, lm, radial.chi_in(j), radial.chi_out(j), ...
        radial.chi_in_signed(j), radial.chi_out_signed(j), sigma);
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

function eta_penalty = loss_efficiency_penalty_breakdown(loss, W2, V3, Delta_h0)
%LOSS_EFFICIENCY_PENALTY_BREAKDOWN Convert zeta source terms to eta penalties.
% Rotor loss is scaled with W2^2 and stator loss with V3^2, matching the
% meanline stage-efficiency estimate used in the main routine.
sources = {'profile','trailing_edge','shock','tip','endwall'};
eta_penalty = struct();
eta_penalty.rotor = struct();
eta_penalty.stator = struct();
eta_penalty.total = 0;
Delta_safe = max(Delta_h0,eps);
for k = 1:numel(sources)
    src = sources{k};
    zr = 0;
    zs = 0;
    if isfield(loss,'breakdown') && isfield(loss.breakdown,'rotor') && isfield(loss.breakdown.rotor,src)
        zr = loss.breakdown.rotor.(src);
    end
    if isfield(loss,'breakdown') && isfield(loss.breakdown,'stator') && isfield(loss.breakdown.stator,src)
        zs = loss.breakdown.stator.(src);
    end
    eta_penalty.rotor.(src) = 0.5*zr*W2^2/Delta_safe;
    eta_penalty.stator.(src) = 0.5*zs*V3^2/Delta_safe;
    eta_penalty.total = eta_penalty.total + eta_penalty.rotor.(src) + eta_penalty.stator.(src);
end
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


function RR = recovery_ratio_from_angles(alpha1, beta2, phi)
%RECOVERY_RATIO_FROM_ANGLES Cumpsty/course compressor recovery-ratio expression.
% alpha1 and beta2 must be provided using the derivation's convention. In the
% default code path alpha1 is signed absolute inlet swirl, while beta2 is the
% positive geometric rotor-exit relative angle given by beta2_RR=-beta2_code.
phi_safe = max(phi,1e-12);
RR = (cos(beta2)^2 * tan(beta2))/phi_safe ...
   - (cos(alpha1)^2 * cos(beta2)^2 * tan(alpha1)*tan(beta2))/max(phi_safe^2,1e-12) ...
   + (cos(alpha1)^2 * tan(alpha1))/phi_safe;
end

