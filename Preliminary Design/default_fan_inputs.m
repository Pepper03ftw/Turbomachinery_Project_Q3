function in = default_fan_inputs()
%DEFAULT_FAN_INPUTS Default inputs for the AE4206 fan meanline code.
% All user-tunable settings, switches, and assumed loss-model inputs live
% here on purpose. The meanline solver should compute zeta_R and zeta_S
% from the enabled loss sub-models rather than from hard-coded guessed row
% loss coefficients.

in = struct();

% Project point
in.altitude_m = 10000;
in.M_flight = 0.78;
in.mdot = 80;                 % kg/s
in.beta_tt_target = 1.60;     % total-to-total pressure ratio
in.N_rpm = 5000;

% Meanline coefficients
% In the recommended pressure-ratio/Mach constrained mode, phi and psi are
% computed from Mx_design, Mrel_tip_design, beta_tt_target, N_rpm and the
% efficiency work guess. The values below then act only as fallbacks/initial
% references for other solve modes.
in.phi = 0.70;
in.psi = 0.30;
in.R   = 0.65;

% Design solution mode
%   'solve_phi_from_geometry' : project sizing assumptions define annulus
%                               geometry first, then continuity returns the
%                               self-consistent phi actually carried by that
%                               annulus. This is useful when size is largely
%                               fixed by the chosen geometry assumptions.
%   'solve_area_from_phi'     : the chosen phi is treated as the design value,
%                               and the annulus area/radii are resized so that
%                               continuity is satisfied at that phi. This is
%                               useful when frontal area is allowed to change.
%   'sweep_project_constraints' :
%                               continuity-first helper. Mx, mdot and hub-to-tip
%                               define r_m, so phi is computed rather than swept.
%   'sweep_fixed_mean_radius_duty' :
%                               Denton/MEANGEN-style helper. A fixed design
%                               mean radius is used, phi and psi are swept,
%                               alpha1 can be fixed at zero, and R is then
%                               computed from the velocity-triangle relation.
%   'project_fixed_candidate' :
%                               internal single-candidate mode used by the
%                               sweeps.
%   'solve_phi_psi_from_mach_pr' :
%                               legacy mode retained for back-comparison.
%in.design.solve_mode = 'solve_phi_from_geometry';
%in.design.solve_mode = 'solve_area_from_phi';
%in.design.solve_mode = 'solve_phi_psi_from_mach_pr';
%in.design.solve_mode = 'sweep_project_constraints';
in.design.solve_mode = 'sweep_fixed_mean_radius_duty';

% Project geometry closure.
% use_continuity_fixed_radius keeps the older continuity-first closure in
% which Mx, mdot and hub-to-tip ratio define the annulus and r_m, so phi is
% computed from Vx/(Omega*r_m).  The Denton/MEANGEN-style mode below instead
% uses a fixed mean radius and sweeps phi directly.
in.design.use_continuity_fixed_radius = true;

% Fixed-mean-radius duty sweep settings.  This follows the MULTALL/MEANGEN
% philosophy: provide a preliminary mean radius, flow coefficient, loading
% coefficient and reaction/angle set, then let MEANGEN/STAGEN/MULTALL refine
% blade shape and annulus details.  If fixed_mean_radius_m is NaN, the code
% uses the continuity radius as a reference and multiplies it by
% fixed_mean_radius_scale.  Set fixed_mean_radius_m to impose a larger radius
% directly, e.g. 0.45 or 0.48 m.
in.design.fixed_mean_radius_m = NaN;       % NaN -> use continuity reference radius
in.design.fixed_mean_radius_scale = 1.00;  % applied only when fixed_mean_radius_m is NaN
in.design.fixed_mean_radius_note = 'Set fixed_mean_radius_m to impose MEANGEN design radius explicitly.';

% Sizing basis suggested by project
in.Mx_design = 0.60;
in.Mrel_tip_design = 1.40;
in.hub_to_tip = 0.30;          % project value; treated as a fixed constraint in sweep mode

% Work and duty-coefficient limits for the pressure-ratio/Mach constrained mode.
% The target pressure ratio fixes the isentropic work. The Euler work is
% estimated as w_is/eta_work_guess, giving psi = Delta_h0/U_m^2.
% Later, eta_work_guess can be replaced by MULTALL feedback or an outer
% fixed-point loop.
in.design.eta_work_guess = 0.90;
% Pressure-ratio selection.
% beta_tt_preferred is the nominal target. beta_tt_min/max are admissible
% values. In auto mode, the code searches the range and selects the best
% feasible design while applying a soft preference for beta_tt_preferred.
in.design.beta_tt_mode = 'achieved_from_phi_psi_eta';   % sweep mode computes achieved beta_tt from psi and losses
in.design.beta_tt_preferred = 1.60;
in.design.beta_tt_min = 1.40;
in.design.beta_tt_max = 1.70;
in.design.beta_tt_search_points = 61;
in.design.beta_tt_penalty_weight = 1.0;   % soft preference only; not a hard constraint by default

% Inlet total-condition convention. The current code default preserves the
% earlier convention: cruise Mach is used to compute inlet total conditions
% from the ISA static atmosphere. Set to 'ambient_static' only if you
% intentionally want fan-face total conditions equal to ISA static values.
in.design.inlet_total_mode = 'freestream_total';

in.design.psi_min = 0.30;
in.design.psi_max = 0.60;
in.design.phi_min = 0.40;
in.design.phi_max = 1.10;

% Duty-coefficient sweep resolution for the fast preliminary search.
in.design.phi_search_points = 15;   % coarse backbone; augmented near mdot target below
in.design.psi_search_points = 11;
in.design.R_search_points = 21;

% Add local phi samples around the value implied by mdot, Mx, RPM and htr.
% This keeps the sweep fast while avoiding a large mdot miss caused only by
% coarse phi spacing. The exact phi is not forced; it is just inserted into
% the candidate list and then ranked by the objective.
in.design.phi_grid_include_mdot_target = true;
in.design.phi_grid_mdot_refine_half_width = 0.08;
in.design.phi_grid_mdot_refine_points = 9;

% Explore the no-IGV option without making it mandatory.
%   'free' : alpha1 follows from phi, psi and R.
%   'zero' : enforce alpha1 = alpha1_deg, so R follows from R=1-phi*tan(alpha1)-psi/2.
%   'both' : search both families and let constraints/objective choose.
in.design.alpha1_mode = 'zero';   % hard no-IGV default for the fixed-radius duty sweep
in.design.alpha1_deg = 0;

% Fixed project quantities are allowed to deviate only by this relative tolerance.
in.design.fixed_tolerance_rel = 0.02;

% Recovery ratio is reported as a soft stability metric, not a hard rejection.
in.design.RR_soft_min = 0.50;
% Recovery-ratio stability is still soft by default, but it is now weighted
% strongly enough that negative RR candidates should not win unless no
% reasonable alternative exists. Set to 0 to ignore RR in ranking.
in.design.RR_penalty_weight = 2.00;
% RR stability equation convention. Default matches the course derivation:
% alpha1 remains signed absolute inlet swirl, while beta2_RR = -beta2_code
% because Vtheta2 = U - Vx*tan(beta2_RR). 'signed' and
% 'positive_magnitudes_legacy' are diagnostic-only sensitivity options.
in.design.RR_angle_convention = 'lecture_signed_alpha_beta2_positive';
% Optional hard rejection switch. With the corrected positive-angle RR
% convention, candidates below RR_soft_min are rejected by default because
% they are poor stability starting points. Set false for exploratory sweeps.
in.design.enforce_RR_as_hard = true;

% Soft duty preference. This prevents a least-loss/least-violation sweep from
% automatically sitting on R = 1 unless that is truly needed. Set the weight
% to 0 to disable it.
in.design.R_preferred = 0.50;
in.design.R_penalty_weight = 0.03;

% Command-window progress messages.
in.design.status_print = true;
in.design.status_every = 1000;

% Meanline mass-flow constraint. In the recommended Mach/PR mode, Vx is
% fixed by Mx and density, while mdot fixes the annulus area. The hub radius
% is therefore solved from area once the tip radius has been selected to
% satisfy Mrel_tip. This makes hub-to-tip ratio an output of this solve mode.
in.design.enforce_mdot_by_hub_radius = true;
in.design.mdot_tolerance_abs = 0.02*in.mdot;   % kg/s, project +/-2% tolerance
in.design.mdot_tolerance_rel = 0.02;           % project +/-2% tolerance
% In the preliminary sweep the project matching quantities are reported and
% optionally penalized, but are not treated as hard feasibility constraints by
% default. This keeps the sweep useful even when a single-stage meanline model
% cannot hit beta_tt, mdot and Mrel_tip simultaneously.
in.design.include_mdot_constraint = true;
in.design.enforce_project_matching_as_hard = false;
in.design.enforce_beta_tt_range_as_hard = false;
% Absolute-error penalties are disabled by default so that values inside
% the +/- fixed_tolerance_rel project band are not punished. Excess penalties
% below activate only after the 2% tolerance band is exceeded.
in.design.mdot_penalty_weight = 0.00;
in.design.Mrel_tip_penalty_weight = 0.00;
in.design.mdot_excess_penalty_weight = 25.0;
in.design.Mrel_tip_excess_penalty_weight = 6.0;
in.design.beta_tt_excess_penalty_weight = 4.0;

% Mean-radius interpretation used after hub and tip radii are known.
%   'arithmetic' : r_m = 0.5*(r_t + r_h)
%   'area'       : r_m = sqrt(0.5*(r_t^2 + r_h^2))
in.design.mean_radius_rule = 'arithmetic';

% Tip relative Mach matching. In the continuity-first baseline, mdot, rho,
% Vx and hub-to-tip ratio define a reference r_m. For the MEANGEN/MULTALL
% starting point, the code may slightly rescale r_m about that continuity
% value to match the project Mrel_tip target. The mass-flow mismatch created
% by this rescaling is reported, but by default is not used to rank/reject
% candidates because the downstream fixed-mdot simulations can change blade
% shape/throttle to carry the prescribed mass flow.
in.design.match_Mrel_tip_by_mean_radius = false;
in.design.Mrel_tip_radius_match_max_rel = 0.20;       % +/- allowed r_m change around continuity radius
in.design.Mrel_tip_radius_match_tol_rel = 0.02;       % same +/-2% tolerance used for project checks
in.design.enforce_Mrel_tip_as_hard = false;           % fixed-radius mode reports this; set true once radius can satisfy it
in.design.relax_mdot_penalty_when_matching_tip_radius = true;
in.design.relax_mdot_penalty_in_fixed_radius_duty = true;


% Feasibility / constraint settings
% These checks are returned as residuals that can later be passed directly
% to an optimizer such as fmincon. A residual <= 0 means the constraint is
% satisfied, while a residual > 0 means the design violates the constraint.
%
% DF_sizing_target is used only to compute the initial solidity estimate.
% After the solidity is iterated for Howell feasibility, the final DF must
% remain within DF_min <= DF <= DF_limit.
in.design.DF_sizing_target = 0.45;
in.design.DF_min = 0.20;
in.design.DF_limit = 0.55;
in.design.enforce_hard_limits = false;  % if true, solver errors on violated checks

% Reaction is not prescribed by the project. It can either be fixed or
% searched to balance rotor/stator DF and Howell turning constraints.
in.design.R_mode = 'fixed';  % sweep mode controls R directly
in.design.R_min = 0.00;
in.design.R_max = 1.00;
% in.design.R_search_points is defined above for the sweep.

% Solidity sizing. The default policy computes the solidity from the target
% DF (0.45), then only increases solidity if needed to satisfy Howell turning
% while remaining inside DF_min <= DF <= DF_limit. This avoids blindly using
% the lowest feasible solidity.
in.design.solidity_mode = 'from_DF_Howell';  % 'fixed' or 'from_DF_Howell'
in.design.solidity_policy = 'DF_target_then_Howell';  % alternatives: 'closest_to_df_target', 'minimum_solidity_feasible'
in.design.solidity_min = 0.80;
in.design.solidity_max = 3.50;
in.design.solidity_search_points = 250;   % was 1200; vectorized search, adequate for preliminary sweep

% Exit condition choice
in.alpha3_deg = 0;

% Plot switches
in.plot.enable_duty_sweep = false;
% Diagnostic plot of beta_tt versus mdot from the preliminary sweep history.
% This is NOT a true fixed-geometry compressor characteristic unless the
% sampled points come from an off-design/fixed-geometry calculation. It is
% Speedline / beta_tt-mdot plotting
% A true compressor speedline requires fixed-geometry off-design points
% (different throttle/back-pressure/mass-flow solutions at the same speed).
% The preliminary design sweep changes duty variables and sometimes geometry,
% so it must not be interpreted as a choke/surge speedline. Keep the old
% diagnostic off by default to avoid a misleading plot.
in.plot.enable_beta_mdot_map = false;
in.plot.beta_mdot_only_feasible = true;
in.plot.beta_mdot_save_png = true;
in.plot.beta_mdot_png_name = 'beta_tt_vs_mdot_preliminary_sweep.png';

% Geometry / blade-count closure
% The active closure below uses solidity from the DF/Howell loop. For each
% integer blade-count combination in the sweep, pitch follows from
% s = 2*pi*r_m/Z and chord follows from c = sigma*s. A fixed t/c assumption
% then sets blade thickness for the Lieblein and loss-model inputs.
%
% The chord values below are only fallbacks/initial guesses for the fixed
% chord mode or for the first Re estimate before the reverse closure updates
% them.
in.rotor_chord = 0.060;       % m, fallback/initial guess only in blade-count closure
in.stator_chord = 0.050;      % m, fallback/initial guess only in blade-count closure
in.rotor_solidity = 1.35;     % c/s, overwritten when solidity_mode='from_DF_Howell'
in.stator_solidity = 1.25;    % c/s, overwritten when solidity_mode='from_DF_Howell'
in.rotor_pitch = in.rotor_chord / in.rotor_solidity;
in.stator_pitch = in.stator_chord / in.stator_solidity;
in.rotor_thickness = 0.10*in.rotor_chord;
in.stator_thickness = 0.10*in.stator_chord;
in.blade_height_factor = 0.9;

% Reverse closure from solidity to chord using integer blade counts.
in.design.chord_closure_mode = 'blade_count_from_solidity';  % 'fixed_chord' or 'blade_count_from_solidity'

% Blade-count handling:
%   'sweep_in_objective'      : old behavior; Z_R/Z_S are swept and can affect the selected objective.
%   'postprocess_table_only'  : recommended behavior; the duty/solidity solution is selected first,
%                               then all blade-count combinations are tabulated afterwards. This avoids
%                               artificial lower/upper Z preferences from weak meanline loss sensitivity.
in.design.blade_count_selection_mode = 'postprocess_table_only';
in.design.blade_count_sweep = true; % only used when blade_count_selection_mode='sweep_in_objective'
in.design.rotor_blade_count_range = 20:40;
in.design.stator_vane_count_range = 20:60;
in.design.rotor_blade_count_default = 30;  % reference only in postprocess_table_only mode
in.design.stator_vane_count_default = 40;  % reference only in postprocess_table_only mode
in.design.t_over_c = 0.10;
in.design.rotor_t_over_c = 0.10;
in.design.stator_t_over_c = 0.10;
in.design.chord_solidity_iterations = 6;

% Leakage defaults
in.rotor_tip_clearance = 0.0005;   % m
in.stator_hub_clearance = 0.0;     % m

% Spanwise / radial evaluation defaults for preliminary hub-mean-tip loss models
% These use the current meanline triangles to create provisional radial
% sections. Once MEANGEN/STAGEN/MULTALL are connected, the generated
% section-resolved angles, Mach numbers, pitches, and blade metal angles can
% replace these estimates without changing the loss-model interfaces.
in.loss_model.radial.vtheta_mode = 'constant_vtheta';  % 'constant_vtheta' or 'free_vortex'
in.loss_model.radial.pitch_mode = 'scale_with_radius'; % blade count fixed: pitch proportional to radius

% Air model
in.air.gamma = 1.4;
in.air.R = 287.05;
in.air.cp = in.air.gamma*in.air.R/(in.air.gamma-1);

% Material for crude centrifugal force estimate
in.material.rho_blade = 4500;

% Howell digitized data
base_dir = fileparts(mfilename('fullpath'));
in.howell_paths.f_beta2 = fullfile(base_dir,'Howell_f_beta2.csv');
in.howell_paths.Phi_Re  = fullfile(base_dir,'Howell_Phi_Re.csv');
in.howell_paths.Psi_sc  = fullfile(base_dir,'Howell_Psi_solidityINVERSE.csv');

% Lieblein digitized incidence/deviation data. The files are resolved
% relative to this default_fan_inputs.m file, so the code is not sensitive
% to the MATLAB current working directory.
in.lieblein.enabled = true;
in.lieblein.profile_family = 'NACA65';      % 'NACA65' or 'DCA'
in.lieblein.default_tmax_over_c = 0.100;    % used if row-specific t/c is not supplied
in.lieblein.rotor_tmax_over_c = NaN;        % NaN -> rotor_thickness/rotor_chord
in.lieblein.stator_tmax_over_c = NaN;       % NaN -> default_tmax_over_c
in.lieblein.Ki_sh_DCA = 0.7;
in.lieblein.Ki_sh_NACA65 = 1.1;
in.lieblein.Kdelta_sh_DCA = 0.75;
in.lieblein.Kdelta_sh_NACA65 = 1.1;

lieblein_dir = fullfile(base_dir,'lieblein_digitized');
in.lieblein_paths.folder = lieblein_dir;
in.lieblein_paths.i0_10 = fullfile(lieblein_dir,'lieblein_i0_10.csv');
in.lieblein_paths.n = fullfile(lieblein_dir,'lieblein_n.csv');
in.lieblein_paths.Ki_t = fullfile(lieblein_dir,'lieblein_Ki_t.csv');
in.lieblein_paths.delta0_10 = fullfile(lieblein_dir,'lieblein_delta0_10.csv');
in.lieblein_paths.m_NACA65 = fullfile(lieblein_dir,'lieblein_m_NACA65.csv');
in.lieblein_paths.Kdelta_t = fullfile(lieblein_dir,'lieblein_Kdelta_t.csv');
in.lieblein_paths.b = fullfile(lieblein_dir,'lieblein_b.csv');

% loss_model.type is a documentation label only. It is printed in outputs so
% the user can record which literature-based configuration was intended, but
% it does NOT switch the active loss equations on/off. That is controlled only
% by the enabled flags below.
in.loss_model.type = 'modular_literature';

% If true, and if the Lieblein digitized charts evaluate successfully, the
% computed optimum incidence/deviation are automatically copied into the
% tip, endwall, and shock stagger calculations. This makes the loss models
% use signed Lieblein metal-angle stagger without manually copying values.
in.loss_model.auto_use_lieblein_stagger = true;
in.loss_model.lieblein_stagger_auto_applied = false;

% Published / assumed sub-model inputs for the modular loss build-up.
% The comments below distinguish between:
%   - directly literature-backed defaults, and
%   - first-pass assumed geometric / empirical inputs to be updated later.

% 1) Denton idealized BL/profile loss (viscous dissipation)
in.loss_model.profile.enabled = true;
in.loss_model.profile.Cd_bl = 0.002;      % literature-backed default
in.loss_model.profile.dV_over_Vbar = 1/sqrt(3); % lecture optimum for first pass

% 2) Denton TE / wake-mixing loss (mixing out loss)
in.loss_model.trailing_edge.enabled = true;
in.loss_model.trailing_edge.use_direct_t_over_s = true;
in.loss_model.trailing_edge.t_over_s = 0.10;    % assumed initial value
in.loss_model.trailing_edge.Cpb = -0.15;        % assumed initial value
in.loss_model.trailing_edge.use_bl_terms = false;
in.loss_model.trailing_edge.theta_over_s = 0.0; % only used if use_bl_terms = true
in.loss_model.trailing_edge.dstar_over_s = 0.0; % only used if use_bl_terms = true

% 3) Shock loss model
% Supported active model:
%   'koch_smith_le_and_passage' : Koch-Smith compressor shock model. It
%                                 separately estimates leading-edge
%                                 bluntness shock entropy and passage-shock
%                                 entropy. The passage-shock component uses
%                                 a representative passage inlet Mach number
%                                 based primarily on maximum suction-surface
%                                 Mach number.
in.loss_model.shock.enabled = true;
in.loss_model.shock.enabled_rotor = true;
in.loss_model.shock.enabled_stator = false;    % switch on only if stator shocks are expected
in.loss_model.shock.model = 'koch_smith_le_and_passage';
in.loss_model.shock.M_crit = 1.0;              % LE bluntness active for M1 > M_crit

% Section selection. Current default evaluates rotor shocks at the tip. The
% following options are supported: 'hub', 'mean', 'tip', 'hub_tip',
% 'hub_mean_tip'. Once MEANGEN/STAGEN/MULTALL section data are available,
% use 'hub_mean_tip' and feed section-resolved Mss/geometry values.
in.loss_model.shock.evaluate_at = 'tip';
in.loss_model.shock.rotor_evaluate_at = 'tip';
in.loss_model.shock.stator_evaluate_at = 'mean';

% Koch-Smith switches
in.loss_model.shock.use_le_bluntness = true;
in.loss_model.shock.use_passage_shock = true;
in.loss_model.shock.Mss_weight = 6.0;          % M_rep = (6*Mss_max + M1)/7
in.loss_model.shock.temperature_reference = 'inlet_static'; % 'inlet_static', 'mean_static', 'surface_static'

% Leading-edge thickness. Koch-Smith Eq. (1) uses t_LE/(b cos beta1), where
% b is tangential spacing. The code accepts direct t_LE/pitch values, or
% estimates t_LE/pitch = (t_LE/chord)*(chord/pitch). Default t_LE/chord=0.008
% is a preliminary thin transonic fan/compressor assumption.
in.loss_model.shock.tLE_over_chord = 0.008;
in.loss_model.shock.rotor_tLE_over_chord = 0.008;
in.loss_model.shock.stator_tLE_over_chord = 0.008;
in.loss_model.shock.tLE_over_pitch = NaN;
in.loss_model.shock.rotor_tLE_over_pitch = NaN;
in.loss_model.shock.stator_tLE_over_pitch = NaN;
in.loss_model.shock.rotor_tLE_over_pitch_hub = NaN;
in.loss_model.shock.rotor_tLE_over_pitch_mean = NaN;
in.loss_model.shock.rotor_tLE_over_pitch_tip = NaN;
in.loss_model.shock.stator_tLE_over_pitch_hub = NaN;
in.loss_model.shock.stator_tLE_over_pitch_mean = NaN;
in.loss_model.shock.stator_tLE_over_pitch_tip = NaN;

% Maximum suction-surface Mach number used by the Koch-Smith passage-shock
% model. Current default estimates Mss_max from the same low-order surface
% velocity split used by the tip/endwall models. Later, set Mss_mode='direct'
% and supply rotor_Mss_max_hub/mean/tip from MULTALL, STAGEN postprocessing,
% or another blade-to-blade calculation.
% Supported: 'estimated', 'direct', 'factor'.
in.loss_model.shock.Mss_mode = 'estimated';
in.loss_model.shock.Mss_factor = 1.0;          % only used when Mss_mode='factor'
in.loss_model.shock.rotor_Mss_max = NaN;
in.loss_model.shock.rotor_Mss_max_hub = NaN;
in.loss_model.shock.rotor_Mss_max_mean = NaN;
in.loss_model.shock.rotor_Mss_max_tip = NaN;
in.loss_model.shock.stator_Mss_max = NaN;
in.loss_model.shock.stator_Mss_max_hub = NaN;
in.loss_model.shock.stator_Mss_max_mean = NaN;
in.loss_model.shock.stator_Mss_max_tip = NaN;

% Stagger approximation used only for Mss_mode='estimated'. Replace with
% section metal angles once blade profiles are generated.
% Stagger is now treated as a signed blade-geometry angle.
% Available modes:
%   'signed_flow_angles'       : signed fallback, 0.5*(chi_in+chi_out)
%   'lieblein_metal_angles'    : uses incidence/deviation as positive-angle chart values
%   'signed_metal_angles'      : uses incidence/deviation as already signed values
%   'user' / 'user_sections'   : prescribed signed stagger(s)
%   'magnitude_mean_flow_angles' : legacy absolute-angle fallback
in.loss_model.shock.stagger_mode = 'signed_flow_angles';
in.loss_model.shock.auto_use_lieblein_stagger = true;
in.loss_model.shock.rotor_stagger_deg = NaN;
in.loss_model.shock.stator_stagger_deg = NaN;
in.loss_model.shock.rotor_stagger_hub_deg = NaN;
in.loss_model.shock.rotor_stagger_mean_deg = NaN;
in.loss_model.shock.rotor_stagger_tip_deg = NaN;
in.loss_model.shock.stator_stagger_hub_deg = NaN;
in.loss_model.shock.stator_stagger_mean_deg = NaN;
in.loss_model.shock.stator_stagger_tip_deg = NaN;

% Optional incidence/deviation inputs for shock.stagger_mode =
% 'lieblein_metal_angles' or 'signed_metal_angles'. Values are in degrees.
% For 'lieblein_metal_angles', these are interpreted as positive-angle
% Lieblein chart values and internally mapped to the signed row convention.
in.loss_model.shock.rotor_incidence_deg = NaN;
in.loss_model.shock.rotor_deviation_deg = NaN;
in.loss_model.shock.stator_incidence_deg = NaN;
in.loss_model.shock.stator_deviation_deg = NaN;
in.loss_model.shock.rotor_incidence_hub_deg = NaN;
in.loss_model.shock.rotor_incidence_mean_deg = NaN;
in.loss_model.shock.rotor_incidence_tip_deg = NaN;
in.loss_model.shock.rotor_deviation_hub_deg = NaN;
in.loss_model.shock.rotor_deviation_mean_deg = NaN;
in.loss_model.shock.rotor_deviation_tip_deg = NaN;
in.loss_model.shock.stator_incidence_hub_deg = NaN;
in.loss_model.shock.stator_incidence_mean_deg = NaN;
in.loss_model.shock.stator_incidence_tip_deg = NaN;
in.loss_model.shock.stator_deviation_hub_deg = NaN;
in.loss_model.shock.stator_deviation_mean_deg = NaN;
in.loss_model.shock.stator_deviation_tip_deg = NaN;

% Legacy weak-shock setting retained only if shock.model is manually set to
% 'weak_shock_entropy'.
in.loss_model.shock.apply_to = 'outlet';

% 4) Tip-leakage loss estimate
% Supported models:
%   'legacy_area_proxy'             : old one-line leakage-area proxy.
%   'denton_hall_unshrouded_optionA': Denton/Hall unshrouded leakage model
%                                      using Appendix-B V_SS/Vx and V_PS/Vx
%                                      approximations, with optional direct
%                                      surface-velocity input.
in.loss_model.tip.enabled_rotor = true;
in.loss_model.tip.enabled_stator = false;
in.loss_model.tip.model = 'denton_hall_unshrouded_optionA';

% Legacy model coefficient. This is kept only for back-comparison with the
% previous code path and is not used by the Denton/Hall model.
in.loss_model.tip.Cd_tip = 0.002;

% Denton/Hall unshrouded model settings.
in.loss_model.tip.Cd_leak = 0.8;                 % leakage discharge coefficient
in.loss_model.tip.Cs_over_c = 1.0;               % blade surface length / chord; update when camber geometry is available
in.loss_model.tip.surface_velocity_mode = 'approx';  % 'approx' or 'direct'
% Stagger is now signed. Set to 'lieblein_metal_angles' after copying the
% computed Lieblein incidence/deviation into the fields below.
in.loss_model.tip.stagger_mode = 'signed_flow_angles'; % see shock.stagger_mode options
in.loss_model.tip.auto_use_lieblein_stagger = true;
in.loss_model.tip.rotor_stagger_deg = NaN;       % only used if stagger_mode = 'user'
in.loss_model.tip.stator_stagger_deg = NaN;      % only used if stagger_mode = 'user'
in.loss_model.tip.rotor_incidence_deg = NaN;
in.loss_model.tip.rotor_deviation_deg = NaN;
in.loss_model.tip.stator_incidence_deg = NaN;
in.loss_model.tip.stator_deviation_deg = NaN;

% Direct surface-velocity placeholders. These allow future MULTALL/StageN
% output to be connected without changing the loss-model logic. Values are
% normalized by Vx. Scalars give the option-A algebraic form; vectors trigger
% numerical integration of the original Denton/Hall integrand over x/Cs.
in.loss_model.tip.rotor_uSS_over_Vx = NaN;
in.loss_model.tip.rotor_uPS_over_Vx = NaN;
in.loss_model.tip.rotor_x_over_Cs = NaN;
in.loss_model.tip.stator_uSS_over_Vx = NaN;
in.loss_model.tip.stator_uPS_over_Vx = NaN;
in.loss_model.tip.stator_x_over_Cs = NaN;

% 5) Endwall boundary-layer loss baseline
% Disabled by default in v8 because the preliminary optimizer can otherwise
% use the simplified chord-dependent endwall term to bias blade-count/chord
% selection. Keep the model available for sensitivity checks, but do not use
% it as a geometry-selection driver.
% Supported models:
%   'hall_denton_cd_baseline' : Hall/Denton constant-C_D endwall BL
%                               dissipation only. This is NOT a full
%                               secondary-flow/corner-separation model.
%   'user_constant'           : row-wise constant zeta placeholders below.
in.loss_model.endwall.enabled_rotor = false;
in.loss_model.endwall.enabled_stator = false;
in.loss_model.endwall.model = 'hall_denton_cd_baseline';
in.loss_model.endwall.secondary_model = 'none';       % reserved for future MULTALL/STAGEN extraction

% Hall/Denton constant dissipation-coefficient assumption for turbulent
% endwall boundary layers.
in.loss_model.endwall.Cd = 0.002;
in.loss_model.endwall.Cs_over_c = 1.0;                % endwall wetted length / chord; update with blade geometry
in.loss_model.endwall.surface_velocity_mode = 'approx';% 'approx', 'direct', or 'tip_shared'
% Stagger is now signed. Set to 'lieblein_metal_angles' after copying the
% computed Lieblein incidence/deviation into the fields below.
in.loss_model.endwall.stagger_mode = 'signed_flow_angles'; % see shock.stagger_mode options
in.loss_model.endwall.auto_use_lieblein_stagger = true;
in.loss_model.endwall.rotor_stagger_deg = NaN;        % only used if stagger_mode = 'user'
in.loss_model.endwall.stator_stagger_deg = NaN;       % only used if stagger_mode = 'user'
in.loss_model.endwall.rotor_incidence_deg = NaN;
in.loss_model.endwall.rotor_deviation_deg = NaN;
in.loss_model.endwall.stator_incidence_deg = NaN;
in.loss_model.endwall.stator_deviation_deg = NaN;

% Endwall radial evaluation. Default now evaluates hub and casing/tip as
% separate endwall surfaces. Supported options: 'mean_combined' (old meanline
% treatment), 'hub', 'tip', 'hub_tip', and 'hub_mean_tip'.
in.loss_model.endwall.evaluate_at = 'hub_tip';
in.loss_model.endwall.rotor_evaluate_at = 'hub_tip';
in.loss_model.endwall.stator_evaluate_at = 'hub_tip';
in.loss_model.endwall.rotor_n_endwalls = 2;           % used only by mean_combined legacy mode
in.loss_model.endwall.stator_n_endwalls = 2;          % used only by mean_combined legacy mode

% Constant-zeta fallback, used only with model = 'user_constant'.
in.loss_model.endwall.rotor_zeta = 0.0;
in.loss_model.endwall.stator_zeta = 0.0;

% Direct endwall edge-velocity placeholders for future MULTALL/StageN data.
% Values are normalized by Vx. Scalars give an algebraic passage average;
% vectors trigger x/Cs integration with linear pressure-to-suction-side
% interpolation across the passage at each x-location.
in.loss_model.endwall.rotor_uSS_over_Vx = NaN;
in.loss_model.endwall.rotor_uPS_over_Vx = NaN;
in.loss_model.endwall.rotor_x_over_Cs = NaN;
in.loss_model.endwall.stator_uSS_over_Vx = NaN;
in.loss_model.endwall.stator_uPS_over_Vx = NaN;
in.loss_model.endwall.stator_x_over_Cs = NaN;

end
