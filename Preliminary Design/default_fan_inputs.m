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
%in.design.solve_mode = 'solve_phi_from_geometry';
in.design.solve_mode = 'solve_area_from_phi';

% Sizing basis suggested by project
in.Mx_design = 0.60;
in.Mrel_tip_design = 1.40;
in.hub_to_tip = 0.30;

% Mean-radius interpretation used after hub and tip radii are known.
%   'arithmetic' : r_m = 0.5*(r_t + r_h)
%   'area'       : r_m = sqrt(0.5*(r_t^2 + r_h^2))
in.design.mean_radius_rule = 'arithmetic';


% Feasibility / constraint settings
% These checks are returned as residuals that can later be passed directly
% to an optimizer such as fmincon. A residual <= 0 means the constraint is
% satisfied, while a residual > 0 means the design violates the constraint.
in.design.DF_limit = 0.45;
in.design.enforce_hard_limits = false;  % if true, solver errors on violated checks

% Exit condition choice
in.alpha3_deg = 0;

% Geometry guesses
in.rotor_chord = 0.060;       % m
in.stator_chord = 0.050;      % m
in.rotor_solidity = 1.35;     % c/s
in.stator_solidity = 1.25;    % c/s
in.rotor_pitch = in.rotor_chord / in.rotor_solidity;
in.stator_pitch = in.stator_chord / in.stator_solidity;
in.rotor_thickness = 0.004;
in.blade_height_factor = 0.9;

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
in.howell_paths.f_beta2 = 'Howell_f_beta2.csv';
in.howell_paths.Phi_Re  = 'Howell_Phi_Re.csv';
in.howell_paths.Psi_sc  = 'Howell_Psi_solidityINVERSE.csv';

% loss_model.type is a documentation label only. It is printed in outputs so
% the user can record which literature-based configuration was intended, but
% it does NOT switch the active loss equations on/off. That is controlled only
% by the enabled flags below.
in.loss_model.type = 'modular_literature';

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
in.loss_model.shock.stagger_mode = 'mean_flow_angles'; % 'mean_flow_angles', 'user', or 'user_sections'
in.loss_model.shock.rotor_stagger_deg = NaN;
in.loss_model.shock.stator_stagger_deg = NaN;
in.loss_model.shock.rotor_stagger_hub_deg = NaN;
in.loss_model.shock.rotor_stagger_mean_deg = NaN;
in.loss_model.shock.rotor_stagger_tip_deg = NaN;
in.loss_model.shock.stator_stagger_hub_deg = NaN;
in.loss_model.shock.stator_stagger_mean_deg = NaN;
in.loss_model.shock.stator_stagger_tip_deg = NaN;

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
in.loss_model.tip.stagger_mode = 'mean_flow_angles'; % 'mean_flow_angles' or 'user'
in.loss_model.tip.rotor_stagger_deg = NaN;       % only used if stagger_mode = 'user'
in.loss_model.tip.stator_stagger_deg = NaN;      % only used if stagger_mode = 'user'

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
% Supported models:
%   'hall_denton_cd_baseline' : Hall/Denton constant-C_D endwall BL
%                               dissipation only. This is NOT a full
%                               secondary-flow/corner-separation model.
%   'user_constant'           : row-wise constant zeta placeholders below.
in.loss_model.endwall.enabled_rotor = true;
in.loss_model.endwall.enabled_stator = true;
in.loss_model.endwall.model = 'hall_denton_cd_baseline';
in.loss_model.endwall.secondary_model = 'none';       % reserved for future MULTALL/STAGEN extraction

% Hall/Denton constant dissipation-coefficient assumption for turbulent
% endwall boundary layers.
in.loss_model.endwall.Cd = 0.002;
in.loss_model.endwall.Cs_over_c = 1.0;                % endwall wetted length / chord; update with blade geometry
in.loss_model.endwall.surface_velocity_mode = 'approx';% 'approx', 'direct', or 'tip_shared'
in.loss_model.endwall.stagger_mode = 'mean_flow_angles'; % 'mean_flow_angles' or 'user'
in.loss_model.endwall.rotor_stagger_deg = NaN;        % only used if stagger_mode = 'user'
in.loss_model.endwall.stator_stagger_deg = NaN;       % only used if stagger_mode = 'user'

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
