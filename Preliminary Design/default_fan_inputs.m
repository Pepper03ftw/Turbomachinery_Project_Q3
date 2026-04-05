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
in.design.DF_limit = 0.60;
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
in.loss_model.profile.dV_over_Vbar = 1/3; % lecture optimum for first pass

% 2) Denton TE / wake-mixing loss (mixing out loss)
in.loss_model.trailing_edge.enabled = true;
in.loss_model.trailing_edge.use_direct_t_over_s = true;
in.loss_model.trailing_edge.t_over_s = 0.10;    % assumed initial value
in.loss_model.trailing_edge.Cpb = -0.15;        % assumed initial value
in.loss_model.trailing_edge.use_bl_terms = false;
in.loss_model.trailing_edge.theta_over_s = 0.0; % only used if use_bl_terms = true
in.loss_model.trailing_edge.dstar_over_s = 0.0; % only used if use_bl_terms = true

% 3) Shock loss from normal-shock entropy relation
in.loss_model.shock.enabled = true;
in.loss_model.shock.M_crit = 1.0;               % apply only above this Mach
in.loss_model.shock.apply_to = 'outlet';        % 'outlet' or 'max'

% 4) Approximate leakage-loss estimate (tip-leakage)
in.loss_model.tip.enabled_rotor = true;
in.loss_model.tip.enabled_stator = false;
in.loss_model.tip.Cd_tip = 0.002;               % first-pass choice, same order as BL Cd

% 5) Endwall / secondary loss placeholder (secondary flows loss)
in.loss_model.endwall.enabled_rotor = false;
in.loss_model.endwall.enabled_stator = false;
in.loss_model.endwall.rotor_zeta = 0.0;         % only used if enabled
in.loss_model.endwall.stator_zeta = 0.0;        % only used if enabled

end
