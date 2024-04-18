# ==============================================================================
# Model description
# Molten Salt Reactor Experiment (MSRE) Model - Steady-State Model
# Primary Loop Thermal Hydraulics Model
# Integrates:
# - Transport of species
# ==============================================================================
# Author(s): Dr. Mauricio Tano, Dr. Samuel Walker
# ==============================================================================
# ==============================================================================
# MODEL PARAMETERS
# ==============================================================================

# Mass flow rate tuning
D_mu = 1e-6 # species molecular diffusivity
D_mu_g = 1e-1 # species molecular diffusivity in gaseous phase

# Numerical scheme parameters
advected_interp_method = 'upwind'
velocity_interp_method = 'rc'

# Dynamic scaling paramters
Sc_t = 0.9

# Average bubble diameter
dp = 0.01

# Blocks
fluid_blocks = 'core lower_plenum upper_plenum down_comer riser pump elbow'

# Coupled transported species parameters.
effective_fission_capture_rate = 1.0

yield_Br91 = 0.00063954
lambda_Br91 = ${fparse log(2)/(0.543) + effective_fission_capture_rate}
Br91_alpha_fg = 1e-5

yield_Kr91 = 0.01178
lambda_Kr91 = ${fparse log(2)/(8.57) + effective_fission_capture_rate}
Kr91_alpha_fg = 1e-5

yield_Rb91 = 0.0293761
lambda_Rb91 = ${fparse log(2)/(58.2) + effective_fission_capture_rate}
Rb91_alpha_fg = 1e-5

yield_Sr91 = 0.00966988
lambda_Sr91 = ${fparse log(2)/(3600*9.65) + effective_fission_capture_rate}
Sr91_alpha_fg = 1e-5

yield_Y91 = 0.00018865
lambda_Y91 = ${fparse log(2)/(24*3600*58.51) + effective_fission_capture_rate}
Y91_alpha_fg = 1e-5

yield_Zr91 = 0.00000168
lambda_Zr91 = ${fparse effective_fission_capture_rate}
Zr91_alpha_fg = 0.0

[GlobalParams]
  rhie_chow_user_object = 'pins_rhie_chow_interpolator'

  two_term_boundary_expansion = true
  advected_interp_method = ${advected_interp_method}
  velocity_interp_method = ${velocity_interp_method}
  u = superficial_vel_x
  v = superficial_vel_y
  porosity = porosity_var
  pressure = pressure

  mixing_length = 'mixing_length'

  block = ${fluid_blocks}
[]

################################################################################
# GEOMETRY
################################################################################

[Mesh]
  coord_type = 'RZ'
  [fmg]
    type = FileMeshGenerator
    file = 'mesh_in.e'
  []
[]

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]

################################################################################
# EQUATIONS: VARIABLES, KERNELS & BCS
################################################################################

[UserObjects]
  [pins_rhie_chow_interpolator]
    type = PINSFVRhieChowInterpolator
    a_u = a_u_var
    a_v = a_v_var
    a_w = a_w_var
  []
[]

[Variables]
  [Br91_f]
    type = MooseVariableFVReal
    block = ${fluid_blocks}
  []
  [Br91_g]
    type = MooseVariableFVReal
    block = ${fluid_blocks}
  []
  [Kr91_f]
    type = MooseVariableFVReal
    block = ${fluid_blocks}
  []
  [Kr91_g]
    type = MooseVariableFVReal
    block = ${fluid_blocks}
  []
  [Rb91_f]
    type = MooseVariableFVReal
    block = ${fluid_blocks}
  []
  [Rb91_g]
    type = MooseVariableFVReal
    block = ${fluid_blocks}
  []
  [Sr91_f]
    type = MooseVariableFVReal
    block = ${fluid_blocks}
  []
  [Sr91_g]
    type = MooseVariableFVReal
    block = ${fluid_blocks}
  []
  [Y91_f]
    type = MooseVariableFVReal
    block = ${fluid_blocks}
  []
  [Y91_g]
    type = MooseVariableFVReal
    block = ${fluid_blocks}
  []
  [Zr91_f]
    type = MooseVariableFVReal
    block = ${fluid_blocks}
  []
  [Zr91_g]
    type = MooseVariableFVReal
    block = ${fluid_blocks}
  []
[]

[FVKernels]
  [Br91_time]
    type = FVFunctorTimeKernel
    variable = 'Br91_f'
    block = ${fluid_blocks}
  []
  [Br91_advection]
    type = PINSFVMassAdvection
    variable = 'Br91_f'
    rho = 'Br91_porous'
    block = ${fluid_blocks}
  []
  [Br91_diffusion]
    type = FVDiffusion
    variable = 'Br91_f'
    coeff = ${D_mu}
    block = ${fluid_blocks}
  []
  [Br91_turb_diffusion]
    type = INSFVMixingLengthScalarDiffusion
    variable = 'Br91_f'
    schmidt_number = ${Sc_t}
    block = ${fluid_blocks}
  []
  [Br91_src]
    type = FVCoupledForce
    variable = 'Br91_f'
    v = fission_source
    coef = ${yield_Br91}
    block = ${fluid_blocks}
  []
  [Br91_decay]
    type = FVReaction
    variable = 'Br91_f'
    rate = ${lambda_Br91}
    block = ${fluid_blocks}
  []
  [Br91_fg]
    type = NSFVMixturePhaseInterface
    variable = 'Br91_f'
    phase_coupled = 'Br91_g'
    alpha = 'hlg_Br91'
  []
  [Br91_gf]
    type = NSFVMixturePhaseInterface
    variable = 'Br91_g'
    phase_coupled = 'Br91_f'
    alpha = 'hlg_Br91'
  []
  [Br91_g_advection]
    type = PINSFVMassAdvection
    variable = 'Br91_g'
    rho = 'Br91_g'
    block = ${fluid_blocks}
  []
  [Br91_g_diffusion]
    type = FVDiffusion
    variable = 'Br91_g'
    coeff = ${D_mu_g}
    block = ${fluid_blocks}
  []
  [Br91_g_decay]
    type = FVReaction
    variable = 'Br91_g'
    rate = ${lambda_Br91}
    block = ${fluid_blocks}
  []

  [Kr91_time]
    type = FVFunctorTimeKernel
    variable = 'Kr91_f'
    block = ${fluid_blocks}
  []
  [Kr91_advection]
    type = PINSFVMassAdvection
    variable = 'Kr91_f'
    rho = 'Kr91_porous'
    block = ${fluid_blocks}
  []
  [Kr91_diffusion]
    type = FVDiffusion
    variable = 'Kr91_f'
    coeff = ${D_mu}
    block = ${fluid_blocks}
  []
  [Kr91_turb_diffusion]
    type = INSFVMixingLengthScalarDiffusion
    variable = 'Kr91_f'
    schmidt_number = ${Sc_t}
    block = ${fluid_blocks}
  []
  [Kr91_src]
    type = FVCoupledForce
    variable = 'Kr91_f'
    v = fission_source
    coef = ${yield_Kr91}
    block = ${fluid_blocks}
  []
  [Kr91_decay]
    type = FVReaction
    variable = 'Kr91_f'
    rate = ${lambda_Kr91}
    block = ${fluid_blocks}
  []
  [Kr91_grow]
    type = FVCoupledForce
    variable = 'Kr91_f'
    v = 'Br91_f'
    coef = ${lambda_Br91}
    block = ${fluid_blocks}
  []
  [Kr91_fg]
    type = NSFVMixturePhaseInterface
    variable = 'Kr91_f'
    phase_coupled = 'Kr91_g'
    alpha = 'hlg_Kr91'
  []
  [Kr91_gf]
    type = NSFVMixturePhaseInterface
    variable = 'Kr91_g'
    phase_coupled = 'Kr91_f'
    alpha = 'hlg_Kr91'
  []
  [Kr91_g_advection]
    type = PINSFVMassAdvection
    variable = 'Kr91_g'
    rho = 'Kr91_g'
    block = ${fluid_blocks}
  []
  [Kr91_g_diffusion]
    type = FVDiffusion
    variable = 'Kr91_g'
    coeff = ${D_mu_g}
    block = ${fluid_blocks}
  []
  [Kr91_g_decay]
    type = FVReaction
    variable = 'Kr91_g'
    rate = ${lambda_Kr91}
    block = ${fluid_blocks}
  []
  [Kr91_grow_g]
    type = FVCoupledForce
    variable = 'Kr91_g'
    v = 'Br91_g'
    coef = ${lambda_Br91}
    block = ${fluid_blocks}
  []

  [Rb91_time]
    type = FVFunctorTimeKernel
    variable = 'Rb91_f'
    block = ${fluid_blocks}
  []
  [Rb91_advection]
    type = PINSFVMassAdvection
    variable = 'Rb91_f'
    rho = 'Rb91_porous'
    block = ${fluid_blocks}
  []
  [Rb91_diffusion]
    type = FVDiffusion
    variable = 'Rb91_f'
    coeff = ${D_mu}
    block = ${fluid_blocks}
  []
  [Rb91_turb_diffusion]
    type = INSFVMixingLengthScalarDiffusion
    variable = 'Rb91_f'
    schmidt_number = ${Sc_t}
    block = ${fluid_blocks}
  []
  [Rb91_src]
    type = FVCoupledForce
    variable = 'Rb91_f'
    v = fission_source
    coef = ${yield_Rb91}
    block = ${fluid_blocks}
  []
  [Rb91_decay]
    type = FVReaction
    variable = 'Rb91_f'
    rate = ${lambda_Rb91}
    block = ${fluid_blocks}
  []
  [Rb91_grow]
    type = FVCoupledForce
    variable = 'Rb91_f'
    v = 'Kr91_f'
    coef = ${lambda_Br91}
    block = ${fluid_blocks}
  []
  [Rb91_fg]
    type = NSFVMixturePhaseInterface
    variable = 'Rb91_f'
    phase_coupled = 'Rb91_g'
    alpha = 'hlg_Rb91'
  []
  [Rb91_gf]
    type = NSFVMixturePhaseInterface
    variable = 'Rb91_g'
    phase_coupled = 'Rb91_f'
    alpha = 'hlg_Rb91'
  []
  [Rb91_g_advection]
    type = PINSFVMassAdvection
    variable = 'Rb91_g'
    rho = 'Rb91_g'
    block = ${fluid_blocks}
  []
  [Rb91_g_diffusion]
    type = FVDiffusion
    variable = 'Rb91_g'
    coeff = ${D_mu_g}
    block = ${fluid_blocks}
  []
  [Rb91_g_decay]
    type = FVReaction
    variable = 'Rb91_g'
    rate = ${lambda_Rb91}
    block = ${fluid_blocks}
  []
  [Rb91_grow_g]
    type = FVCoupledForce
    variable = 'Rb91_g'
    v = 'Kr91_g'
    coef = ${lambda_Br91}
    block = ${fluid_blocks}
  []

  [Sr91_time]
    type = FVFunctorTimeKernel
    variable = 'Sr91_f'
    block = ${fluid_blocks}
  []
  [Sr91_advection]
    type = PINSFVMassAdvection
    variable = 'Sr91_f'
    rho = 'Sr91_porous'
    block = ${fluid_blocks}
  []
  [Sr91_diffusion]
    type = FVDiffusion
    variable = 'Sr91_f'
    coeff = ${D_mu}
    block = ${fluid_blocks}
  []
  [Sr91_turb_diffusion]
    type = INSFVMixingLengthScalarDiffusion
    variable = 'Sr91_f'
    schmidt_number = ${Sc_t}
    block = ${fluid_blocks}
  []
  [Sr91_src]
    type = FVCoupledForce
    variable = 'Sr91_f'
    v = fission_source
    coef = ${yield_Sr91}
    block = ${fluid_blocks}
  []
  [Sr91_decay]
    type = FVReaction
    variable = 'Sr91_f'
    rate = ${lambda_Sr91}
    block = ${fluid_blocks}
  []
  [Sr91_grow]
    type = FVCoupledForce
    variable = 'Sr91_f'
    v = 'Rb91_f'
    coef = ${lambda_Br91}
    block = ${fluid_blocks}
  []
  [Sr91_fg]
    type = NSFVMixturePhaseInterface
    variable = 'Sr91_f'
    phase_coupled = 'Sr91_g'
    alpha = 'hlg_Sr91'
  []
  [Sr91_gf]
    type = NSFVMixturePhaseInterface
    variable = 'Sr91_g'
    phase_coupled = 'Sr91_f'
    alpha = 'hlg_Sr91'
  []
  [Sr91_g_advection]
    type = PINSFVMassAdvection
    variable = 'Sr91_g'
    rho = 'Sr91_g'
    block = ${fluid_blocks}
  []
  [Sr91_g_diffusion]
    type = FVDiffusion
    variable = 'Sr91_g'
    coeff = ${D_mu_g}
    block = ${fluid_blocks}
  []
  [Sr91_g_decay]
    type = FVReaction
    variable = 'Sr91_g'
    rate = ${lambda_Sr91}
    block = ${fluid_blocks}
  []
  [Sr91_grow_g]
    type = FVCoupledForce
    variable = 'Sr91_g'
    v = 'Rb91_g'
    coef = ${lambda_Br91}
    block = ${fluid_blocks}
  []

  [Y91_time]
    type = FVFunctorTimeKernel
    variable = 'Y91_f'
    block = ${fluid_blocks}
  []
  [Y91_advection]
    type = PINSFVMassAdvection
    variable = 'Y91_f'
    rho = 'Y91_porous'
    block = ${fluid_blocks}
  []
  [Y91_diffusion]
    type = FVDiffusion
    variable = 'Y91_f'
    coeff = ${D_mu}
    block = ${fluid_blocks}
  []
  [Y91_turb_diffusion]
    type = INSFVMixingLengthScalarDiffusion
    variable = 'Y91_f'
    schmidt_number = ${Sc_t}
    block = ${fluid_blocks}
  []
  [Y91_src]
    type = FVCoupledForce
    variable = 'Y91_f'
    v = fission_source
    coef = ${yield_Y91}
    block = ${fluid_blocks}
  []
  [Y91_decay]
    type = FVReaction
    variable = 'Y91_f'
    rate = ${lambda_Y91}
    block = ${fluid_blocks}
  []
  [Y91_grow]
    type = FVCoupledForce
    variable = 'Y91_f'
    v = 'Sr91_f'
    coef = ${lambda_Br91}
    block = ${fluid_blocks}
  []
  [Y91_fg]
    type = NSFVMixturePhaseInterface
    variable = 'Y91_f'
    phase_coupled = 'Y91_g'
    alpha = 'hlg_Y91'
  []
  [Y91_gf]
    type = NSFVMixturePhaseInterface
    variable = 'Y91_g'
    phase_coupled = 'Y91_f'
    alpha = 'hlg_Y91'
  []
  [Y91_g_advection]
    type = PINSFVMassAdvection
    variable = 'Y91_g'
    rho = 'Y91_g'
    block = ${fluid_blocks}
  []
  [Y91_g_diffusion]
    type = FVDiffusion
    variable = 'Y91_g'
    coeff = ${D_mu_g}
    block = ${fluid_blocks}
  []
  [Y91_g_decay]
    type = FVReaction
    variable = 'Y91_g'
    rate = ${lambda_Y91}
    block = ${fluid_blocks}
  []
  [Y91_grow_g]
    type = FVCoupledForce
    variable = 'Y91_g'
    v = 'Sr91_g'
    coef = ${lambda_Br91}
    block = ${fluid_blocks}
  []

  [Zr91_time]
    type = FVFunctorTimeKernel
    variable = 'Zr91_f'
    block = ${fluid_blocks}
  []
  [Zr91_advection]
    type = PINSFVMassAdvection
    variable = 'Zr91_f'
    rho = 'Zr91_porous'
    block = ${fluid_blocks}
  []
  [Zr91_diffusion]
    type = FVDiffusion
    variable = 'Zr91_f'
    coeff = ${D_mu}
    block = ${fluid_blocks}
  []
  [Zr91_turb_diffusion]
    type = INSFVMixingLengthScalarDiffusion
    variable = 'Zr91_f'
    schmidt_number = ${Sc_t}
    block = ${fluid_blocks}
  []
  [Zr91_src]
    type = FVCoupledForce
    variable = 'Zr91_f'
    v = fission_source
    coef = ${yield_Zr91}
    block = ${fluid_blocks}
  []
  [Zr91_decay]
    type = FVReaction
    variable = 'Zr91_f'
    rate = ${lambda_Zr91}
    block = ${fluid_blocks}
  []
  [Zr91_grow]
    type = FVCoupledForce
    variable = 'Zr91_f'
    v = 'Y91_f'
    coef = ${lambda_Br91}
    block = ${fluid_blocks}
  []
  [Zr91_fg]
    type = NSFVMixturePhaseInterface
    variable = 'Zr91_f'
    phase_coupled = 'Zr91_g'
    alpha = 'hlg_Zr91'
  []
  [Zr91_gf]
    type = NSFVMixturePhaseInterface
    variable = 'Zr91_g'
    phase_coupled = 'Zr91_f'
    alpha = 'hlg_Zr91'
  []
  [Zr91_g_advection]
    type = PINSFVMassAdvection
    variable = 'Zr91_g'
    rho = 'Zr91_g'
    block = ${fluid_blocks}
  []
  [Zr91_g_diffusion]
    type = FVDiffusion
    variable = 'Zr91_g'
    coeff = ${D_mu_g}
    block = ${fluid_blocks}
  []
  [Zr91_g_decay]
    type = FVReaction
    variable = 'Zr91_g'
    rate = ${lambda_Zr91}
    block = ${fluid_blocks}
  []
  [Zr91_grow_g]
    type = FVCoupledForce
    variable = 'Zr91_g'
    v = 'Y91_g'
    coef = ${lambda_Br91}
    block = ${fluid_blocks}
  []
[]

[FVBCs]
  [no-slip-u]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom right loop_boundary'
    variable = superficial_vel_x
    function = 0
  []
  [no-slip-v]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom right loop_boundary'
    variable = superficial_vel_y
    function = 0
  []

  [sym_Br91_f]
    type = INSFVSymmetryScalarBC
    boundary = 'left'
    variable = 'Br91_f'
  []
  [sym_Br91_g]
    type = INSFVSymmetryScalarBC
    boundary = 'left'
    variable = 'Br91_g'
  []
  [Br91_g_top]
    type = FVDirichletBC
    boundary = 'top'
    variable = 'Br91_g'
    value = 0.0
  []

  [sym_Kr91_f]
    type = INSFVSymmetryScalarBC
    boundary = 'left'
    variable = 'Kr91_f'
  []
  [sym_Kr91_g]
    type = INSFVSymmetryScalarBC
    boundary = 'left'
    variable = 'Kr91_g'
  []
  [Kr91_g_top]
    type = FVDirichletBC
    boundary = 'top'
    variable = 'Kr91_g'
    value = 0.0
  []

  [sym_Rb91_f]
    type = INSFVSymmetryScalarBC
    boundary = 'left'
    variable = 'Rb91_f'
  []
  [sym_Rb91_g]
    type = INSFVSymmetryScalarBC
    boundary = 'left'
    variable = 'Rb91_g'
  []
  [Rb91_g_top]
    type = FVDirichletBC
    boundary = 'top'
    variable = 'Rb91_g'
    value = 0.0
  []

  [sym_Sr91_f]
    type = INSFVSymmetryScalarBC
    boundary = 'left'
    variable = 'Sr91_f'
  []
  [sym_Sr91_g]
    type = INSFVSymmetryScalarBC
    boundary = 'left'
    variable = 'Sr91_g'
  []
  [Sr91_g_top]
    type = FVDirichletBC
    boundary = 'top'
    variable = 'Sr91_g'
    value = 0.0
  []

  [sym_Y91_f]
    type = INSFVSymmetryScalarBC
    boundary = 'left'
    variable = 'Y91_f'
  []
  [sym_Y91_g]
    type = INSFVSymmetryScalarBC
    boundary = 'left'
    variable = 'Y91_g'
  []
  [Y91_g_top]
    type = FVDirichletBC
    boundary = 'top'
    variable = 'Y91_g'
    value = 0.0
  []

  [sym_Zr91_f]
    type = INSFVSymmetryScalarBC
    boundary = 'left'
    variable = 'Zr91_f'
  []
  [sym_Zr91_g]
    type = INSFVSymmetryScalarBC
    boundary = 'left'
    variable = 'Zr91_g'
  []
  [Zr91_g_top]
    type = FVDirichletBC
    boundary = 'top'
    variable = 'Zr91_g'
    value = 0.0
  []
[]

# ==============================================================================
# AUXVARIABLES AND AUXKERNELS
# ==============================================================================

[AuxVariables]
  [fission_source]
    type = MooseVariableFVReal
    initial_condition = 1e-6
  []
  [mixing_length]
    type = MooseVariableFVReal
    initial_condition = 1e-6
  []
  [superficial_vel_x]
    type = PINSFVSuperficialVelocityVariable
    initial_condition = 1e-6
  []
  [superficial_vel_y]
    type = PINSFVSuperficialVelocityVariable
    initial_condition = 1e-6
  []
  [pressure]
    type = INSFVPressureVariable
    initial_condition = 1.0
  []
  [a_u_var]
    type = MooseVariableFVReal
    initial_condition = 1e-6
  []
  [a_v_var]
    type = MooseVariableFVReal
    initial_condition = 1e-6
  []
  [porosity_var]
    type = MooseVariableFVReal
    initial_condition = 1.0
  []
  [T_fluid]
    type = MooseVariableFVReal
    initial_condition = 1.0
  []
  [void_f]
    type = MooseVariableFVReal
  []
  [interface_area_concentration_var]
    type = MooseVariableFVReal
    initial_condition = 1.0
  []
[]

[AuxKernels]
  [interface_area_concentration_aux]
    type = FunctorAux
    variable = interface_area_concentration_var
    functor = 'interface_area_concentration'
    execute_on = 'timestep_end'
    block = ${fluid_blocks}
  []
[]

################################################################################
# MATERIALS
################################################################################

[FunctorMaterials]
  [Br91_mat]
    type = ADParsedFunctorMaterial
    expression = 'Br91_f / porosity_var'
    functor_names = 'Br91_f porosity_var'
    functor_symbols = 'Br91_f porosity_var'
    property_name = 'Br91_porous'
  []
  [Kr91_mat]
    type = ADParsedFunctorMaterial
    expression = 'Kr91_f / porosity_var'
    functor_names = 'Kr91_f porosity_var'
    functor_symbols = 'Kr91_f porosity_var'
    property_name = 'Kr91_porous'
  []
  [Rb91_mat]
    type = ADParsedFunctorMaterial
    expression = 'Rb91_f / porosity_var'
    functor_names = 'Rb91_f porosity_var'
    functor_symbols = 'Rb91_f porosity_var'
    property_name = 'Rb91_porous'
  []
  [Sr91_mat]
    type = ADParsedFunctorMaterial
    expression = 'Sr91_f / porosity_var'
    functor_names = 'Sr91_f porosity_var'
    functor_symbols = 'Sr91_f porosity_var'
    property_name = 'Sr91_porous'
  []
  [Y91_mat]
    type = ADParsedFunctorMaterial
    expression = 'Y91_f / porosity_var'
    functor_names = 'Y91_f porosity_var'
    functor_symbols = 'Y91_f porosity_var'
    property_name = 'Y91_porous'
  []
  [Zr91_mat]
    type = ADParsedFunctorMaterial
    expression = 'Zr91_f / porosity_var'
    functor_names = 'Zr91_f porosity_var'
    functor_symbols = 'Zr91_f porosity_var'
    property_name = 'Zr91_porous'
  []

  [interface_area_concentration]
    type = ADParsedFunctorMaterial
    expression = 'void_f * 6 / ${dp}'
    functor_names = 'void_f'
    functor_symbols = 'void_f'
    property_name = 'interface_area_concentration'
  []

  [hlg_Br91]
    type = ADParsedFunctorMaterial
    expression = 'interface_area_concentration * ${Br91_alpha_fg}'
    functor_names = 'interface_area_concentration'
    functor_symbols = 'interface_area_concentration'
    property_name = 'hlg_Br91'
  []
  [hlg_Kr91]
    type = ADParsedFunctorMaterial
    expression = 'interface_area_concentration * ${Kr91_alpha_fg}'
    functor_names = 'interface_area_concentration'
    functor_symbols = 'interface_area_concentration'
    property_name = 'hlg_Kr91'
  []
  [hlg_Rb91]
    type = ADParsedFunctorMaterial
    expression = 'interface_area_concentration * ${Rb91_alpha_fg}'
    functor_names = 'interface_area_concentration'
    functor_symbols = 'interface_area_concentration'
    property_name = 'hlg_Rb91'
  []
  [hlg_Sr91]
    type = ADParsedFunctorMaterial
    expression = 'interface_area_concentration * ${Sr91_alpha_fg}'
    functor_names = 'interface_area_concentration'
    functor_symbols = 'interface_area_concentration'
    property_name = 'hlg_Sr91'
  []
  [hlg_Y91]
    type = ADParsedFunctorMaterial
    expression = 'interface_area_concentration * ${Y91_alpha_fg}'
    functor_names = 'interface_area_concentration'
    functor_symbols = 'interface_area_concentration'
    property_name = 'hlg_Y91'
  []
  [hlg_Zr91]
    type = ADParsedFunctorMaterial
    expression = 'interface_area_concentration * ${Zr91_alpha_fg}'
    functor_names = 'interface_area_concentration'
    functor_symbols = 'interface_area_concentration'
    property_name = 'hlg_Zr91'
  []
[]
# ==============================================================================
# POSTPROCESSORS
# ==============================================================================
[Postprocessors]
  [fission_source_integral]
    type = ElementAverageValue
    variable = fission_source
  []
[]
################################################################################
# EXECUTION / SOLVE
################################################################################

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -sub_pc_factor_shift_type'
  petsc_options_value = ' lu       NONZERO'
  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  nl_rel_tol = 1e-6
  #line_search = l2
  nl_max_its = 100

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    optimal_iterations = 20
    iteration_window = 2
    growth_factor = 2
    cutback_factor = 0.5
  []

  end_time = 1e10
  steady_state_detection = true
[]

[Debug]
  show_var_residual_norms = true
[]

################################################################################
# SIMULATION OUTPUTS
################################################################################

[Outputs]
  csv = true
  exodus = true
  print_linear_converged_reason = false
  print_linear_residuals = false
  print_nonlinear_converged_reason = false
[]
