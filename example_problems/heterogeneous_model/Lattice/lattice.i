# key lattice parameter before applying volume preserve
# There can be a way to provide a universal conversion function of
# geting the l_fuel and r_fuel after applying volume preserve
# but for now, let me try an easy way first

l_lat           = 5.08  # length of a lattice edge
r_fuel_original = 0.508 # radius of the round end
r_fuel          = 0.535385 # radius of the round end after applying the volume preservce with Th_div=2
l_fuel_original = 3.048 # length of the fuel channel (including round end)
l_fuel          = ${fparse l_fuel_original-2*r_fuel_original+ 2* r_fuel} # ength of the fuel channel (including round end) after applying the volume preservce with Th_div=2

Th_div   = 2

[Mesh]
  # final_generator = left_part_3
  [Corner_bottom_left_gen_2]
    type = ConcentricCircleMeshGenerator
    num_sectors = ${Th_div}
    radii = ${r_fuel}
    rings = '1 1'
    has_outer_square = on
    pitch = 2
    portion = top_right
    preserve_volumes = off
    # smoothing_max_it = 3
  []
  [Corner_bottom_left_2]
    type = BlockDeletionGenerator
    input = Corner_bottom_left_gen_2
    block = '2'
    new_boundary = '333'
  []
  [Corner_bottom_left_3]
    type = RenameBoundaryGenerator
    input = Corner_bottom_left_2
    old_boundary = '1 2'
    new_boundary = '600 1503'
  []
  ####$$$$####
  [Corner_top_right_gen_3]
    type = ConcentricCircleMeshGenerator
    num_sectors = ${Th_div}
    radii = ${r_fuel}
    rings = '1 1'
    has_outer_square = on
    pitch = 2
    portion = bottom_left
    preserve_volumes = off
    #smoothing_max_it = 3
  []
  [Corner_top_right_del_3]
    type = BlockDeletionGenerator
    input = Corner_top_right_gen_3
    block = '2'
    new_boundary = '444'
  []
  [Corner_top_right_del_4]
    type = RenameBoundaryGenerator
    input = Corner_top_right_del_3
    old_boundary = '1 2'
    new_boundary = '1502 600'
  []
  # [Corner_top_right_del_5]
  #   type = RenameBoundaryGenerator
  #   input = Corner_top_right_del_4
  #   old_boundary = '600 1502'
  #   new_boundary = 'cleft cright'
  # []
  [Corner_top_right_3]
    type = TransformGenerator
    transform = TRANSLATE
    vector_value = '${fparse (l_lat-l_fuel_original)/2+r_fuel_original} ${fparse (l_lat-l_fuel_original)/2+r_fuel_original}  0'
    input = Corner_top_right_del_4
  []
  ####$$$$####
  [connect_two_circles]
    type = FillBetweenSidesetsGenerator
    input_mesh_1 = 'Corner_bottom_left_3'
    input_mesh_2 = 'Corner_top_right_3'
    boundary_1 = '333'
    boundary_2 = '444'
    begin_side_boundary_id = 20001
    end_side_boundary_id = 20002
    num_layers = ${Th_div}
    keep_inputs = true
    use_quad_elements = true
    block_id = 200
  []
  ####################################################
  [Corner_top_left_gen_4]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 4
    num_sectors_per_side = '${Th_div} ${Th_div} ${Th_div} ${Th_div}'
    background_intervals = 1
    background_block_ids = 200
    polygon_size = ${fparse sin(45/180*pi)*(l_lat-l_fuel)/2}
    quad_center_elements=true
    preserve_volumes = on
  []
  [Corner_top_left_gen_del_pre_4]
    type = PlaneDeletionGenerator
    point = '0 0 0'
    normal = '-1 0 0'
    input = Corner_top_left_gen_4
    new_boundary = 400
  []
  [Corner_top_left_gen_del_pre_2_4]
    type = PlaneDeletionGenerator
    point = '0 0 0'
    normal = '0 1 0'
    input = Corner_top_left_gen_del_pre_4
    new_boundary = 400
  []
  [Corner_top_left_gen_del_pre_2_5]
    type = BoundaryDeletionGenerator
    input = Corner_top_left_gen_del_pre_2_4
    boundary_names = '1'    
  []
  [Corner_top_left_4]
    type = TransformGenerator
    transform = TRANSLATE
    vector_value = '0 ${fparse (l_lat-l_fuel_original)/2+r_fuel_original}  0'
    input = Corner_top_left_gen_del_pre_2_5
  []
  [Left_corner_pre]
    type = StitchedMeshGenerator
    inputs = 'connect_two_circles Corner_top_left_4'
    clear_stitched_boundary_ids = true
    stitch_boundaries_pairs = '20002 10000'
  []
  [Left_corner_pre_1]
    type = RenameBoundaryGenerator
    input = Left_corner_pre
    old_boundary = '400'
    new_boundary = '600'
  []
  [Left_combine_pre]
    type = TransformGenerator
    transform = ROTATE
    vector_value ='0 0 90'
    input = Left_corner_pre_1
  []
  [Left_combine]
    type = TransformGenerator
    transform = TRANSLATE
    vector_value = '${fparse (l_lat-l_fuel_original)/2+r_fuel_original} -${fparse (l_lat+l_fuel)/2-r_fuel}  0'
   # vector_value = '${fparse (l_lat-l_fuel_original)/2+r_fuel_original} -${fparse (l_fuel_original+r_fuel_original)}  0'   
    input = Left_combine_pre
  []
  [left_part_1]
    type = FillBetweenSidesetsGenerator
    input_mesh_1 = 'Left_corner_pre_1'
    input_mesh_2 = 'Left_combine'
    boundary_1 = '1503'
    boundary_2 = '1502'
    begin_side_boundary_id = 400
    end_side_boundary_id = 20000
    num_layers = ${Th_div}
    keep_inputs = true
    use_quad_elements = true
    block_id = '1'
  []
  [left_part_2]
    type = FillBetweenSidesetsGenerator
    input_mesh_1 = 'Left_corner_pre_1'
    input_mesh_2 = 'Left_combine'
    boundary_1 = '20001'
    boundary_2 = '20001'
    begin_side_boundary_id = 20000
    end_side_boundary_id = 20000
    input_boundary_1_id = 20000
    input_boundary_2_id = 20000
    num_layers = ${Th_div}
    keep_inputs = false
    use_quad_elements = true
    block_id = 200
  []
  [rename_left_part_1]
    type = RenameBoundaryGenerator
    input = left_part_1
    old_boundary = 20001
    new_boundary = 20000
  []
  [left_part_3]
    type = StitchedMeshGenerator
    inputs = 'rename_left_part_1 left_part_2'
    clear_stitched_boundary_ids = true
    stitch_boundaries_pairs = '20000 20000'
  []

######################################################
  [right_mirror_1]
    type = TransformGenerator
    transform = ROTATE
    vector_value ='180 0 0'
    input = left_part_3
  []
   [right_mirror_2]
    type = TransformGenerator
    transform = TRANSLATE
    vector_value = '${fparse l_lat} ${fparse -(l_lat-l_fuel_original)/2-2*r_fuel_original}  0'
    input = right_mirror_1
  []
  [Final_connection_1]
    type = FillBetweenSidesetsGenerator
    input_mesh_1 = 'left_part_3'
    input_mesh_2 = 'right_mirror_2'
    boundary_1 = '1503'
    boundary_2 = '1502'
    begin_side_boundary_id = 400
    end_side_boundary_id = 20001
    num_layers = ${Th_div}
    keep_inputs = true
    use_quad_elements = true
    block_id = 1
  []
  [Final_connection_2]
    type = FillBetweenSidesetsGenerator
    input_mesh_1 = 'left_part_3'
    input_mesh_2 = 'right_mirror_2'
    boundary_1 = '20001'
    boundary_2 = '20001'
    begin_side_boundary_id = 20001
    end_side_boundary_id = 20001
    input_boundary_1_id = 20001
    input_boundary_2_id = 20001
    num_layers = ${Th_div}
    keep_inputs = false
    use_quad_elements = true
    block_id = 200
  []
  [Final_connection_3]
    type = StitchedMeshGenerator
    inputs = 'Final_connection_1 Final_connection_2'
    clear_stitched_boundary_ids = true
    stitch_boundaries_pairs = '20001 20001'
  []
  [Final_connection_4]
    type = FillBetweenSidesetsGenerator
    input_mesh_1 = 'left_part_3'
    input_mesh_2 = 'right_mirror_2'
    boundary_1 = '1502'
    boundary_2 = '1503'
    input_boundary_1_id = 20002
    input_boundary_2_id = 20002
    begin_side_boundary_id = 400
    end_side_boundary_id = 20002
    num_layers = ${Th_div}
    keep_inputs = false
    use_quad_elements = true
    block_id = 1
  []
  [rename_Final_connection_3]
    type = RenameBoundaryGenerator
    input = Final_connection_3
    old_boundary = '1503 1502'
    new_boundary = '20002 20002'
  []
  [Final_connection_6]
    type = StitchedMeshGenerator
    inputs = 'rename_Final_connection_3 Final_connection_4'
    clear_stitched_boundary_ids = true
    stitch_boundaries_pairs = '20002 20002'
  []
  [lattice]
    type = RenameBoundaryGenerator
    input = Final_connection_6
    old_boundary = '20003 400'
    new_boundary = '600 600'
  []
  [extrude]
     type = AdvancedExtruderGenerator
     input = lattice
     heights = '170.0'
     num_layers = '17 '
     direction = '0 0 1'
     top_boundary = 2000
     bottom_boundary = 2001
     subdomain_swaps = '1 1 200 200'
  []
  [scale]
    type = TransformGenerator
    input = extrude
    transform = SCALE
    vector_value = '1e-2 1e-2 1e-2'
  []
[]

[Variables]
  [./u]
  [../]
[]

[Kernels]
  [./dt]
    type = TimeDerivative
    variable = u
  [../]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[BCs]
  [./outer]
    type = DirichletBC
    variable = u
    boundary = '10000'
    value = 1
  [../]
  [./inner]
    type = DirichletBC
    variable = u
    boundary = '2000'
    value = 0
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  exodus = true
[]