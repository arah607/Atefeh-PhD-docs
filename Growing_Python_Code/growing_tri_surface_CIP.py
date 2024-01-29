#!/usr/bin/env python
import os
import numpy as np

from aether.diagnostics import set_diagnostics_on
from aether.indices import define_problem_type, get_ne_radius
from aether.geometry import define_node_geometry, define_1d_elements, define_rad_from_file, define_rad_from_geom
from aether.geometry import import_ply_triangles, make_data_grid
from aether.growtree import grow_tree , smooth_1d_tree
from aether.exports import export_1d_elem_geometry, export_node_geometry, export_elem_field, export_data_geometry
from aether.exports import export_1d_elem_field, export_triangle_elements, export_triangle_nodes

########################################################################################################

def main():

    set_diagnostics_on(False)
    define_problem_type('grow_tree')  # sets up the array indices for a 1D geometry
    list_chck=[]

    # RUL_airway = 15          # 1d element number that supplies the RUL
    RUL_airway = 534#478#515          # 1d element number that supplies the RUL
    # RML_airway = 511#515
    RLL_airway = 510#515
    LUL_airway = 480#515
    LLL_airway = 481#515
    # list_chck.append(RUL_airway)
    # # list_chck.append(RML_airway)
    # list_chck.append(RLL_airway)
    # list_chck.append(LUL_airway)
    # list_chck.append(LLL_airway)


    Status_list = ['RUL_mapping','RML_mapping','RLL_mapping','LUL_mapping','LLL_mapping']

    Status_list1 = ['RUL_surf','RML_surf','RLL_surf','LUL_surf','LLL_surf']

    angle_max = 60.0         # maximuRUL_airwaym allowed branching angle (to parent direction)
    angle_min = 20.0         # minimum allowed branching angle (to parent direction)
    branch_fraction = 0.4    # fraction of distance to COFM to branch
    length_limit = 1.0       # limit on terminal branches
    min_length = 1.2         # minimum length of elements
    rotation_limit = 180.0   # limit on angle of rotation between planes (not working)
    peel = 5.0               # proportional gap between surface and data grid (%scaling, not distance)
    # spacing = 5.48           # distance between grid points in x,y,z directions
    spacing = 5.5           # distance between grid points in x,y,z directions

    # subject = 'P2BRP268-H12816'
    # input_directory = '../geometry'
    subject = 'Graph_15814W_G2_testtt'
    input_directory = '/home/arah607/Documents/working_demos/check_node_inside_surface/Data/COPD_gene/18874D'
    output_directory = '/home/arah607/Documents/working_demos/check_node_inside_surface/Data/COPD_gene/18874D'
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # template = os.path.join(input_directory, subject + '_UpperAirway')
    template = os.path.join(input_directory, subject)
    define_node_geometry(template)
    define_1d_elements(template)
    define_rad_from_file(template)



    ply_name = 'RUL_surf'
    group_name = ply_name
    import_ply_triangles(os.path.join(input_directory, ply_name))
    export_triangle_nodes(os.path.join(output_directory, ply_name), group_name)
    export_triangle_elements(os.path.join(output_directory, ply_name), group_name)

    make_data_grid( [0], peel, spacing, os.path.join(output_directory, ply_name), group_name)
    export_data_geometry( os.path.join(output_directory, ply_name), group_name, 0)
    #for i in range(len(list_chck)):
    grow_tree([0] , RUL_airway, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit, True, 'RUL_mapping')
    smooth_1d_tree( 1, length_limit)

    # ply_name = 'RUL_surf'
    # group_name = ply_name
    # import_ply_triangles(os.path.join(input_directory, ply_name))
    # export_triangle_nodes(os.path.join(output_directory, ply_name), group_name)
    # export_triangle_elements(os.path.join(output_directory, ply_name), group_name)
    # 
    # make_data_grid( [0], peel, spacing, os.path.join(output_directory, ply_name), group_name)
    # export_data_geometry( os.path.join(output_directory, ply_name), group_name, 0)
    # #for i in range(len(list_chck)):
    # grow_tree([0] , RUL_airway, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit, True, 'RUL_mapping')
    # smooth_1d_tree( 1, length_limit)
    # 
    # for i in range(list_chck):
    #     ply_name = Status_list1[i]
    #     group_name = ply_name
    #     import_ply_triangles(os.path.join(input_directory, ply_name))
    #     export_triangle_nodes(os.path.join(output_directory, ply_name), group_name)
    #     export_triangle_elements(os.path.join(output_directory, ply_name), group_name)
    # 
    #     make_data_grid( [0], peel, spacing, os.path.join(output_directory, ply_name), group_name)
    #     export_data_geometry( os.path.join(output_directory, ply_name), group_name, 0)
    #     #for i in range(len(list_chck)):
    #     grow_tree([0] , list_chck[i], angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit, True, Status_list[i])
    #     smooth_1d_tree( 1, length_limit)


    # ply_name = 'RML_surf'
    # group_name = ply_name
    # import_ply_triangles(os.path.join(input_directory, ply_name))
    # export_triangle_nodes(os.path.join(output_directory, ply_name), group_name)
    # export_triangle_elements(os.path.join(output_directory, ply_name), group_name)
    #
    # make_data_grid( [0], peel, spacing, os.path.join(output_directory, ply_name), group_name)
    # export_data_geometry( os.path.join(output_directory, ply_name), group_name, 0)
    # #for i in range(len(list_chck)):
    # grow_tree([0] , LUL_airway, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit, True, 'RUL_mapping')
    # smooth_1d_tree( 1, length_limit)




    #
    #
    # ply_name = 'LUL_surf'
    # group_name = ply_name
    # import_ply_triangles(os.path.join(input_directory, ply_name))
    # export_triangle_nodes(os.path.join(output_directory, ply_name), group_name)
    # export_triangle_elements(os.path.join(output_directory, ply_name), group_name)
    #
    # make_data_grid( [0], peel, spacing, os.path.join(output_directory, ply_name), group_name)
    # export_data_geometry( os.path.join(output_directory, ply_name), group_name, 0)
    # #for i in range(len(list_chck)):
    # grow_tree([0] , LUL_airway, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit, True, 'LUL_mapping')
    # smooth_1d_tree( 1, length_limit)
    #
    #
    # ply_name = 'RML_surf'
    # group_name = ply_name
    # import_ply_triangles(os.path.join(input_directory, ply_name))
    # export_triangle_nodes(os.path.join(output_directory, ply_name), group_name)
    # export_triangle_elements(os.path.join(output_directory, ply_name), group_name)
    #
    # make_data_grid( [0], peel, spacing, os.path.join(output_directory, ply_name), group_name)
    # export_data_geometry( os.path.join(output_directory, ply_name), group_name, 0)
    # #for i in range(len(list_chck)):
    # grow_tree([0] , RML_airway, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit, True, 'RML_mapping')
    # smooth_1d_tree( 1, length_limit)

    filename = 'grown'
    group_name = '1d_tree'
    export_node_geometry(os.path.join(output_directory, filename), group_name)
    export_1d_elem_geometry(os.path.join(output_directory, filename), group_name)

    order_system = 'fit'  # fit the radii between read-in values and min_rad at order 1
    order_options = 'all' # apply to all branches (but will only update ones with radius=0)
    start_at = 'inlet'    # required, but previous option should make obsolete?
    min_rad = 0.2         # radius of order 1 branches
    h_ratio = 0.0         # doesn't matter for the 'fit' option
    define_rad_from_geom(order_system, h_ratio, start_at, min_rad)
    ne_radius = get_ne_radius()
    field_name = 'radius'
    export_1d_elem_field(ne_radius, os.path.join(output_directory, filename + '_radius'), group_name, field_name)


if __name__ == '__main__':
    main()
