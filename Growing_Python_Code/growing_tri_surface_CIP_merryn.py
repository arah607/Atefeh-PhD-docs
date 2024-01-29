#!/usr/bin/env python
import os
import numpy as np

from aether.diagnostics import set_diagnostics_on
from aether.indices import define_problem_type, get_ne_radius
from aether.geometry import define_node_geometry, define_1d_elements, define_rad_from_file, define_rad_from_geom
from aether.geometry import evaluate_ordering, import_ply_triangles, make_data_grid
from aether.growtree import grow_tree, smooth_1d_tree
from aether.exports import export_1d_elem_geometry, export_node_geometry, export_elem_field, export_data_geometry
from aether.exports import export_1d_elem_field, export_triangle_elements, export_triangle_nodes

########################################################################################################

def main():
    
    set_diagnostics_on(False)
    define_problem_type('grow_tree')  # sets up the array indices for a 1D geometry

    lobe_parameters = { 'RUL' : [ 500, 4.88 ], 'RML' : [ 525, 4.93], 'RLL' : [533, 5.15], 'LUL' : ['LULA', 5.0], 'LLL' : ['LLLA', 5.0] }
    lobes = ['RUL', 'RML', 'RLL']

    angle_max = 60.0         # maximum allowed branching angle (to parent direction)
    angle_min = 20.0         # minimum allowed branching angle (to parent direction)
    branch_fraction = 0.4    # fraction of distance to COFM to branch
    length_limit = 1.0       # limit on terminal branches
    min_length = 1.2         # minimum length of elements
    rotation_limit = 180.0   # limit on angle of rotation between planes (not working)
    peel = 5.0               # proportional gap between surface and data grid (%scaling, not distance)

    subject = 'Graph_15814W_G2_testtt'
    input_directory = '/home/arah607/Documents/working_demos/check_node_inside_surface/Data/COPD_gene/test_merryn_files'
    output_directory = '/home/arah607/Documents/working_demos/check_node_inside_surface/Data/COPD_gene/test_merryn_files'
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    template = os.path.join(input_directory, subject )
    define_node_geometry(template)
    print ('Nodes read')
    define_1d_elements(template)
    print ('Elements read')
    define_rad_from_file(template)
    print ('Radii read')
    filename = 'CIP'
    group_name = '1d_tree'
    export_node_geometry(os.path.join(output_directory, filename), group_name)
    export_1d_elem_geometry(os.path.join(output_directory, filename), group_name)
    ne_radius = get_ne_radius()
    field_name = 'radius'
    export_1d_elem_field(ne_radius, os.path.join(output_directory, filename + '_radius'), group_name, field_name)
    print ('CIP tree exported')

    for lobe in lobes:
        ply_name = lobe + '_surf'
        group_name = ply_name
        vessel = lobe_parameters[lobe][0]
        map_name = lobe + '_mapping'
        spacing = lobe_parameters[lobe][1]
        print ("Growing tree for", lobe, " from element", vessel, "with spacing", spacing)
        import_ply_triangles(os.path.join(input_directory, ply_name))
        print ('Surfaces imported')
        export_triangle_nodes(os.path.join(output_directory, ply_name), group_name)
        export_triangle_elements(os.path.join(output_directory, ply_name), group_name)
        print ('Triangles exported, now going to make the grid....')
        make_data_grid( [0], peel, spacing, os.path.join(output_directory, ply_name), group_name)
        print ('Grid made')
        export_data_geometry( os.path.join(output_directory, ply_name), group_name, 0)
        print ('Data exported, now going to grow the tree .... ')
        grow_tree( [0], vessel, angle_max, angle_min, branch_fraction, length_limit, min_length, rotation_limit, True, map_name)
        print (lobe + ' tree grown')
    smooth_1d_tree( 1, length_limit)
    print ('Tree smoothed')

    filename = 'grown'
    group_name = '1d_tree'
    export_node_geometry(os.path.join(output_directory, filename), group_name)
    export_1d_elem_geometry(os.path.join(output_directory, filename), group_name)
    print ('Tree exported')

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
