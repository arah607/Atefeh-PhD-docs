
from math import cos
from math import floor
import perfusion_Clark2011 as per_clark
import os
import shutil
from aether.diagnostics import set_diagnostics_on
from aether.indices import perfusion_indices, get_ne_radius
from aether.filenames import read_geometry_main, get_filename
from aether.geometry import append_units, define_node_geometry, define_1d_elements, define_rad_from_geom, add_matching_mesh
from aether.exports import export_1d_elem_geometry, export_node_geometry, export_1d_elem_field, export_node_field, export_terminal_perfusion
from aether.pressure_resistance_flow import evaluate_prq
#Import numpy
import numpy as np
import matplotlib.pyplot as plt
from switch import Switch
from math import pi
from math import exp
from math import sqrt
#import importlib

if __name__ == '__main__':

    for i in range(3):
        # P_rv = 11929.443538438134
        # P_la = 202.33952809536052
        per_clark.perfusion(i)

