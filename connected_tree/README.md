Project Name: Connected Tree

This code aims to build a connected tree for the pulmonary arteries based on the centerline points in the trunk area created by a 3D digitizer and CIP data for one branch.

The code contains two functions:

    buildTree_trunk: This function creates the graph for the trunk and some main pulmonary arteries.
    buildTree_all_branches: This function creates the connected tree based on the CIP data.

Requirements

To run this code, you need to have the following:

    Python 3.x
    Two CSV files as input:
        test_one_branch.csv: This file contains the centerline points in the trunk area created by a 3D digitizer.
        modify_trunk_mainpul_points_27April.csv: This file contains the CIP data for one branch as an example.
