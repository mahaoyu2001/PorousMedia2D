"""Export element stress and area data from completed Abaqus ODB files.

Run this script from the repository root with the Abaqus Python interpreter::

    abaqus python "python script/PostProcessing/ExtractLE_Vol.py"

Each output row contains ``S11 S22 S12 EVOL`` for one CPE3 element. The
MATLAB scripts use EVOL as the element-area weight. This file deliberately
uses Python 2.7-compatible syntax because Abaqus 2022 ships an older Python
interpreter on many installations.
"""
from __future__ import print_function

import os
import sys
import time

import numpy as np
import odbAccess


# ============================= USER SETTINGS =============================
# Run one base name at a time, or list several completed analyses here.
MODEL_NAME_LIST = ['2DVoid_Vol5']
LABEL_LIST = ['a', 'b', 'c', 'd', 'e']

INSTANCE_NAME = 'COMPOSITE-CUT'
STEP_NAME = 'Step-1'
FIELD_NAME = 'S'
STRESS_COMPONENT_IDS = (0, 1, 3)  # Abaqus order -> S11, S22, S12 in 2-D.

# ``S2`` is retained for compatibility with the MATLAB plotting scripts:
# frame 1 is the initial frame and frame 2 is the final static-load frame.
OUTPUT_TAG = 'S2'
# ========================================================================


def locate_project_directory():
    """Return the repository root in GUI and no-GUI Abaqus sessions."""
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        script_argument = sys.argv[0] if len(sys.argv) else ''
        if script_argument.lower().endswith('.py'):
            script_dir = os.path.dirname(os.path.abspath(script_argument))
        else:
            script_dir = os.path.abspath(os.path.join(
                os.getcwd(), 'python script', 'PostProcessing'))
    return os.path.abspath(os.path.join(script_dir, '..', '..'))


def extract_component_array(field_outputs, region, field_name, component_ids):
    """Return selected components from an element integration-point field."""
    values = field_outputs[field_name].getSubset(region=region).values
    rows = []
    for value in values:
        rows.append([value.data[index] for index in component_ids])
    return np.asarray(rows, dtype=float)


def extract_scalar_array(field_outputs, region, field_name):
    """Return one scalar value per element integration point."""
    values = field_outputs[field_name].getSubset(region=region).values
    return np.asarray([value.data for value in values], dtype=float)


def export_one_odb(project_dir, output_dir, job_name):
    """Export the final frame of one ODB and return the output file path."""
    odb_path = os.path.join(project_dir, job_name + '.odb')
    if not os.path.isfile(odb_path):
        raise IOError("ODB file not found: '{}'".format(odb_path))

    odb = odbAccess.openOdb(path=odb_path, readOnly=True)
    try:
        if STEP_NAME not in odb.steps.keys():
            raise KeyError("Step '{}' not found in '{}'".format(
                STEP_NAME, odb_path))
        step = odb.steps[STEP_NAME]
        if len(step.frames) < 2:
            raise ValueError("Step '{}' in '{}' has no result frame".format(
                STEP_NAME, odb_path))

        instance_key = INSTANCE_NAME.upper()
        if instance_key not in odb.rootAssembly.instances.keys():
            raise KeyError("Instance '{}' not found in '{}'".format(
                instance_key, odb_path))
        instance = odb.rootAssembly.instances[instance_key]
        fields = step.frames[-1].fieldOutputs

        stress = extract_component_array(
            fields, instance, FIELD_NAME, STRESS_COMPONENT_IDS)
        element_area = extract_scalar_array(fields, instance, 'EVOL')
        if stress.shape[0] != element_area.shape[0]:
            raise ValueError(
                "S and EVOL contain different numbers of values in '{}'"
                .format(odb_path))

        data = np.column_stack((stress, element_area))
        output_path = os.path.join(
            output_dir, '{}_{}.txt'.format(job_name, OUTPUT_TAG))
        np.savetxt(output_path, data, fmt='%.10e', delimiter=' ',
                   header='S11 S22 S12 EVOL')
        return output_path
    finally:
        odb.close()


def main():
    started = time.time()
    project_dir = locate_project_directory()
    output_dir = os.path.join(project_dir, 'output', 'FieldValue')
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    exported = 0
    for base_name in MODEL_NAME_LIST:
        for label in LABEL_LIST:
            job_name = '{}_{}'.format(base_name, label)
            output_path = export_one_odb(project_dir, output_dir, job_name)
            exported += 1
            print('Exported: {}'.format(output_path))

    print('Finished: {} file(s), {:.2f} s'.format(
        exported, time.time() - started))


if __name__ == '__main__':
    main()
