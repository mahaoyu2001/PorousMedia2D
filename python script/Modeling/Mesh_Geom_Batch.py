"""Mesh models produced by ``Matrix2D_Geom_Batch.py`` with CPE3 elements."""
from abaqus import *
from abaqusConstants import *
import mesh
import regionToolset
import time

try:
    RUN_CONFIG
except NameError:
    RUN_CONFIG = {}


def option(name, default):
    return RUN_CONFIG.get(name, default)


PART_NAME = option('PART_NAME', 'Composite-cut')
MODEL_NAME_LIST = option('MODEL_NAME_LIST', ['2DVoid_Vol1_r0'])
LABEL_LIST = option('LABEL_LIST', ['a', 'b', 'c', 'd', 'e'])
SUFFIX = option('SUFFIX', '_')
MESH_SIZE = option('MESH_SIZE', 0.0004)
DEVIATION_FACTOR = option('DEVIATION_FACTOR', 0.1)
MIN_SIZE_FACTOR = option('MIN_SIZE_FACTOR', 0.1)
ELEMENT_CODE = option('ELEMENT_CODE', CPE3)


def target_names():
    return [base + SUFFIX + label for label in LABEL_LIST
            for base in MODEL_NAME_LIST]


for model_name in target_names():
    started = time.time()
    part = mdb.models[model_name].parts[PART_NAME]
    # Deleting any existing mesh makes this stage safe to repeat from the
    # geometry checkpoint.
    part.deleteMesh()
    region = regionToolset.Region(faces=part.faces)
    element_type = mesh.ElemType(elemCode=ELEMENT_CODE, elemLibrary=STANDARD)
    part.setMeshControls(regions=part.faces, elemShape=TRI)
    part.setElementType(regions=region, elemTypes=(element_type,))
    part.seedPart(size=MESH_SIZE, deviationFactor=DEVIATION_FACTOR,
                  minSizeFactor=MIN_SIZE_FACTOR)
    part.generateMesh()
    print('Meshed {}: {} elements; time = {:.2f} s'.format(
        model_name, len(part.elements), time.time() - started))
