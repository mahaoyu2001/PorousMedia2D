"""Create 2-D perforated RVE geometry from pore-centre coordinate files.

This stage is executed by ``Run2D_Geom_Batch.py``. Each input text file has
two whitespace-separated columns (x and y), one circular-pore centre per row.
Coordinates and VOID_RADIUS use the same consistent unit system as RVE_SIZE.
"""
from abaqus import *
from abaqusConstants import *
import numpy as np
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
RVE_SIZE = option('RVE_SIZE', (1.0, 1.0))
VOID_RADIUS = option('VOID_RADIUS', 0.005)
MATERIAL_FILE = option('MATERIAL_FILE', 'MaterialData/Ceramic.txt')
MODEL_DATA_DIR = option('MODEL_DATA_DIR', 'ModelData')
REPLACE_EXISTING_MODELS = option('REPLACE_EXISTING_MODELS', False)


def target_names():
    return [base + SUFFIX + label for label in LABEL_LIST
            for base in MODEL_NAME_LIST]


def read_material(path):
    """Read ``key value`` material data and validate required properties."""
    material = {}
    with open(path, 'r') as material_file:
        for line in material_file:
            fields = line.split()
            if len(fields) == 2:
                material[fields[0]] = float(fields[1])
    for required in ('Density', 'E', 'v'):
        if required not in material:
            raise ValueError("Material file '{}' lacks '{}'".format(path, required))
    return material


def reference_point_after_creation(assembly, point):
    """Create a reference point and return its repository object."""
    assembly.ReferencePoint(point=point)
    return assembly.referencePoints[max(assembly.referencePoints.keys())]


def create_model(model_name):
    """Create one material/geometry model named after its coordinate file."""
    started = time.time()
    if model_name in mdb.models.keys():
        if not REPLACE_EXISTING_MODELS:
            raise ValueError("Model '{}' already exists. Set "
                             "REPLACE_EXISTING_MODELS=True to replace it."
                             .format(model_name))
        del mdb.models[model_name]

    l1, l2 = RVE_SIZE
    material_data = read_material(MATERIAL_FILE)
    model = mdb.Model(name=model_name)
    matrix_material = model.Material(name='Material-matrix')
    matrix_material.Elastic(table=((material_data['E'], material_data['v']),))
    matrix_material.Density(table=((material_data['Density'],),))
    model.HomogeneousSolidSection(name='Section-matrix',
                                  material='Material-matrix', thickness=None)

    sketch = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    sketch.setPrimaryObject(option=STANDALONE)
    sketch.rectangle(point1=(0.0, 0.0), point2=(l1, l2))
    matrix = model.Part(name='Matrix', dimensionality=TWO_D_PLANAR,
                        type=DEFORMABLE_BODY)
    matrix.BaseShell(sketch=sketch)
    sketch.unsetPrimaryObject()
    del model.sketches['__profile__']
    matrix.Set(faces=matrix.faces, name='Set-Matrix')
    matrix.SectionAssignment(region=matrix.sets['Set-Matrix'],
                             sectionName='Section-matrix', offset=0.0,
                             offsetType=MIDDLE_SURFACE, offsetField='',
                             thicknessAssignment=FROM_SECTION)

    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)
    assembly.Instance(name='Matrix', part=matrix, dependent=ON)

    data_file = MODEL_DATA_DIR + '/' + model_name + '.txt'
    points = np.atleast_2d(np.loadtxt(data_file, dtype=float))
    if points.shape[1] != 2:
        raise ValueError("Point file '{}' must contain two columns.".format(data_file))
    particle_instances = []
    for index, point in enumerate(points):
        sketch = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        sketch.setPrimaryObject(option=STANDALONE)
        sketch.CircleByCenterPerimeter(center=(0.0, 0.0),
                                       point1=(VOID_RADIUS, 0.0))
        particle_name = 'particle-{}'.format(index + 1)
        particle = model.Part(name=particle_name, dimensionality=TWO_D_PLANAR,
                              type=DEFORMABLE_BODY)
        particle.BaseShell(sketch=sketch)
        sketch.unsetPrimaryObject()
        del model.sketches['__profile__']
        instance_name = 'particle-1-{}'.format(index + 1)
        assembly.Instance(name=instance_name, part=particle, dependent=ON)
        assembly.translate(instanceList=(instance_name,),
                           vector=(point[0], point[1], 0.0))
        particle_instances.append(assembly.instances[instance_name])

    if not particle_instances:
        raise ValueError("Point file '{}' contains no void centers.".format(data_file))
    assembly.InstanceFromBooleanMerge(name='AllParticles',
                                      instances=tuple(particle_instances),
                                      keepIntersections=ON,
                                      originalInstances=DELETE, domain=GEOMETRY)
    assembly.InstanceFromBooleanCut(name=PART_NAME,
                                    instanceToBeCut=assembly.instances['Matrix'],
                                    cuttingInstances=(assembly.instances['AllParticles-1'],),
                                    originalInstances=SUPPRESS)
    assembly.features.changeKey(fromName=PART_NAME + '-1', toName=PART_NAME)

    rp0 = reference_point_after_creation(assembly, (0.0, 0.0, 0.0))
    rp1 = reference_point_after_creation(assembly, (l1, l2, 0.0))
    assembly.Set(referencePoints=(rp0,), name='Set-rp0')
    assembly.Set(referencePoints=(rp1,), name='Set-rp1')

    final_part = model.parts[PART_NAME]
    final_instance_name = PART_NAME
    for part_name in list(model.parts.keys()):
        if part_name != PART_NAME:
            del model.parts[part_name]
    for instance_name in list(assembly.instances.keys()):
        if instance_name != final_instance_name:
            try:
                assembly.deleteFeatures((instance_name,))
            except Exception:
                pass

    void_fraction = 1.0 - final_part.getArea(final_part.faces) / (l1 * l2)
    print('Created {}: void fraction = {:.6g}; time = {:.2f} s'.format(
        model_name, void_fraction, time.time() - started))


for _model_name in target_names():
    create_model(_model_name)
