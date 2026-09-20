"""One-command, no-GUI workflow for geometry, mesh, load and analysis.

Edit only the USER SETTINGS block below, then run this file with Abaqus/CAE.
"""
from abaqus import mdb, openMdb
from abaqusConstants import OFF
import os
import sys


# ============================= USER SETTINGS =============================
MODEL_NAME_LIST = ['2DVoid_Vol5']
LABEL_LIST = ['a', 'b', 'c', 'd', 'e']

# Geometry and material
RVE_SIZE = (1.0, 1.0)       # (L1, L2)
VOID_RADIUS = 0.005
MATERIAL_FILE = 'MaterialData/Ceramic.txt'
MODEL_DATA_DIR = 'ModelData'

# Mesh and static analysis
# This is the mesh size used for the reported calculations. Increasing it
# reduces the cost of a workflow smoke test but changes the statistics.
MESH_SIZE = 0.0004
FORCE_Y = 1.0                # Actual total force = FORCE_Y * L1
JOB_CPUS = 4
JOB_DOMAINS = 4

# Execution switches
WRITE_INPUT = True            # Write .inp files before submitting jobs.
SAVE_CAE = True               # Save the complete CAE model database after setup.
CAE_FILE_NAME = None          # None -> '<single ModelName>.cae', including all labels.
OVERWRITE_CAE_FILE = False    # Safety: do not replace an existing .cae by default.
# Recovery settings. ``all`` builds from scratch, ``mesh`` resumes after
# geometry creation, and ``load`` resumes after meshing. Checkpoints group
# all labels belonging to the one base name in MODEL_NAME_LIST.
START_STAGE = 'all'           # Choose: 'all', 'mesh', or 'load'.
SAVE_CHECKPOINTS = True       # Save <ModelName>_geometry.cae and _mesh.cae.
OVERWRITE_CHECKPOINT_FILES = True
SUBMIT_JOBS = False           # Set True to submit after all models are built.
WAIT_FOR_COMPLETION = True    # Only used when SUBMIT_JOBS is True.
REPLACE_EXISTING_MODELS = False
REPLACE_EXISTING_JOBS = True
# ========================================================================


# Abaqus cae noGUI does not define __file__.  The fallback supports the
# documented command executed from the project root.
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    script_argument = sys.argv[0] if len(sys.argv) else ''
    if script_argument.lower().endswith('.py'):
        SCRIPT_DIR = os.path.dirname(os.path.abspath(script_argument))
    else:
        SCRIPT_DIR = os.path.abspath(os.path.join(
            os.getcwd(), 'python script', 'Modeling'))
PROJECT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
if not os.path.isfile(os.path.join(SCRIPT_DIR, 'Matrix2D_Geom_Batch.py')):
    raise IOError("Cannot locate batch scripts in '{}'. Run the documented "
                  "command from the project root.".format(SCRIPT_DIR))
os.chdir(PROJECT_DIR)

RUN_CONFIG = {
    'MODEL_NAME_LIST': MODEL_NAME_LIST,
    'LABEL_LIST': LABEL_LIST,
    'SUFFIX': '_',
    'PART_NAME': 'Composite-cut',
    'RVE_SIZE': RVE_SIZE,
    'VOID_RADIUS': VOID_RADIUS,
    'MATERIAL_FILE': MATERIAL_FILE,
    'MODEL_DATA_DIR': MODEL_DATA_DIR,
    'MESH_SIZE': MESH_SIZE,
    'FORCE_Y': FORCE_Y,
    'JOB_CPUS': JOB_CPUS,
    'JOB_DOMAINS': JOB_DOMAINS,
    'WRITE_INPUT': WRITE_INPUT,
    'REPLACE_EXISTING_MODELS': REPLACE_EXISTING_MODELS,
    'REPLACE_EXISTING_JOBS': REPLACE_EXISTING_JOBS,
}


def target_names():
    return [base + RUN_CONFIG['SUFFIX'] + label for label in LABEL_LIST
            for base in MODEL_NAME_LIST]


def run_stage(file_name):
    path = os.path.join(SCRIPT_DIR, file_name)
    print('\n===== Running {} ====='.format(file_name))
    with open(path, 'r') as stage_file:
        code = compile(stage_file.read(), path, 'exec')
    exec(code, globals(), globals())


def require_single_base_model(feature_name):
    if len(MODEL_NAME_LIST) != 1:
        raise ValueError("{} requires exactly one item in MODEL_NAME_LIST. "
                         "Run each base model separately.".format(feature_name))


def checkpoint_path(stage):
    require_single_base_model('Checkpoint recovery')
    return os.path.join(PROJECT_DIR,
                        '{}_{}.cae'.format(MODEL_NAME_LIST[0], stage))


def save_database(path, overwrite, description):
    if os.path.exists(path):
        if not overwrite:
            raise IOError("{} already exists: {}".format(description, path))
        os.remove(path)
    mdb.saveAs(pathName=path)
    print('Saved {}: {}'.format(description, path))


def open_checkpoint(stage):
    global mdb
    path = checkpoint_path(stage)
    if not os.path.isfile(path):
        raise IOError("{} checkpoint not found: {}. Run with START_STAGE='all' "
                      "first.".format(stage, path))
    # Rebind the module-level handle so subsequently executed stage scripts
    # always operate on the database that was just opened.
    mdb = openMdb(pathName=path)
    print('Opened {} checkpoint: {}'.format(stage, path))


if START_STAGE == 'all':
    run_stage('Matrix2D_Geom_Batch.py')
    if SAVE_CHECKPOINTS:
        save_database(checkpoint_path('geometry'), OVERWRITE_CHECKPOINT_FILES,
                      'geometry checkpoint')
    run_stage('Mesh_Geom_Batch.py')
    if SAVE_CHECKPOINTS:
        save_database(checkpoint_path('mesh'), OVERWRITE_CHECKPOINT_FILES,
                      'mesh checkpoint')
elif START_STAGE == 'mesh':
    open_checkpoint('geometry')
    run_stage('Mesh_Geom_Batch.py')
    if SAVE_CHECKPOINTS:
        save_database(checkpoint_path('mesh'), OVERWRITE_CHECKPOINT_FILES,
                      'mesh checkpoint')
elif START_STAGE == 'load':
    open_checkpoint('mesh')
else:
    raise ValueError("START_STAGE must be 'all', 'mesh', or 'load'; got '{}'"
                     .format(START_STAGE))

run_stage('Load2D_static_Geom_Batch.py')

if SAVE_CAE:
    if CAE_FILE_NAME is None:
        if len(MODEL_NAME_LIST) != 1:
            raise ValueError("Set CAE_FILE_NAME explicitly when MODEL_NAME_LIST "
                             "contains more than one base model name.")
        cae_file_name = MODEL_NAME_LIST[0] + '.cae'
    else:
        cae_file_name = CAE_FILE_NAME
    cae_path = os.path.join(PROJECT_DIR, cae_file_name)
    save_database(cae_path, OVERWRITE_CAE_FILE, 'final CAE database')

if SUBMIT_JOBS:
    submitted_jobs = []
    for job_name in target_names():
        job = mdb.jobs[job_name]
        job.submit(consistencyChecking=OFF)
        submitted_jobs.append(job)
        print('Submitted {}'.format(job_name))
    if WAIT_FOR_COMPLETION:
        for job in submitted_jobs:
            job.waitForCompletion()
            print('Finished {}'.format(job.name))

print('\nWorkflow complete. Models/jobs: {}'.format(', '.join(target_names())))
