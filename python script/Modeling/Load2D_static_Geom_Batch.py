"""Apply the reported static loading and create one Abaqus job per RVE."""
from abaqus import *
from abaqusConstants import *
from caeModules import *
import regionToolset

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
FORCE_Y = option('FORCE_Y', 1.0)
INITIAL_INCREMENT = option('INITIAL_INCREMENT', 1.0)
MAX_INCREMENT = option('MAX_INCREMENT', 1.0)
MIN_INCREMENT = option('MIN_INCREMENT', 1.0e-8)
MAX_NUM_INCREMENTS = option('MAX_NUM_INCREMENTS', 10000)
JOB_CPUS = option('JOB_CPUS', 4)
JOB_DOMAINS = option('JOB_DOMAINS', 4)
WRITE_INPUT = option('WRITE_INPUT', True)
REPLACE_EXISTING_JOBS = option('REPLACE_EXISTING_JOBS', True)


def target_names():
    return [base + SUFFIX + label for label in LABEL_LIST
            for base in MODEL_NAME_LIST]


def make_boundary_sets(model, l1, l2):
    """Create left/right and bottom/top edge sets using a scale-aware box."""
    part = model.parts[PART_NAME]
    tolerance = 0.0001 * min(l1, l2)
    boundaries = (
        ('Line0-1', part.edges.getByBoundingBox(-tolerance, 0.0, 0.0,
                                                 tolerance, l2, 0.0)),
        ('Line1-1', part.edges.getByBoundingBox(l1 - tolerance, 0.0, 0.0,
                                                 l1 + tolerance, l2, 0.0)),
        ('Line0-2', part.edges.getByBoundingBox(0.0, -tolerance, 0.0,
                                                 l1, tolerance, 0.0)),
        ('Line1-2', part.edges.getByBoundingBox(0.0, l2 - tolerance, 0.0,
                                                 l1, l2 + tolerance, 0.0)),
    )
    for name, edges in boundaries:
        if len(edges) == 0:
            raise ValueError("Could not find boundary '{}' in model '{}'"
                             .format(name, model.name))
        part.Set(edges=edges, name=name)


for model_name in target_names():
    model = mdb.models[model_name]
    assembly = model.rootAssembly
    rp0 = assembly.sets['Set-rp0']
    rp1 = assembly.sets['Set-rp1']
    # Do not depend on RP-2: its name/id can vary after CAE edits.
    l1, l2 = RVE_SIZE

    model.StaticStep(name='Step-1', previous='Initial', nlgeom=OFF,
                     maxNumInc=MAX_NUM_INCREMENTS,
                     initialInc=INITIAL_INCREMENT, minInc=MIN_INCREMENT,
                     maxInc=MAX_INCREMENT)
    model.FieldOutputRequest(name='F-Output-1', createStepName='Step-1',
                             variables=('LE', 'PE', 'S', 'TRIAX', 'U', 'EVOL'),
                             numIntervals=1, exteriorOnly=OFF)
    model.HistoryOutputRequest(name='H-Output-energy', createStepName='Step-1',
                               variables=('ALLAE', 'ALLCD', 'ALLDMD', 'ALLFD',
                                          'ALLIE', 'ALLPD', 'ALLSE', 'ALLVD',
                                          'ALLWK', 'ETOTAL'), numIntervals=1)
    model.HistoryOutputRequest(name='H-Output-rp0', createStepName='Step-1',
                               variables=('U1', 'U2', 'U3', 'UR1', 'UR2', 'UR3',
                                          'RF1', 'RF2', 'RF3', 'RM1', 'RM2', 'RM3'),
                               numIntervals=1, region=rp0, sectionPoints=DEFAULT,
                               rebar=EXCLUDE)
    model.HistoryOutputRequest(name='H-Output-rp1', createStepName='Step-1',
                               variables=('U1', 'U2', 'U3', 'UR1', 'UR2', 'UR3',
                                          'RF1', 'RF2', 'RF3', 'RM1', 'RM2', 'RM3'),
                               numIntervals=1, region=rp1, sectionPoints=DEFAULT,
                               rebar=EXCLUDE)

    make_boundary_sets(model, l1, l2)
    # Opposite boundary displacement components are tied to reference points.
    # This reproduces the constraint convention used in the reported models.
    for dof in range(2):
        component = dof + 1
        model.Equation(name='EquationConstraint0-{}'.format(component),
                       terms=((1.0, PART_NAME + '.Line0-{}'.format(component), component),
                              (-1.0, 'Set-rp0', component)))
        model.Equation(name='EquationConstraint1-{}'.format(component),
                       terms=((1.0, PART_NAME + '.Line1-{}'.format(component), component),
                              (-1.0, 'Set-rp1', component)))

    model.ConcentratedForce(name='Load-1', createStepName='Step-1', region=rp1,
                            cf2=FORCE_Y * l1, distributionType=UNIFORM,
                            field='', localCsys=None)
    model.DisplacementBC(name='Load-0', createStepName='Step-1', region=rp0,
                         u1=0.0, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET,
                         ur3=UNSET, amplitude=UNSET, fixed=OFF,
                         distributionType=UNIFORM, fieldName='', localCsys=None)

    if model_name in mdb.jobs.keys():
        if not REPLACE_EXISTING_JOBS:
            raise ValueError("Job '{}' already exists.".format(model_name))
        del mdb.jobs[model_name]
    mdb.Job(name=model_name, model=model_name, userSubroutine='',
            numCpus=JOB_CPUS, numDomains=JOB_DOMAINS,
            numThreadsPerMpiProcess=1, activateLoadBalancing=False)
    if WRITE_INPUT:
        mdb.jobs[model_name].writeInput(consistencyChecking=OFF)
    print('Load and job created for {}'.format(model_name))
