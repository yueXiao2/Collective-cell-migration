import sys
from os import environ
from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])

import CompuCellSetup

sim,simthread = CompuCellSetup.getCoreSimulationObjects()
        
# add extra attributes here      
CompuCellSetup.initializeSimulationObjects(sim,simthread)
# Definitions of additional Python-managed fields go here
        
#Add Python steppables here
steppableRegistry=CompuCellSetup.getSteppableRegistry()
from Epithelial1Steppables import ConstraintInitializerSteppable
CsteppableInstance=ConstraintInitializerSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(CsteppableInstance)

from Epithelial1Steppables import Epithelial1Steppable
steppableInstance=Epithelial1Steppable(sim,_frequency=1)
steppableRegistry.registerSteppable(steppableInstance)

from Epithelial1Steppables import GrowthSteppable
GsteppableInstance=GrowthSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(GsteppableInstance)

from Epithelial1Steppables import MitosisSteppable
MitosisSteppableInstance=MitosisSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(MitosisSteppableInstance)

#....................................................................................        
#Create extra player fields here or add attributes
#I added the following two lines from the manual file for Python scripting 3.6 page 21.  
dim=sim.getPotts().getCellFieldG().getDim()
VelocityVectorField = simthread.createVectorFieldCellLevelPy("V field")
VelocityScalarField=simthread.createScalarFieldCellLevelPy("V magnitude")
PolarityScalarField=simthread.createScalarFieldCellLevelPy("P magnitude")
PolarityAngleScalarField=simthread.createScalarFieldCellLevelPy("P angle")
ForceScalarField=simthread.createScalarFieldCellLevelPy("F magnitude")
ForceAngleScalarField=simthread.createScalarFieldCellLevelPy("F angle")
ForceXScalarField=simthread.createScalarFieldCellLevelPy("Fx")
ForceYScalarField=simthread.createScalarFieldCellLevelPy("Fy")

from Epithelial1Steppables import VectorFieldCellLevelVisualizationSteppable # This line inserter automatically (Right-click on "Epithelial1Steppables.py"; select Add steppables..).
instanceOfVectorFieldCellLevelVisualizationSteppable=VectorFieldCellLevelVisualizationSteppable(_simulator=sim,_frequency=1)# This line inserter automatically (Right-click on "Epithelial1Steppables.py"; select Add steppables..).
instanceOfVectorFieldCellLevelVisualizationSteppable.setVectorField(VelocityVectorField)#I added the following two lines from the manual file for Python scripting 3.6 page 21. 
steppableRegistry.registerSteppable(instanceOfVectorFieldCellLevelVisualizationSteppable)# This line inserter automatically (Right-click on "Epithelial1Steppables.py"; select Add steppables..).

from Epithelial1Steppables import ScalarFieldVisualizationSteppableV
instanceOfScalarFieldVisualizationSteppableV=ScalarFieldVisualizationSteppableV(_simulator=sim,_frequency=1)
instanceOfScalarFieldVisualizationSteppableV.setScalarField(VelocityScalarField)
steppableRegistry.registerSteppable(instanceOfScalarFieldVisualizationSteppableV)

from Epithelial1Steppables import ScalarFieldVisualizationSteppableP
instanceOfScalarFieldVisualizationSteppableP=ScalarFieldVisualizationSteppableP(_simulator=sim,_frequency=1)
instanceOfScalarFieldVisualizationSteppableP.setScalarField(PolarityScalarField)
steppableRegistry.registerSteppable(instanceOfScalarFieldVisualizationSteppableP)

from Epithelial1Steppables import ScalarFieldVisualizationSteppablePA
instanceOfScalarFieldVisualizationSteppablePA=ScalarFieldVisualizationSteppablePA(_simulator=sim,_frequency=1)
instanceOfScalarFieldVisualizationSteppablePA.setScalarField(PolarityAngleScalarField)
steppableRegistry.registerSteppable(instanceOfScalarFieldVisualizationSteppablePA)

from Epithelial1Steppables import ScalarFieldVisualizationSteppableF
instanceOfScalarFieldVisualizationSteppableF=ScalarFieldVisualizationSteppableF(_simulator=sim,_frequency=1)
instanceOfScalarFieldVisualizationSteppableF.setScalarField(ForceScalarField)
steppableRegistry.registerSteppable(instanceOfScalarFieldVisualizationSteppableF)

from Epithelial1Steppables import ScalarFieldVisualizationSteppableFA
instanceOfScalarFieldVisualizationSteppableFA=ScalarFieldVisualizationSteppableFA(_simulator=sim,_frequency=1)
instanceOfScalarFieldVisualizationSteppableFA.setScalarField(ForceAngleScalarField)
steppableRegistry.registerSteppable(instanceOfScalarFieldVisualizationSteppableFA)

from Epithelial1Steppables import ScalarFieldVisualizationSteppableFx
instanceOfScalarFieldVisualizationSteppableFx=ScalarFieldVisualizationSteppableFx(_simulator=sim,_frequency=1)
instanceOfScalarFieldVisualizationSteppableFx.setScalarField(ForceXScalarField)
steppableRegistry.registerSteppable(instanceOfScalarFieldVisualizationSteppableFx)

from Epithelial1Steppables import ScalarFieldVisualizationSteppableFy
instanceOfScalarFieldVisualizationSteppableFy=ScalarFieldVisualizationSteppableFy(_simulator=sim,_frequency=1)
instanceOfScalarFieldVisualizationSteppableFy.setScalarField(ForceYScalarField)
steppableRegistry.registerSteppable(instanceOfScalarFieldVisualizationSteppableFy)


#.................................................................................... 
        
CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        