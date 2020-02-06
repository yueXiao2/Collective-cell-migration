from PySteppables import *
import CompuCell
import sys
#import matplotlib.pyplot as plt
import random # Mine
#from numpy import *
#from math import *
import numpy
import numpy as np
import math
import os.path
from XMLUtils import dictionaryToMapStrStr as d2mss
import random as ra
# 17-Sep-2019
##########################################################################################################
#     v = (dCOM_x/dt, dCOM_y/dt) 
#     dp/dt = -Beta p + Gamma V    (in x and y directions)
#    
#     COM = COM + v dt             (in x and y directions)               
#     p = p + dp/dt dt
#   
#     F = Fmax (p/|p|)(|p|^n/[|p|^n + Alpha^n])   (in x and y directions, i.e. lambdaVecX and lambdaVecY)
#
#     Force F (applied to a cell) is calculated for all cells. Then, F is pluged into energy function below.
#
#                           Area elasticity       Contractility(actin)       Adhesion(E-cadherin)   
#     Potts energy = SUM_AllCells_i[Ai - At]^2 + SUM_AllCells_i[Li - Lt]^2 - SUM_AllCells_i[Contact_i]   
#
##########################################################################################################
#------------- Global parameters for scanning -------------------
n_Hill = 10
Alpha = 1
Beta = 0.1
Gamma = 1
F_max = -1000
Polarization_threshold_UnstableSteadyPolarity = 0.9   
Initial_Polarization = 4 # This should be > unstable steady polarity (see single-cell simulations).
initial_position_edge_x = 0
Noise_Angle = 0 # Noise magnitude for direction of polarity (in degrees). 45 degrees generated some backward motion. 
MaintainPolarization_LeadingCells = 1 # 1: if leading cells maintain the Initial_Polarization for the entire period of simulation.
val_Dimension_y = 0
cell_size = 10 # Eeverycell is 10*10 pixels.
#EdgeBoxPolarity = 6 #cell.dict['EdgeBox']
Edge_Box_width_x = 3

val_LambdaVolume = 70
MCS_Closed_Edge_cells_Zero_F = 2900/2    # Number of MCS with: barrier is on and zero force (cell polarity is not linked to F).       50; 500         
MCS_Closed_Edge_cells_NonZero_F = 600/2 # Number of MCS with: barrier is on and cell can apply force (cell polarity is linked to F). 20; 500 
mcs_T1 = MCS_Closed_Edge_cells_Zero_F + MCS_Closed_Edge_cells_NonZero_F + 10  # T1: SUM+10
mcs_T2 = 3710/2  # T2: SUM+100
mcs_T3 = 3960/2  # T3: SUM+300
mcs_T4 = 4160/2 # T4: SUM+500
#---------------------------------------------------------------
class Epithelial1Steppable(SteppableBasePy):
    
    #...........................
    #...........................
    #...........................
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    
    #...........................
    #...........................
    #...........................
    def start(self):
        # any code in the start function runs before MCS=0        
        for cell in self.cellList:
            #if cell.type != 3:
            if cell.type == 1 or cell.type == 2:
                cell.dict['n'] = n_Hill 
                cell.dict['alpha'] = Alpha 
                cell.dict['beta'] = Beta 
                cell.dict['gamma'] = Gamma  
                cell.dict['Fmax'] = F_max
                cell.dict['P_threshold'] = Polarization_threshold_UnstableSteadyPolarity # This is the unstable steady polarity.  
                
                cell.dict['PreviousComX'] = cell.xCOM     
                cell.dict['PreviousComY'] = cell.yCOM 
                
                cell.dict['velocityX'] = 0  
                cell.dict['velocityY'] = 0
                cell.dict['VelocityMagnitude'] = 0
                cell.dict['VelocityAngle'] = 0               
                 
                cell.dict['PolarityX'] = 0   
                cell.dict['PolarityY'] = 0 
                cell.dict['PolarityMagnitude'] = 0
                cell.dict['PolarityAngle'] = 0 
                
                cell.lambdaVecX = 0
                cell.lambdaVecY = 0                 
                cell.dict['ForceMagnitude'] = 0
                cell.dict['ForceAngle'] = 0
                
                cell.dict['EdgeRegion'] = 0 # If the cell is in the edge region of the sheet.
                cell.dict['EndRegion'] = 0  # If the cell is in the end region of the sheet.
                cell.dict['EdgeBox'] = 0  # If the cell is in the most front box for which density > (0.5 * density of a non-moving box).
                
                cell.dict['initial_position_edge_x'] = 0  # For calculating the displacement of edge cells (1 layer) form the initial position of the edge.
                
                cell.dict['IsConnectedToSheet'] = 0 
                cell.dict['IsConnectedToMedium'] = 0 # If the cell is close to the free area.
                
                
                cell.dict['comXbeforePro'] = 0
                cell.dict['comYbeforePro'] = 0

        #=============================== Print in File ==========================================      
        #.........Read XML parameters......
        
        pottsXMLData=self.simulator.getCC3DModuleData("Potts")
        if pottsXMLData:
            temperatureElement=pottsXMLData.getFirstElement("Temperature")
            if temperatureElement:
                val_Temperature = (temperatureElement.getText())  
            Dimension_y=pottsXMLData.getFirstElement("Dimensions",d2mss({"x":"1500"}))
            if Dimension_y:
                val_Dimension_y = Dimension_y.getAttributeAsDouble("y")
           
        LambdaVolumeXMLData=self.simulator.getCC3DModuleData("Plugin","Volume") 
        if LambdaVolumeXMLData:   
            Volume_Pol = LambdaVolumeXMLData.getFirstElement("VolumeEnergyParameters",d2mss({"CellType":"Polarized"})) 
            if Volume_Pol:
                val_LambdaVolume=Volume_Pol.getAttributeAsDouble("LambdaVolume")
        
        contactXMLData=self.simulator.getCC3DModuleData("Plugin","Contact") # get <Plugin Name="Contact"> section of XML file
        if contactXMLData:   # check if we were able to successfully get the section from simulator
            ContactEnergy_PolPol = contactXMLData.getFirstElement("Energy",d2mss({"Type1":"Polarized","Type2":"Polarized"})) # get <Energy Type1="Polarized" Type2="Unpolarized"> element            
            if ContactEnergy_PolPol: # check if the attempt was succesful
                val_ContactEnergy_PolPol=(ContactEnergy_PolPol.getText()) # get value of <Energy Type1="Polarized" Type2="Polarized"> element and convert it into float
        
        FPP_XMLData=self.simulator.getCC3DModuleData("Plugin","FocalPointPlasticity") 
        if FPP_XMLData:   
            FPP_PolPol = FPP_XMLData.getFirstElement("Parameters",d2mss({"Type1":"Polarized","Type2":"Polarized"}))          
            if FPP_PolPol: 
                FPP_Lambda = FPP_PolPol.getFirstElement("Lambda")
                if FPP_Lambda:
                   val_FPP_Lambda = (FPP_Lambda.getText())
                
        LambdaSurfaceXMLData=self.simulator.getCC3DModuleData("Plugin","Surface") 
        if LambdaSurfaceXMLData:   
            SurfaceEnergy_Pol = LambdaSurfaceXMLData.getFirstElement("SurfaceEnergyParameters",d2mss({"CellType":"Polarized","TargetSurface":"0"})) 
            if SurfaceEnergy_Pol:
                val_LambdaSurface=SurfaceEnergy_Pol.getAttributeAsDouble("LambdaSurface")       
        #..................................
        File_path = 'C:\Users\qq871\Desktop\Multi'        
        if Alpha == 1: 
           File_path = 'C:\Users\qq871\Desktop\Multi\ParameterScan\Alpha1'            
        if Alpha == 0.5: 
           File_path = 'C:\Users\qq871\Desktop\Multi\ParameterScan\Alpha05' 
        FileGeneral_path = os.path.join(File_path, "Epithelial_Multi.txt") # Path to the general file "Epithelial_Multi.txt".
        txt_f = open(FileGeneral_path, "a") # Opens the file "Epithelial_Multi.txt" to write in.
        
        Cells_str = ''
        Cells_str += "---------------------------------------------------\n"
        Cells_str += "The initial properties of all cells are: \n" 
        txt_f.write(Cells_str) # Writes the content of "Cells_str" in the file "Epithelial_Multi.txt".
        for cell in self.cellList:
            if cell.type == 1 or cell.type == 2:
                Cells_str = ''
                Cells_str += "Cell ID " + str(cell.id) + "; " + "COM = (" + str(round(cell.xCOM, 3)) + "," + str(round(cell.yCOM, 3)) + "); " + "V = (" + str(round(cell.dict['velocityX'], 3)) + "," + str(round(cell.dict['velocityY'], 3)) + "); " + "P = (" + str(round(cell.dict['PolarityX'], 3)) + "," + str(round(cell.dict['PolarityY'], 3)) + "); " + "|P| = " + str(round(cell.dict['PolarityMagnitude'], 3)) + "; " + "F = (" + str(round(cell.lambdaVecX, 3)) + "," + str(round(cell.lambdaVecY, 3)) + "); " +  "\n"
                txt_f.write(Cells_str) # Writes the content of "Cells_str" in the file "Epithelial_Multi.txt".
        Cells_str = "---------------------------------------------------\n\n"
        txt_f.write(Cells_str) 
        
        FileMCSF_path = os.path.join(File_path, "MatlabMCSF.txt")
        txt_Matlab_MCS_F = open(FileMCSF_path, "a")
        txt_Matlab_MCS_F.write(str(MCS_Closed_Edge_cells_Zero_F) + " " + str(MCS_Closed_Edge_cells_NonZero_F) + "\n")
        
        File_Pvalues_path = os.path.join(File_path, "MatlabPvalues.txt")
        txt_Matlab_Pvalues = open(File_Pvalues_path, "a")  
        val_LambdaVolume = 70
        #txt_Matlab_Pvalues.write("Fmax = " + str(F_max) + "\n" + "Temperature = " + str(val_Temperature) + "\n" + "Lambda_Volume = " + str(val_LambdaVolume) + "\n" + "Contact PoLPol = " + str(val_ContactEnergy_PolPol) +  "\n" + "LambdaSurface = " + str(val_LambdaSurface) + "\n" + "Lambda_Spring(FPP) = " + str(val_FPP_Lambda) + "\n")
        txt_Matlab_Pvalues.write("Edge_Box_width_x = " + str(Edge_Box_width_x) + "\n" + "Dimension_y = " + str(val_Dimension_y) + " pixels \n" + "Maintain polarization for leading cells (1/0): " + str(MaintainPolarization_LeadingCells) + "\n" + "Beta = " + str(Beta)+ "\n" +"Alpha = " + str(Alpha)+ "\n" + "Fmax = " + str(F_max) + "\n" + "Noise_Angle = " + str(Noise_Angle) + "\n" +"Temperature = " + str(val_Temperature) + "\n" + "Lambda_Volume = " + str(val_LambdaVolume) + "\n" + "Contact PoLPol = " + str(val_ContactEnergy_PolPol) +  "\n" + "LambdaSurface = " + str(val_LambdaSurface) + "\n" + "Lambda_Spring(FPP) = No Spring")
        
        txt_f.close() # Closes file "Epithelial_Multi.txt", when I finished the writing (it must be closed).   
        txt_Matlab_MCS_F.close() 
        txt_Matlab_Pvalues.close()
        #======================================================================================== 
        pass
        
    #...........................
    #...........................
    #...........................    
    def step(self,mcs):          
        #************************************************************************* 
        #***************************XML Info************************************** 
        #************************************************************************* 
        pottsXMLData=self.simulator.getCC3DModuleData("Potts")
        if pottsXMLData:
            Dimension_y=pottsXMLData.getFirstElement("Dimensions",d2mss({"x":"1500"}))
            if Dimension_y:
                val_Dimension_y = Dimension_y.getAttributeAsDouble("y")
        #************************************************************************* 
        #************************************************************************* 
        #*************************************************************************  
        MCS_Barrier_removal = (MCS_Closed_Edge_cells_Zero_F + MCS_Closed_Edge_cells_NonZero_F) # Number of MCS with: barrier is removed and cell can apply force (cell polarity is linked to F).        
        MCS_polarizing_edge = MCS_Barrier_removal + 1        
        #======================================================================================== 
        #=============================Remove the barrier only==================================== 
        #======================================================================================== 
        if mcs == MCS_Barrier_removal:
            for cell in self.cellList:
                if cell.type == 1 or cell.type ==  2:
                    cell.dict['initVolume'] = cell.volume
                if cell.type == self.BARRIER:
                    self.deleteCell(cell)# To be able to delete a cell, I must put PixelTracker Plugin in the CC3DML.
        #======================================================================================== 
        #==============After barrier removal: find Cells at the Edge and polarize them=========== 
        #======================================================================================== 
        if mcs == MCS_polarizing_edge: 
            #-------- Estimate an edge region ---------
            COM_X_max = 0 
            for cell in self.cellList:
                if cell.type == 1 or cell.type == 2: 
                    if cell.xCOM >= COM_X_max:
                        COM_X_max = cell.xCOM
            
            Edge_region = COM_X_max - (10*3) # Edge_region = COM_X_max - (3 * cell_size). Edge region: 3 layers from the edge..              
            for cell in self.cellList:
                if cell.type == 1 or cell.type == 2: 
                    if cell.xCOM >= Edge_region:
                        cell.dict['EdgeRegion'] = 1
            #------------------------------------------
            #......Find cells at the edge (first layer only) and polarise them ..........
#             for cell in self.cellList:
#                 if cell.type == 1 or cell.type == 2: 
#                     if cell.dict['EdgeRegion'] == 1:
#                         for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell): 
#                             if neighbor == None: # Only the first layer is plarised.
#                                 cell.type = 1 
            
            for cell in self.cellList: 
                if cell.type == 1 or cell.type == 2: 
                    if cell.dict['EdgeRegion'] == 1: # The first three layers are plarised.
                        cell.type = 1             
            #............................................................................
            #******** Initialise the polarised cells *********
            for cell in self.cellList:
                cell.dict['PreviousComX'] = cell.xCOM     
                cell.dict['PreviousComY'] = cell.yCOM  
                if cell.type == 1: 
                    cell.dict['PolarityMagnitude'] = Initial_Polarization                   
#                     cell.dict['PolarityAngle'] = ra.uniform(-1,1) * 60 # This is a random angle in degrees.
#                     cell.dict['PolarityX'] = cell.dict['PolarityMagnitude'] * math.cos(math.radians(cell.dict['PolarityAngle']))
#                     cell.dict['PolarityY'] = cell.dict['PolarityMagnitude'] * math.sin(math.radians(cell.dict['PolarityAngle']))                    
                    cell.dict['PolarityX'] = Initial_Polarization
                    cell.dict['PolarityY'] = 0
                    cell.dict['PolarityAngle'] = 0
                    
                    cell.lambdaVecX = cell.dict['Fmax'] * (cell.dict['PolarityX'] / cell.dict['PolarityMagnitude']) * ((cell.dict['PolarityMagnitude']**cell.dict['n'])/((cell.dict['PolarityMagnitude']**cell.dict['n']) + (cell.dict['alpha']**cell.dict['n'])))# Fx = (Px/|Px|)(|Px|^n/[|Px|^n + alpha^n])
                    cell.lambdaVecY = cell.dict['Fmax'] * (cell.dict['PolarityY'] / cell.dict['PolarityMagnitude']) * ((cell.dict['PolarityMagnitude']**cell.dict['n'])/((cell.dict['PolarityMagnitude']**cell.dict['n']) + (cell.dict['alpha']**cell.dict['n'])))
                    cell.dict['ForceMagnitude'] = math.sqrt(cell.lambdaVecX**2 + cell.lambdaVecY**2)
                    if cell.dict['ForceMagnitude'] == 0:
                        cell.dict['ForceAngle'] = 0
                    else:
                        if cell.lambdaVecY < 0:
                            cell.dict['ForceAngle'] = math.degrees(math.acos(((-1) * cell.lambdaVecX) / cell.dict['ForceMagnitude'])) 
                        if cell.lambdaVecY > 0:
                            cell.dict['ForceAngle'] = (-1) * math.degrees(math.acos(((-1) * cell.lambdaVecX) / cell.dict['ForceMagnitude'])) 
                        if cell.lambdaVecY == 0:
                            if cell.lambdaVecX < 0:
                                cell.dict['ForceAngle'] = 0
                            if cell.lambdaVecX > 0:
                                cell.dict['ForceAngle'] = 180         
#             for cell in self.cellList:
#                 cell.dict['PreviousComX'] = cell.xCOM     
#                 cell.dict['PreviousComY'] = cell.yCOM                 
#                 if cell.type == 1:                    
#                     cell.dict['PolarityX'] = cell.dict['P_threshold']   
#                     cell.dict['PolarityY'] = 0
#                     cell.dict['PolarityMagnitude'] = math.sqrt(cell.dict['PolarityX']**2 + cell.dict['PolarityY']**2)
#                     cell.dict['PolarityAngle'] = 0 
#                     cell.lambdaVecX = cell.dict['Fmax'] * (cell.dict['PolarityX'] / cell.dict['PolarityMagnitude']) * ((cell.dict['PolarityMagnitude']**cell.dict['n'])/((cell.dict['PolarityMagnitude']**cell.dict['n']) + (cell.dict['alpha']**cell.dict['n'])))# Fx = (Px/|Px|)(|Px|^n/[|Px|^n + alpha^n])
#                     cell.lambdaVecY = 0
#                     cell.dict['ForceMagnitude'] = math.sqrt(cell.lambdaVecX**2 + cell.lambdaVecY**2)
#                     cell.dict['ForceAngle'] = 0 
#                 if cell.type == 2:
#                     cell.dict['PolarityX'] = 0   
#                     cell.dict['PolarityY'] = 0 
#                     cell.dict['PolarityMagnitude'] = 0
#                     cell.dict['PolarityAngle'] = 0  
#                     cell.lambdaVecX = 0
#                     cell.lambdaVecY = 0 
#                     cell.dict['ForceMagnitude'] = 0
#                     cell.dict['ForceAngle'] = 0  
            #****************************************************      
        #========================================================================================= 
        #====After barrier removal: Stop migration if polarization has reached the end boundary===
        #====Before barrier removal: cells do nothing (just jiggle in their position)=============       
        #========================================================================================= 
        #.......If polarity has reached the end region?......
        Polarity_reached_end = 0
        if mcs > MCS_polarizing_edge:            
            End_region = (10 * 1) 
            for cell in self.cellList:
                if cell.type == 1 or cell.type == 2:
                    if cell.dict['EndRegion'] == 1:
                        Polarity_reached_end = 1 
               
            if Polarity_reached_end == 0:
               for cell in self.cellList:
                   if cell.type == 1 or cell.type == 2:
                       if cell.xCOM <= End_region:
                           if cell.type == 1:
                               cell.dict['EndRegion'] = 1
                               Polarity_reached_end = 1  
        #....................................................
#         #*******Stop migration (after barrier removal)********
        if ((mcs > MCS_polarizing_edge) and (Polarity_reached_end == 1)): #Polarization has reached the end of lattice.
            P_end = 1
            for cell in self.cellList:
                cell.dict['PreviousComX'] = cell.xCOM     
                cell.dict['PreviousComY'] = cell.yCOM 
                if cell.type == 1 or cell.type == 2:                   
                    cell.dict['PolarityX'] = 0   
                    cell.dict['PolarityY'] = 0
                    cell.dict['PolarityMagnitude'] = 0
                    cell.dict['PolarityAngle'] = 0 
                    cell.lambdaVecX = 0
                    cell.lambdaVecY = 0
                    cell.dict['ForceMagnitude'] = 0
                    cell.dict['ForceAngle'] = 0                     
        #*******Stop migration (before barrier removal)********
        initial_position_edge_x_temp = 0 # Saves the initial position of the edge (in x axis) attached to the barrier.
        if (mcs < MCS_Closed_Edge_cells_Zero_F): # Before removing barrier when cells have no force.
            for cell in self.cellList:
                cell.dict['PreviousComX'] = cell.xCOM     
                cell.dict['PreviousComY'] = cell.yCOM 
                if cell.type == 1 or cell.type == 2:                   
                    cell.dict['PolarityX'] = 0   
                    cell.dict['PolarityY'] = 0
                    cell.dict['PolarityMagnitude'] = 0
                    cell.dict['PolarityAngle'] = 0 
                    cell.lambdaVecX = 0
                    cell.lambdaVecY = 0
                    cell.dict['ForceMagnitude'] = 0
                    cell.dict['ForceAngle'] = 0 
                    if cell.xCOM >= initial_position_edge_x_temp:
                        initial_position_edge_x_temp = cell.xCOM
            for cell in self.cellList: 
                if cell.type == 1 or cell.type == 2:
                   cell.dict['initial_position_edge_x'] = initial_position_edge_x_temp
        min_P_angle = 180
        max_P_angle = -180
        #*******************************************************                
        #========================================================================================      
        #==========================After barrier removal: cells migrate==========================
        #========================================================================================
        #if ((mcs > MCS_Closed_Edge_cells_Zero_F) and (Polarity_reached_end == 0) and (mcs != MCS_polarizing_edge)):
        cells_generate_force = 0
        cells_generate_force_barrier_On = 0
        if (MCS_Closed_Edge_cells_Zero_F < mcs < MCS_Barrier_removal):# Cells apply force while barrier is on.
            cells_generate_force = 1
            cells_generate_force_barrier_On = 1
        if ((mcs > MCS_polarizing_edge) and (Polarity_reached_end == 0)): # Cells apply force after barrier removal and while polarization has not reached the end.
            cells_generate_force = 1                         
            
        if (cells_generate_force == 1):                         
            File_path = 'C:\Users\qq871\Desktop\Multi'
            if Alpha == 1:
                File_path = 'C:\Users\qq871\Desktop\Multi\ParameterScan\Alpha1'
            if Alpha == 0.5: 
                File_path = 'C:\Users\qq871\Desktop\Multi\ParameterScan\Alpha05'    
            Cells_str = ''
            Cells_str += "At MCS = " + str(mcs)+"\n"
            FileGeneral_path = os.path.join(File_path, "Epithelial_Multi.txt")  
            txt_f = open(FileGeneral_path, "a") 
            txt_f.write(Cells_str) # Writes the content of "Cells_str" in the file "Epithelial_Multi.txt".

            File_NCOMX_path = os.path.join(File_path, "MatlabNCOMX.txt")
            txt_Matlab_NCOMx = open(File_NCOMX_path, "a")  

            FileCOMX_path = os.path.join(File_path, "MatlabCOMX.txt")
            txt_Matlab_COMx = open(FileCOMX_path, "a")  
            FileCOMY_path = os.path.join(File_path, "MatlabCOMY.txt")
            txt_Matlab_COMy = open(FileCOMY_path, "a")
            
            FileVelocity_path = os.path.join(File_path, "MatlabVelocity.txt")
            txt_Matlab_V = open(FileVelocity_path, "a")                 
            FileVelocityAngle_path = os.path.join(File_path, "MatlabVelocityAngle.txt")
            txt_Matlab_VA = open(FileVelocityAngle_path, "a")            
            FileOrderV_PCells_path = os.path.join(File_path, "MatlabOrderVPCells.txt")
            txt_Matlab_OV_PCells = open(FileOrderV_PCells_path, "a")           
            FileOrderV_AllCells_path = os.path.join(File_path, "MatlabOrderVAllCells.txt")
            txt_Matlab_OV_AllCells = open(FileOrderV_AllCells_path, "a")  
            
            FilePolarity_path = os.path.join(File_path, "MatlabPolarity.txt")
            txt_Matlab_P = open(FilePolarity_path,"a")             
            FilePolarityAngle_path = os.path.join(File_path, "MatlabPolarityAngle.txt")
            txt_Matlab_PA = open(FilePolarityAngle_path, "a")            
            FilePolarityAngleEdge_path = os.path.join(File_path, "MatlabPolarityAngleEdge.txt")
            txt_Matlab_PA_Edge = open(FilePolarityAngleEdge_path, "a")            
            FileOrderP_PCells_path = os.path.join(File_path, "MatlabOrderPPCells.txt")
            txt_Matlab_OP_PCells = open(FileOrderP_PCells_path, "a")           
            FileOrderP_AllCells_path = os.path.join(File_path, "MatlabOrderPAllCells.txt")
            txt_Matlab_OP_AllCells = open(FileOrderP_AllCells_path, "a")  
            FilePolarityAngleMinMax_path = os.path.join(File_path, "MatlabPolarityAngleMinMax.txt")
            txt_Matlab_PA_MinMax = open(FilePolarityAngleMinMax_path, "a")           
            
            FileForce_path = os.path.join(File_path, "MatlabForce.txt")
            txt_Matlab_F = open(FileForce_path, "a") 
            FileForceX_path = os.path.join(File_path, "MatlabForceX.txt")
            txt_Matlab_FX = open(FileForceX_path, "a") 
            FileForceY_path = os.path.join(File_path, "MatlabForceY.txt")
            txt_Matlab_FY = open(FileForceY_path, "a")             
            FileForceAngle_path = os.path.join(File_path, "MatlabForceAngle.txt")
            txt_Matlab_FA = open(FileForceAngle_path, "a")
            FileOrderF_PCells_path = os.path.join(File_path, "MatlabOrderFPCells.txt")
            txt_Matlab_OF_PCells = open(FileOrderF_PCells_path, "a")           
            FileOrderF_AllCells_path = os.path.join(File_path, "MatlabOrderFAllCells.txt")
            txt_Matlab_OF_AllCells = open(FileOrderF_AllCells_path, "a")

            FileNp_path = os.path.join(File_path, "MatlabNp.txt")
            txt_Matlab_Np = open(FileNp_path, "a")

            FileVEdgeMean_path = os.path.join(File_path, "MatlabVEdgeMean.txt")
            txt_Matlab_VEdge_Mean = open(FileVEdgeMean_path, "a") 
            FileVEdgeSTD_path = os.path.join(File_path, "MatlabVEdgeSTD.txt")
            txt_Matlab_VEdge_STD = open(FileVEdgeSTD_path, "a")            
        
            FileDisplEdgeMean_path = os.path.join(File_path, "MatlabDisplEdgeMean.txt")
            txt_Matlab_DisplEdge_Mean = open(FileDisplEdgeMean_path, "a") 
            FileDisplEdgeSTD_path = os.path.join(File_path, "MatlabDisplEdgeSTD.txt")
            txt_Matlab_DisplEdge_STD = open(FileDisplEdgeSTD_path, "a") 
            
            FileAreaMeanT1_path = os.path.join(File_path, "MatlabAreaMeanT1.txt")
            txt_Matlab_Area_Mean_T1 = open(FileAreaMeanT1_path, "a") 
            FileAreaSTDT1_path = os.path.join(File_path, "MatlabAreaSTDT1.txt")
            txt_Matlab_Area_STD_T1 = open(FileAreaSTDT1_path, "a")
            
            
            FileSurfaceMeanT1_path = os.path.join(File_path, "MatlabSurfaceMeanT1.txt")
            txt_Matlab_Surface_Mean_T1 = open(FileSurfaceMeanT1_path, "a") 
            FileSurfaceSTDT1_path = os.path.join(File_path, "MatlabSurfaceSTDT1.txt")
            txt_Matlab_Surface_STD_T1 = open(FileSurfaceSTDT1_path, "a") 
           
       
            FileEccentricityMeanT1_path = os.path.join(File_path, "MatlabEccentricityMeanT1.txt")
            txt_Matlab_Eccentricity_Mean_T1 = open(FileEccentricityMeanT1_path, "a") 
            FileEccentricitySTDT1_path = os.path.join(File_path, "MatlabEccentricitySTDT1.txt")
            txt_Matlab_Eccentricity_STD_T1 = open(FileEccentricitySTDT1_path, "a") 
            
      
            FileDensityT1_path = os.path.join(File_path, "MatlabDensityT1.txt")
            txt_Matlab_Density_T1 = open(FileDensityT1_path, "a") 
            

            FileVelocityMeanT1_path = os.path.join(File_path, "MatlabVelocityMeanT1.txt")
            txt_Matlab_Velocity_Mean_T1 = open(FileVelocityMeanT1_path, "a") 
            FileVelocitySTDT1_path = os.path.join(File_path, "MatlabVelocitySTDT1.txt")
            txt_Matlab_Velocity_STD_T1 = open(FileVelocitySTDT1_path, "a") 
            
                        
            FileVelocityXMeanT1_path = os.path.join(File_path, "MatlabVelocityXMeanT1.txt")
            txt_Matlab_VelocityX_Mean_T1 = open(FileVelocityXMeanT1_path, "a") 
            FileVelocityXSTDT1_path = os.path.join(File_path, "MatlabVelocityXSTDT1.txt")
            txt_Matlab_VelocityX_STD_T1 = open(FileVelocityXSTDT1_path, "a") 
           

            FilePolarityMeanT1_path = os.path.join(File_path, "MatlabPolarityMeanT1.txt")
            txt_Matlab_Polarity_Mean_T1 = open(FilePolarityMeanT1_path, "a") 
            FilePolaritySTDT1_path = os.path.join(File_path, "MatlabPolaritySTDT1.txt")
            txt_Matlab_Polarity_STD_T1 = open(FilePolaritySTDT1_path, "a") 
            

            FileForceMeanT1_path = os.path.join(File_path, "MatlabForceMeanT1.txt")
            txt_Matlab_Force_Mean_T1 = open(FileForceMeanT1_path, "a") 
            FileForceSTDT1_path = os.path.join(File_path, "MatlabForceSTDT1.txt")
            txt_Matlab_Force_STD_T1 = open(FileForceSTDT1_path, "a") 
            
            FilebinX_path = os.path.join(File_path, "MatlabbinX.txt")
            txt_Matlab_bin = open(FilebinX_path, "a") 

        
            FileMaxCOMX_path = os.path.join(File_path, "MatlabMaxCOMX.txt")
            txt_Matlab_MAXCOM_x = open(FileMaxCOMX_path, "a")
             
            
            txt_Matlab_COMx.write(str(mcs) + " ")
            txt_Matlab_COMy.write(str(mcs) + " ")
            txt_Matlab_V.write(str(mcs) + " ")
#            txt_Matlab_VA.write(str(mcs) + " ")
            txt_Matlab_P.write(str(mcs) + " ")
#             txt_Matlab_PA.write("MCS = " + str(mcs) + "\n " + "[Cell ID: P_angle] = ")
            txt_Matlab_PA_MinMax.write(str(mcs) + " ")
#             if mcs > (MCS_polarizing_edge + 1):
#                 txt_Matlab_PA_Edge.write(str(mcs) + " ")
            txt_Matlab_F.write(str(mcs) + " ")
#            txt_Matlab_FA.write("MCS = " + str(mcs) + "\n " + "[Cell ID: F_angle] = ")            
            txt_Matlab_FX.write("MCS = " + str(mcs) + "\n " + "[Cell ID: FX] = ")  
            txt_Matlab_FY.write("MCS = " + str(mcs) + "\n " + "[Cell ID: FY] = ")       
            txt_Matlab_NCOMx.write(str(mcs) + " ")
            #-------- Number of polarized cells ---------
            n_polarised_cells = 0
            for cell in self.cellList:
                if cell.type == 1:
                    n_polarised_cells = n_polarised_cells + 1 
            txt_Matlab_Np.write(str(mcs) + " " + str(n_polarised_cells) + "\n")
            #--------------------------------------------
            #.......... Average velocity and displacement of edge ........            
            VelocitiesEdge = [ ]
            Mean_Velocity_Edge = 0
            std_Velocity_Edge = 0           
            DisplsEdge = [ ]
            Mean_Displ_Edge = 0
            std_Displ_Edge = 0   
            if cells_generate_force_barrier_On == 1: # While barrier is On.
                COM_X_max = 0# Max COM at every MCS.
                for cell in self.cellList:
                    if cell.type == 1 or cell.type == 2: 
                        if cell.xCOM >= COM_X_max: # unpolarised cell is on the edge layer, attached to the barrier.
                            COM_X_max = cell.xCOM

                for cell in self.cellList:
                   if cell.type == 1 or cell.type == 2:
                       if cell.xCOM >= (COM_X_max - 3): # Make sure to consider all unpolarised cells attached to the barrier. 
                            VelocitiesEdge.append(cell.dict['velocityX'])
                            DisplsEdge.append((cell.xCOM - cell.dict['initial_position_edge_x']))   
            
            if cells_generate_force_barrier_On == 0: # While barrier has been removed.
                for cell in self.cellList:
                    cell.dict['IsConnectedToMedium'] = 0
                    cell.dict['IsConnectedToSheet'] = 0 
                    cell.dict['EdgeRegion'] = 0
                    cell.dict['EdgeBox'] = 0   
                    
                for cell in self.cellList:
                    if cell.type == 1 or cell.type == 2: 
                        for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell):
                            if neighbor == None: 
                                cell.dict['IsConnectedToMedium'] = 1
                            if neighbor != None:
                                cell.dict['IsConnectedToSheet'] = 1             
                
                for cell in self.cellList:
                    if cell.type == 1 or cell.type == 2: 
                        if cell.dict['IsConnectedToSheet'] == 1 and cell.dict['IsConnectedToMedium'] == 1:
                            cell.dict['EdgeRegion'] = 1
                        else:
                            cell.dict['EdgeRegion'] = 0      
                            
#                 for cell in self.cellList:
#                     if cell.type == 1:
#                         Ref_COM = 0
#                         if cell.dict['IsConnectedToSheet'] == 1 and cell.dict['IsConnectedToMedium'] == 1:
#                             cell.dict['EdgeRegion'] = 1
#                             Ref_COM = cell.xCOM 
#                             for cell in self.cellList:
#                                 if cell.xCOM >= (Ref_COM - 3):
#                                     cell.dict['EdgeRegion'] = 1                            
#                         else:
#                             cell.dict['EdgeRegion'] = 0  
                            
                COM_X_max_Box = 0 # Max COM at every MCS.
                for cell in self.cellList:
                    if cell.type == 1 or cell.type == 2: 
                        if cell.xCOM >= COM_X_max_Box:
                            COM_X_max_Box = cell.xCOM
                
                N_Cells_Box_y_axis = (val_Dimension_y/cell_size) - 2 # 2 is for top and bottom wall cells. Thus, N_Cells_Box_y_axis = 30.
                N_Cells_Box_x_axis = Edge_Box_width_x # width of edge box.      
                N_Cells_NonMoving_Box = N_Cells_Box_x_axis * N_Cells_Box_y_axis           
                Box_size_x = N_Cells_Box_x_axis * cell_size 
                
              
                    
                n_Boxes_x = 0 # number of boxes, each with size of Box_size_x.        
                if (COM_X_max_Box % Box_size_x) == 0:
                    n_Boxes_x = int(COM_X_max_Box // Box_size_x)
                if (COM_X_max_Box % Box_size_x) != 0:
                    n_Boxes_x = int((COM_X_max_Box // Box_size_x) + 1)     
                                
                Densities_Boxes = []  
                i_Box = n_Boxes_x # sarting from the edge box.
                Edge_Box_Found = 0
                while i_Box > 0 and Edge_Box_Found == 0: 
                    N_Cells_per_Box = 0
                    del VelocitiesEdge[:] 
                    del DisplsEdge[:]
                    for cell in self.cellList:
                        if cell.type == 1 or cell.type == 2: 
                            if (Box_size_x * (i_Box-1)) < cell.xCOM < (Box_size_x * i_Box):
                                Medium_Tag = 0
                                Sheet_Tag = 0
                                for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell):
                                    if neighbor != None:
                                        Sheet_Tag = 1
                                if Sheet_Tag == 1: # excludes single cells
                                    N_Cells_per_Box = N_Cells_per_Box + 1
                                                  
                    if N_Cells_per_Box > (N_Cells_NonMoving_Box * 0.5):
                        Edge_Box_Found = 1 
                        for cell in self.cellList:
                            if cell.type == 1 or cell.type == 2: 
                                if (Box_size_x * (i_Box-1)) < cell.xCOM < (Box_size_x * i_Box):
                                    VelocitiesEdge.append(cell.dict['velocityX'])
                                    DisplsEdge.append((cell.xCOM - cell.dict['initial_position_edge_x'])) 
                                    cell.dict['EdgeBox'] = 1                     
                        
                    i_Box = i_Box - 1 
            txt_Matlab_VEdge_Mean.write(str(mcs) + " " + str(np.mean(VelocitiesEdge)) + "\n")   
            txt_Matlab_VEdge_STD.write(str(mcs) + " " + str(np.std(VelocitiesEdge)) + "\n")            
            txt_Matlab_DisplEdge_Mean.write(str(mcs) + " " + str(np.mean(DisplsEdge)) + "\n")   
            txt_Matlab_DisplEdge_STD.write(str(mcs) + " " + str(np.std(DisplsEdge)) + "\n")            
#             if (mcs == 80):
#                 txt_Matlab_Test.write(str(mcs) + ": \n")
#                 L = len(VelocitiesEdge)
#                 txt_Matlab_Test.write('L(VelocitiesEdge) = ' + str(L) + " \n \n")
#                 for i in range(0, (L)):
#                     txt_Matlab_Test.write(str(VelocitiesEdge[i]) + "\n")   
            #.............................................................  
            #==========Area, surface, Eccentricity, Density, Velocity=========
            if (mcs >= mcs_T1 and mcs % 10 == 0):
                COM_X_max = 0 # Max COM at every 10 MCS.
                for cell in self.cellList:
                    if cell.type == 1 or cell.type == 2: 
                        if cell.xCOM >= COM_X_max:
                            COM_X_max = cell.xCOM
            
                Bin_size_x = 20 # Bin size in the x direction is 3 cells wide (every cell is 10*10) for averaging area and surface of the cells in every bin.             
                n_bins = 0 # number of bins each with size of Bin_size_x.        
                if (COM_X_max % Bin_size_x) == 0:
                    n_bins = int(COM_X_max // Bin_size_x)
                if (COM_X_max % Bin_size_x) != 0:
                    n_bins = int((COM_X_max // Bin_size_x) + 1) 
              
                Areas = [] 
                Surfaces = [] 
                Eccentricities = [] 
                Densities = []   
                Velocities = []
                VelocitiesX = []
                Polarities = []
                Forces = []       
                for i_bin in range(1, (n_bins+1)): # For [1, n_bins]: For every bin in x_direction, calculates average area and surface of the cells in that bin.
                    N_Cells_Densities = 0
                    del Areas[:] 
                    del Surfaces[:] 
                    del Eccentricities[:] 
                    del Densities[:]   
                    del Velocities[:]
                    del VelocitiesX[:]
                    del Polarities[:]
                    del Forces[:] 
#                    if mcs == mcs_T4: 
#                         if i_bin == n_bins:
#                             txt_Matlab_Test.write("At MCS = " + str(mcs) + ": \n")
#                             L = len(Velocities)
#                             txt_Matlab_Test.write('L(Velocities[:]) = ' + str(L) + " \n \n")
#                             for i in range(0, (L)):
#                                 txt_Matlab_Test.write(str(Velocities[i]) + "\n") 
#                     Areas.clear()
#                     Surfaces.clear() 
#                     Eccentricities.clear() 
#                     Densities.clear()   
#                     Velocities.clear()
#                     Polarities.clear()
#                     Forces.clear() 
                    for cell in self.cellList:
                        if cell.type == 1 or cell.type == 2: 
                            if (Bin_size_x * (i_bin-1)) < cell.xCOM < (Bin_size_x * i_bin):
                                Areas.append(cell.volume)
                                Surfaces.append(cell.surface)
                                Eccentricities.append(cell.ecc)
                                N_Cells_Densities = N_Cells_Densities + 1
                                if cell.dict['VelocityMagnitude'] != 0:
                                   Velocities.append(cell.dict['VelocityMagnitude'])
                                if cell.dict['velocityX'] != 0:
                                   VelocitiesX.append(cell.dict['velocityX'])                              
                                if cell.dict['PolarityMagnitude'] != 0:
                                   Polarities.append(cell.dict['PolarityMagnitude'])
                                if cell.dict['ForceMagnitude'] != 0:
                                   Forces.append(cell.dict['ForceMagnitude'])

                    txt_Matlab_Area_Mean_T1.write(str(np.mean(Areas)) + " ") 
                    txt_Matlab_Area_STD_T1.write(str(np.std(Areas)) + " ") 
                    #txt_Matlab_Surface_Mean_T1.write(str(i_bin * Bin_size_x) + " " + str(np.mean(Surfaces)) + "\n") 
                    #txt_Matlab_Surface_STD_T1.write(str(i_bin * Bin_size_x) + " " + str(np.std(Surfaces)) + "\n") 
                    #txt_Matlab_Eccentricity_Mean_T1.write(str(i_bin * Bin_size_x) + " " + str(np.mean(Eccentricities)) + "\n") 
                    #txt_Matlab_Eccentricity_STD_T1.write(str(i_bin * Bin_size_x) + " " + str(np.std(Eccentricities)) + "\n") 
                    txt_Matlab_Density_T1.write(str(N_Cells_Densities) + " ")   
                    txt_Matlab_Velocity_Mean_T1.write(str(np.mean(Velocities)) + " ") 
                    txt_Matlab_Velocity_STD_T1.write(str(np.std(Velocities)) + " ")  
                    txt_Matlab_VelocityX_Mean_T1.write(str(np.mean(VelocitiesX)) + " ") 
                    txt_Matlab_VelocityX_STD_T1.write(str(np.std(VelocitiesX)) + " ")                      
                    txt_Matlab_Polarity_Mean_T1.write(str(np.mean(Polarities)) + " ") 
                    txt_Matlab_Polarity_STD_T1.write(str(np.std(Polarities)) + " ") 
                    txt_Matlab_Force_Mean_T1.write(str(np.mean(Forces)) + " ") 
                    txt_Matlab_Force_STD_T1.write(str(np.std(Forces)) + " ")
                    if (mcs == mcs_T1):
                        txt_Matlab_bin.write(str(i_bin * Bin_size_x) + "\n")
                    
                txt_Matlab_Area_Mean_T1.write("\n") 
                txt_Matlab_Area_STD_T1.write("\n") 
                #txt_Matlab_Surface_Mean_T1.write(str(i_bin * Bin_size_x) + " " + str(np.mean(Surfaces)) + "\n") 
                #txt_Matlab_Surface_STD_T1.write(str(i_bin * Bin_size_x) + " " + str(np.std(Surfaces)) + "\n") 
                #txt_Matlab_Eccentricity_Mean_T1.write(str(i_bin * Bin_size_x) + " " + str(np.mean(Eccentricities)) + "\n") 
                #txt_Matlab_Eccentricity_STD_T1.write(str(i_bin * Bin_size_x) + " " + str(np.std(Eccentricities)) + "\n") 
                txt_Matlab_Density_T1.write("\n")   
                txt_Matlab_Velocity_Mean_T1.write("\n") 
                txt_Matlab_Velocity_STD_T1.write("\n")  
                txt_Matlab_VelocityX_Mean_T1.write("\n")  
                txt_Matlab_VelocityX_STD_T1.write("\n")                       
                txt_Matlab_Polarity_Mean_T1.write("\n") 
                txt_Matlab_Polarity_STD_T1.write("\n")  
                txt_Matlab_Force_Mean_T1.write("\n")  
                txt_Matlab_Force_STD_T1.write("\n")                              
            #======================================================          
            sum_V_x_PCells = 0 # For order parameter
            sum_V_y_PCells = 0 # For order parameter
            sum_magnitude_V_PCells = 0 # For order parameter
            OrderV_MCS_PCells = 0 # For order parameter
            sum_V_x_AllCells = 0 # For order parameter
            sum_V_y_AllCells = 0 # For order parameter
            sum_magnitude_V_AllCells = 0 # For order parameter
            OrderV_MCS_AllCells = 0 # For order parameter 
  
            sum_P_x_PCells = 0 # For order parameter
            sum_P_y_PCells = 0 # For order parameter
            sum_magnitude_P_PCells = 0 # For order parameter
            OrderP_MCS_PCells = 0 # For order parameter
            sum_P_x_AllCells = 0 # For order parameter
            sum_P_y_AllCells = 0 # For order parameter
            sum_magnitude_P_AllCells = 0 # For order parameter
            OrderP_MCS_AllCells = 0 # For order parameter 

            sum_F_x_PCells = 0 # For order parameter
            sum_F_y_PCells = 0 # For order parameter
            sum_magnitude_F_PCells = 0 # For order parameter
            OrderF_MCS_PCells = 0 # For order parameter
            sum_F_x_AllCells = 0 # For order parameter
            sum_F_y_AllCells = 0 # For order parameter
            sum_magnitude_F_AllCells = 0 # For order parameter
            OrderF_MCS_AllCells = 0 # For order parameter 
              
            N_COMX = 0 
            for cell in self.cellList:
                if cell.type == 1 or cell.type == 2:
                    dpdtX = 0
                    dpdtY = 0           
                
                    Cells_str = ''
                    #Cells_str += "ID: " + str(cell.id) + "; " + "P/UP = " + str(cell.type) + "; " + "Pre_COM = (" + str(round(cell.dict['PreviousComX'], 3)) + "," + str(round(cell.dict['PreviousComY'], 3)) + "); " + "Pre_V = (" + str(round(cell.dict['velocityX'], 3)) + "," + str(round(cell.dict['velocityY'], 3)) + "); " + "Pre_|V| = " + str(round(cell.dict['VelocityMagnitude'], 3)) + "; " + "Pre_V_Angle = " + str(round(cell.dict['VelocityAngle'], 3)) + "; " + "Pre_P = (" + str(round(cell.dict['PolarityX'], 3)) + "," + str(round(cell.dict['PolarityY'], 3)) + "); " + "Pre_|P| = " + str(round(cell.dict['PolarityMagnitude'], 3)) + "; " + "PreP_Angle = " + str(round(cell.dict['PolarityAngle'], 3)) + "; " + "Pre_F = (" + str(round(cell.lambdaVecX, 3)) + "," + str(round(cell.lambdaVecY, 3)) + "); " + "Pre_|F| = " + str(round(cell.dict['ForceMagnitude'], 3)) + "; " + "PreF_Angle = " + str(round(cell.dict['ForceAngle'], 3)) + " \n"
                                      
                    txt_Matlab_COMx.write(str(round(cell.dict['PreviousComX'], 3)) + " ")   
                    txt_Matlab_COMy.write(str(round(cell.dict['PreviousComY'], 3)) + " ") 
                    txt_Matlab_V.write(str(round(cell.dict['VelocityMagnitude'], 3)) + " ") 
#                    txt_Matlab_VA.write(str(round(cell.dict['VelocityAngle'], 3)) + " ")                   
                    txt_Matlab_P.write(str(round(cell.dict['PolarityMagnitude'], 3)) + " ") 
#                     txt_Matlab_PA.write(str(round(cell.dict['PolarityAngle'], 3)) + " ")
                    if mcs > (MCS_polarizing_edge + 1):
                        if cell.type == 1:
                            txt_Matlab_PA.write(str(round(cell.dict['PolarityAngle'], 3)) + " ")
                            txt_Matlab_VA.write(str(round(cell.dict['VelocityAngle'], 3)) + " ")    
                            txt_Matlab_FA.write(str(round(cell.dict['ForceAngle'], 3)) + " ")
                            if cell.dict['EdgeRegion'] == 1: 
                                txt_Matlab_PA_Edge.write(str(round(cell.dict['PolarityAngle'], 3)) + " ")  
#                    if cell.type == 1:
#                        txt_Matlab_PA.write("[" + str(cell.id) + ": " + str(round(cell.dict['PolarityAngle'], 3)) + "]; ") 
                    txt_Matlab_F.write(str(round(cell.dict['ForceMagnitude'], 3)) + " ")
#                    txt_Matlab_FA.write(str(round(cell.dict['ForceAngle'], 3)) + " ") 
#                    if cell.type == 1:
#                        txt_Matlab_FA.write("[" + str(cell.id) + ": " + str(round(cell.dict['ForceAngle'], 3)) + "]; ") 
#                        txt_Matlab_FX.write("[" + str(cell.id) + ": " + str(round(cell.lambdaVecX, 3)) + "]; ") 
#                        txt_Matlab_FY.write("[" + str(cell.id) + ": " + str(round(cell.lambdaVecY, 3)) + "]; ") 
                    txt_Matlab_FX.write(str(round(cell.lambdaVecX, 3)) + " ") 
                    txt_Matlab_FX.write(str(round(cell.lambdaVecY, 3)) + " ") 

                    N_COMX = N_COMX + 1
                   
                    cell.dict['velocityX'] = cell.xCOM - cell.dict['PreviousComX'] #Vx = comX - previousComX
                    cell.dict['velocityY'] = cell.yCOM - cell.dict['PreviousComY'] #Vx = comX - previousComX
                    cell.dict['VelocityMagnitude'] = math.sqrt(cell.dict['velocityX']**2 + cell.dict['velocityY']**2)
                   
                    sum_magnitude_V_AllCells = sum_magnitude_V_AllCells + cell.dict['VelocityMagnitude']
                    sum_V_x_AllCells = sum_V_x_AllCells + cell.dict['velocityX'] 
                    sum_V_y_AllCells = sum_V_y_AllCells + cell.dict['velocityY']                  
                    if cell.type == 1: 
                        sum_magnitude_V_PCells = sum_magnitude_V_PCells + cell.dict['VelocityMagnitude']
                        sum_V_x_PCells = sum_V_x_PCells + cell.dict['velocityX'] 
                        sum_V_y_PCells = sum_V_y_PCells + cell.dict['velocityY']                      
                   
                    if cell.dict['VelocityMagnitude'] == 0:
                        cell.dict['VelocityAngle'] = 0
                    else:
                        if cell.dict['velocityY'] > 0:
                            cell.dict['VelocityAngle'] = math.degrees(math.acos(cell.dict['velocityX'] / cell.dict['VelocityMagnitude']))
                        if cell.dict['velocityY'] < 0:
                            cell.dict['VelocityAngle'] = (-1) * math.degrees(math.acos(cell.dict['velocityX'] / cell.dict['VelocityMagnitude']))
                        if cell.dict['velocityY'] == 0:
                            if cell.dict['velocityX'] > 0:
                                cell.dict['VelocityAngle'] = 0
                            if cell.dict['velocityX'] < 0:
                                cell.dict['VelocityAngle'] = 180                
               
                    dpdtX = (-1) * cell.dict['beta'] * cell.dict['PolarityX'] + cell.dict['gamma'] * cell.dict['velocityX'] #dP/dt_x = -Beta Px + Gamma Vx
                    dpdtY = (-1) * cell.dict['beta'] * cell.dict['PolarityY'] + cell.dict['gamma'] * cell.dict['velocityY'] #dP/dt_y = -Beta Py + Gamma Vy
                                  
                    cell.dict['PreviousComX'] = cell.xCOM # previousComX = Xcom
                    cell.dict['PreviousComY'] = cell.yCOM # previousComY = Ycom                                       
                   
                    cell.dict['PolarityX'] = cell.dict['PolarityX'] + dpdtX # PolarityX = PolarityX + dP/dt
                    cell.dict['PolarityY'] = cell.dict['PolarityY'] + dpdtY # PolarityX = PolarityX + dP/dt  
                   
                    cell.dict['PolarityMagnitude'] = math.sqrt(cell.dict['PolarityX']**2 + cell.dict['PolarityY']**2) # with new polarity
                   #--------------- Polarity angle --------------------                   
                    if Noise_Angle == 0:
                        if cell.dict['PolarityMagnitude'] == 0:
                            cell.dict['PolarityAngle'] = 0
                        else:
                            if cell.dict['PolarityY'] > 0:
                                cell.dict['PolarityAngle'] = math.degrees(math.acos(cell.dict['PolarityX'] / cell.dict['PolarityMagnitude']))
                            if cell.dict['PolarityY'] < 0:
                                cell.dict['PolarityAngle'] = (-1) * math.degrees(math.acos(cell.dict['PolarityX'] / cell.dict['PolarityMagnitude']))
                            if cell.dict['PolarityY'] == 0:
                                if cell.dict['PolarityX'] > 0:
                                    cell.dict['PolarityAngle'] = 0
                                if cell.dict['PolarityX'] < 0:
                                    cell.dict['PolarityAngle'] = 180
                                   
                        if cell.dict['EdgeRegion'] == 1:
                            cell.dict['PolarityMagnitude'] = Initial_Polarization
                            cell.dict['PolarityX'] = Initial_Polarization
                            cell.dict['PolarityY'] = 0
                            cell.dict['PolarityAngle'] = 0                            
                       
                    if Noise_Angle > 0:
                        if cell.type == 2: # White noise is not added to the polarity direction of unpolarized cells.
                            if cell.dict['PolarityMagnitude'] == 0:
                                cell.dict['PolarityAngle'] = 0
                            else:
                                if cell.dict['PolarityY'] > 0:
                                    cell.dict['PolarityAngle'] = math.degrees(math.acos(cell.dict['PolarityX'] / cell.dict['PolarityMagnitude']))
                                if cell.dict['PolarityY'] < 0:
                                    cell.dict['PolarityAngle'] = (-1) * math.degrees(math.acos(cell.dict['PolarityX'] / cell.dict['PolarityMagnitude']))
                                if cell.dict['PolarityY'] == 0:
                                    if cell.dict['PolarityX'] > 0:
                                        cell.dict['PolarityAngle'] = 0
                                    if cell.dict['PolarityX'] < 0:
                                        cell.dict['PolarityAngle'] = 180
                            
                            if cell.dict['EdgeRegion'] == 1:
                                cell.dict['PolarityMagnitude'] = Initial_Polarization
                                cell.dict['PolarityX'] = Initial_Polarization
                                cell.dict['PolarityY'] = 0
                                cell.dict['PolarityAngle'] = 0           
                                        
                        if cell.type == 1: # Add white noise to the polarity direction of polarized cells.
                            if cell.dict['EdgeRegion'] == 0:  # Just for cells within the sheet, to avoid double-update of the angle of edge cells.  
#                                 if cell.dict['EdgeBox'] == 1:
#                                     if mcs > 3000:
#                                         cell.dict['PolarityMagnitude'] = EdgeBoxPolarity
                                RanAngle = (Noise_Angle * ra.uniform(-1,1)) # a random rotation angle in [-Noise_Angle, Noise_Angle].
                                dAngledt = (-1)* cell.dict['beta'] * cell.dict['PolarityAngle'] + RanAngle #dAngle/dt = -Beta(Angle - Angle_0) + rnd_angle.
                                cell.dict['PolarityAngle'] = cell.dict['PolarityAngle'] + dAngledt
                                if cell.dict['PolarityAngle'] > 180:
                                    cell.dict['PolarityAngle'] = -180 + (cell.dict['PolarityAngle'] - 180)
                                if cell.dict['PolarityAngle'] < -180:
                                    cell.dict['PolarityAngle'] = 180 + (180 + cell.dict['PolarityAngle'])
                                cell.dict['PolarityX'] = (math.cos(math.radians(cell.dict['PolarityAngle'])) * cell.dict['PolarityX']) - (math.sin(math.radians(cell.dict['PolarityAngle'])) * cell.dict['PolarityY']) #P2x = cos(theta)*P1x - sin(theta)*P1y;
                                cell.dict['PolarityY'] = (math.sin(math.radians(cell.dict['PolarityAngle'])) * cell.dict['PolarityX']) + (math.cos(math.radians(cell.dict['PolarityAngle'])) * cell.dict['PolarityY']) #P2y = sin(theta)*P1x + cos(theta)*P1y;            
                            if (cell.dict['EdgeRegion'] == 1 and MaintainPolarization_LeadingCells == 1):
                                cell.dict['PolarityMagnitude'] = Initial_Polarization  
#                                 if cell.dict['EdgeBox'] == 1:
#                                     if mcs > 3000:
#                                         cell.dict['PolarityMagnitude'] = EdgeBoxPolarity
                                RanAngle = (Noise_Angle * ra.uniform(-1,1)) # a random rotation angle in [-EtaAngle, EtaAngle].
                                dAngledt = (-1)* cell.dict['beta'] * cell.dict['PolarityAngle'] + RanAngle #dAngle/dt = -Beta(Angle-Angle_0) + rnd_angle.
                                cell.dict['PolarityAngle'] = cell.dict['PolarityAngle'] + dAngledt
                                if cell.dict['PolarityAngle'] > 180:
                                    cell.dict['PolarityAngle'] = -180 + (cell.dict['PolarityAngle'] - 180)
                                if cell.dict['PolarityAngle'] < -180:
                                    cell.dict['PolarityAngle'] = 180 + (180 + cell.dict['PolarityAngle'])
                                cell.dict['PolarityX'] = math.cos(math.radians(cell.dict['PolarityAngle'])) * cell.dict['PolarityMagnitude']  
                                cell.dict['PolarityY'] = math.sin(math.radians(cell.dict['PolarityAngle'])) * cell.dict['PolarityMagnitude']                     
#                   cell.dict['PolarityAngle'] = ra.uniform(-1,1) * 60 # This is a random angle in degrees.
#                   cell.dict['PolarityX'] = cell.dict['PolarityMagnitude'] * math.cos(math.radians(cell.dict['PolarityAngle']))
#                   cell.dict['PolarityY'] = cell.dict['PolarityMagnitude'] * math.sin(math.radians(cell.dict['PolarityAngle'])) 
                   #---------------------------------------------------
                    if cell.dict['PolarityAngle'] >= max_P_angle:
                        max_P_angle = cell.dict['PolarityAngle']
                    if cell.dict['PolarityAngle'] <= min_P_angle:
                        min_P_angle = cell.dict['PolarityAngle']                 
                    sum_magnitude_P_AllCells = sum_magnitude_P_AllCells + cell.dict['PolarityMagnitude']
                    sum_P_x_AllCells = sum_P_x_AllCells + cell.dict['PolarityX'] 
                    sum_P_y_AllCells = sum_P_y_AllCells + cell.dict['PolarityY']                  
                    if cell.type == 1: 
                        sum_magnitude_P_PCells = sum_magnitude_P_PCells + cell.dict['PolarityMagnitude']
                        sum_P_x_PCells = sum_P_x_PCells + cell.dict['PolarityX'] 
                        sum_P_y_PCells = sum_P_y_PCells + cell.dict['PolarityY']                   
                   
                    if cell.dict['PolarityMagnitude'] == 0:
                        cell.lambdaVecX = 0
                        cell.lambdaVecY = 0
                    else:  
                        cell.lambdaVecX = (cell.dict['Fmax'] * (cell.dict['PolarityX'] / cell.dict['PolarityMagnitude'])) * ((cell.dict['PolarityMagnitude']**cell.dict['n'])/((cell.dict['PolarityMagnitude']**cell.dict['n']) + (cell.dict['alpha']**cell.dict['n'])))# Fx = (Px/|P|)(|P|^n/[|P|^n + alpha^n])
                        cell.lambdaVecY = (cell.dict['Fmax'] * (cell.dict['PolarityY'] / cell.dict['PolarityMagnitude'])) * ((cell.dict['PolarityMagnitude']**cell.dict['n'])/((cell.dict['PolarityMagnitude']**cell.dict['n']) + (cell.dict['alpha']**cell.dict['n'])))# Fy = (Py/|P|)(|P|^n/[|P|^n + alpha^n])
                                                                                  
                    cell.dict['ForceMagnitude'] = math.sqrt(cell.lambdaVecX**2 + cell.lambdaVecY**2)           
                    if cell.dict['ForceMagnitude'] == 0:
                        cell.dict['ForceAngle'] = 0
                    else:
                        if cell.lambdaVecY < 0:
                            cell.dict['ForceAngle'] = math.degrees(math.acos(((-1) * cell.lambdaVecX) / cell.dict['ForceMagnitude'])) 
                        if cell.lambdaVecY > 0:
                            cell.dict['ForceAngle'] = (-1) * math.degrees(math.acos(((-1) * cell.lambdaVecX) / cell.dict['ForceMagnitude'])) 
                        if cell.lambdaVecY == 0:
                            if cell.lambdaVecX < 0:
                                cell.dict['ForceAngle'] = 0
                            if cell.lambdaVecX > 0:
                                cell.dict['ForceAngle'] = 180                    
         
                    sum_magnitude_F_AllCells = sum_magnitude_F_AllCells + cell.dict['ForceMagnitude']
                    sum_F_x_AllCells = sum_F_x_AllCells + cell.lambdaVecX 
                    sum_F_y_AllCells = sum_F_y_AllCells + cell.lambdaVecY                  
                    if cell.type == 1: 
                        sum_magnitude_F_PCells = sum_magnitude_F_PCells + cell.dict['ForceMagnitude']
                        sum_F_x_PCells = sum_F_x_PCells + cell.lambdaVecX  
                        sum_F_y_PCells = sum_F_y_PCells + cell.lambdaVecY  
         
                    Cells_str += "ID:" + str(cell.id) + "; " + "P/UP=" + str(cell.type) + "; " + "COM=(" + str(round(cell.dict['PreviousComX'], 3)) + "," + str(round(cell.dict['PreviousComY'], 3)) + "); " + "P_reached_end = " + str(Polarity_reached_end) + "; " + "V=(" + str(round(cell.dict['velocityX'], 3)) + "," + str(round(cell.dict['velocityY'], 3)) + "); " + "|V|=" + str(round(cell.dict['VelocityMagnitude'], 3)) + "; " + "V_Angle=" + str(round(cell.dict['VelocityAngle'], 3)) + "; " + "dP/dt=(" + str(round(dpdtX, 3)) + "," + str(round(dpdtY, 3)) + "); " + "P=(" + str(round(cell.dict['PolarityX'], 3)) + "," + str(round(cell.dict['PolarityY'], 3)) + "); " + "|P|=" + str(round(cell.dict['PolarityMagnitude'], 3)) + "; " + "P_Angle=" + str(round(cell.dict['PolarityAngle'], 3)) + "; " + "F=(" + str(round(cell.lambdaVecX, 3)) + "," + str(round(cell.lambdaVecY, 3)) + "); " + "|F|=" + str(round(cell.dict['ForceMagnitude'], 3)) + "; "+ "F_Angle=" + str(round(cell.dict['ForceAngle'], 3)) + ";" + "\n"
                    txt_f.write(Cells_str) # I write "Cells_str" in the file.           
                   
                    if cell.type == 2: 
                        if abs(cell.dict['PolarityMagnitude']) >= cell.dict['P_threshold']:
                            cell.type = 1  
                   
#                     if cell.type == 1: 
#                         if abs(cell.dict['PolarityMagnitude']) < cell.dict['P_threshold']:
#                             cell.type = 2         
            
            #if (mcs > MCS_polarizing_edge): # Barrier has been removed, there are polarised cells. So sum_magnitude_V_PCells is not zero.
            if (mcs > MCS_Closed_Edge_cells_Zero_F):
                if (sum_magnitude_V_PCells != 0):# When there are polarised cells.
                    OrderV_MCS_PCells = math.sqrt(sum_V_x_PCells**2 + sum_V_y_PCells**2)/sum_magnitude_V_PCells
            txt_Matlab_OV_PCells.write(str(mcs) + " " + str(round(OrderV_MCS_PCells, 3)) + "\n")           
            if (mcs > MCS_Closed_Edge_cells_Zero_F):
                if (sum_magnitude_V_AllCells != 0):
                    OrderV_MCS_AllCells = math.sqrt(sum_V_x_AllCells**2 + sum_V_y_AllCells**2)/sum_magnitude_V_AllCells
            txt_Matlab_OV_AllCells.write(str(mcs) + " " + str(round(OrderV_MCS_AllCells, 3)) + "\n")
            
            #if (mcs > MCS_polarizing_edge): # Barrier has been removed, there are polarised cells. So sum_magnitude_V_PCells is not zero.
            if (mcs > MCS_Closed_Edge_cells_Zero_F):
                if (sum_magnitude_P_PCells != 0): # When there are polarised cells.
                    OrderP_MCS_PCells = math.sqrt(sum_P_x_PCells**2 + sum_P_y_PCells**2)/sum_magnitude_P_PCells
            txt_Matlab_OP_PCells.write(str(mcs) + " " + str(round(OrderP_MCS_PCells, 3)) + "\n")           
            if (mcs > MCS_Closed_Edge_cells_Zero_F):
                if (sum_magnitude_P_AllCells != 0):
                    OrderP_MCS_AllCells = math.sqrt(sum_P_x_AllCells**2 + sum_P_y_AllCells**2)/sum_magnitude_P_AllCells
            txt_Matlab_OP_AllCells.write(str(mcs) + " " + str(round(OrderP_MCS_AllCells, 3)) + "\n")
            
            #if (mcs > MCS_polarizing_edge): # Barrier has been removed, there are polarised cells. So sum_magnitude_V_PCells is not zero.
            if (mcs > MCS_Closed_Edge_cells_Zero_F):
                if (sum_magnitude_F_PCells != 0):# When there are polarised cells.
                    OrderF_MCS_PCells = math.sqrt(sum_F_x_PCells**2 + sum_F_y_PCells**2)/sum_magnitude_F_PCells
            txt_Matlab_OF_PCells.write(str(mcs) + " " + str(round(OrderF_MCS_PCells, 3)) + "\n")           
            
            if (mcs > MCS_Closed_Edge_cells_Zero_F):
                if (sum_magnitude_F_AllCells != 0):
                    OrderF_MCS_AllCells = math.sqrt(sum_F_x_AllCells**2 + sum_F_y_AllCells**2)/sum_magnitude_F_AllCells
            txt_Matlab_OF_AllCells.write(str(mcs) + " " + str(round(OrderF_MCS_AllCells, 3)) + "\n")
            
            COM_X_MAX = 0 # Max COM for plotting.
            
            for cell in self.cellList:
                if cell.type == 1 or cell.type == 2: 
                    if cell.xCOM >= COM_X_MAX: 
                        COM_X_MAX = cell.xCOM
                
                        
            txt_Matlab_MAXCOM_x.write(str(COM_X_MAX) + "\n")
            
            
        
            
            
            txt_Matlab_PA_MinMax.write("Max polarity angle = " + str(round(max_P_angle, 3)) + "\n" + "Min polarity angle = " + str(round(min_P_angle, 3)) + "\n")

            txt_Matlab_COMx.write("\n")
            txt_Matlab_COMy.write("\n")
            txt_Matlab_V.write("\n")
#            txt_Matlab_VA.write("\n")
            txt_Matlab_P.write("\n")
#             txt_Matlab_PA.write("\n")
#             if mcs > (MCS_polarizing_edge + 1):
#                 txt_Matlab_PA_Edge.write("\n")
            txt_Matlab_PA_MinMax.write("\n")
            txt_Matlab_F.write("\n")
#            txt_Matlab_FA.write("\n") 
            txt_Matlab_FX.write("\n") 
            txt_Matlab_FY.write("\n")          
            txt_Matlab_NCOMx.write(str(N_COMX) + "\n")
           
            txt_f.close() 
            txt_Matlab_COMx.close()
            txt_Matlab_COMy.close()
            txt_Matlab_V.close()   
            txt_Matlab_VA.close()  
            txt_Matlab_OV_PCells.close()
            txt_Matlab_OV_AllCells.close() 
            txt_Matlab_P.close()   
            txt_Matlab_PA_Edge.close()
            txt_Matlab_PA.close()
            txt_Matlab_PA_MinMax.close()
            txt_Matlab_OP_PCells.close()
            txt_Matlab_OP_AllCells.close()   
            txt_Matlab_F.close()   
            txt_Matlab_FA.close()
            txt_Matlab_FX.close()
            txt_Matlab_FY.close()
            txt_Matlab_OF_PCells.close()
            txt_Matlab_OF_AllCells.close()
            txt_Matlab_Np.close()
            #txt_Matlab_VEdge.close()
            txt_Matlab_VEdge_Mean.close()
            txt_Matlab_VEdge_STD.close()
            txt_Matlab_DisplEdge_Mean.close()
            txt_Matlab_DisplEdge_STD.close()
            txt_Matlab_NCOMx.close()           
            txt_Matlab_Area_Mean_T1.close()             
            txt_Matlab_Area_STD_T1.close()   
              
            txt_Matlab_Surface_Mean_T1.close() 
            txt_Matlab_Surface_STD_T1.close() 
            
            txt_Matlab_Eccentricity_Mean_T1.close() 
            txt_Matlab_Eccentricity_STD_T1.close() 
            
            txt_Matlab_Density_T1.close() 
           
            txt_Matlab_Velocity_Mean_T1.close() 
            txt_Matlab_Velocity_STD_T1.close() 
            
            txt_Matlab_VelocityX_Mean_T1.close() 
            txt_Matlab_VelocityX_STD_T1.close() 
             
            txt_Matlab_Polarity_Mean_T1.close() 
            txt_Matlab_Polarity_STD_T1.close() 
                       
            txt_Matlab_Force_Mean_T1.close() 
            txt_Matlab_Force_STD_T1.close() 
            
            txt_Matlab_bin.close()
            
            txt_Matlab_MAXCOM_x.close()
            
        #========================================================================================            
        pass
    
    #...........................
    #...........................
    #...........................    
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        
#===================Velocity Vector Field======================        
from PySteppables import *
import CompuCell
import sys
from PlayerPython import *
import CompuCellSetup
from math import *
from PlayerPython import * 
class VectorFieldCellLevelVisualizationSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)          
    def setVectorField(self,_field): #I added this function from the manual file for Python scripting 3.6 page 21. 
        self.vectorField=_field
    def step(self,mcs):
        clearVectorCellLevelField(self.vectorField)
        for cell in self.cellList:
            if cell.type==1 or cell.type==2:
                insertVectorIntoVectorCellLevelField(self.vectorField, cell, cell.dict['velocityX'], cell.dict['velocityY'], 0.0)
    def finish(self):
        return        
#================================================================
#--------------------Velocity magnitude--------------------------
from PySteppables import *
import CompuCell
import sys
from PlayerPython import *
import CompuCellSetup
from math import *
class ScalarFieldVisualizationSteppableV(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)   
    def setScalarField(self,_field):
        self.scalarField=_field 
    def step(self,mcs):
        clearScalarValueCellLevel(self.scalarField)
        for cell in self.cellList:
            if cell.type==1 or cell.type==2:
                fillScalarValueCellLevel(self.scalarField,cell,cell.dict['VelocityMagnitude'])       
    def finish(self):
        return
#--------------------------------------------------------------------
#*************************Polarity magnitude*************************
from PySteppables import *
import CompuCell
import sys
from PlayerPython import *
import CompuCellSetup
from math import *
class ScalarFieldVisualizationSteppableP(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)   
    def setScalarField(self,_field):
        self.scalarField=_field 
    def step(self,mcs):
        clearScalarValueCellLevel(self.scalarField)
        for cell in self.cellList:
            if cell.type==1 or cell.type==2:
                fillScalarValueCellLevel(self.scalarField,cell,cell.dict['PolarityMagnitude'])       
    def finish(self):
        return
#**********************************************************************    
#...........................Polarity angle.........................
from PySteppables import *
import CompuCell
import sys
from PlayerPython import *
import CompuCellSetup
from math import *
class ScalarFieldVisualizationSteppablePA(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)   
    def setScalarField(self,_field):
        self.scalarField=_field 
    def step(self,mcs):
        clearScalarValueCellLevel(self.scalarField)
        for cell in self.cellList:
            if cell.type==1 or cell.type==2:
                fillScalarValueCellLevel(self.scalarField,cell,cell.dict['PolarityAngle'])       
    def finish(self):
        return
#......................................................................... 
#*************************Force magnitude*************************
from PySteppables import *
import CompuCell
import sys
from PlayerPython import *
import CompuCellSetup
from math import *
class ScalarFieldVisualizationSteppableF(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)   
    def setScalarField(self,_field):
        self.scalarField=_field 
    def step(self,mcs):
        clearScalarValueCellLevel(self.scalarField)
        for cell in self.cellList:
            if cell.type==1 or cell.type==2:
                fillScalarValueCellLevel(self.scalarField,cell,cell.dict['ForceMagnitude'])       
    def finish(self):
        return
#**********************************************************************    
#...........................Force angle.........................
from PySteppables import *
import CompuCell
import sys
from PlayerPython import *
import CompuCellSetup
from math import *
class ScalarFieldVisualizationSteppableFA(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)   
    def setScalarField(self,_field):
        self.scalarField=_field 
    def step(self,mcs):
        clearScalarValueCellLevel(self.scalarField)
        for cell in self.cellList:
            if cell.type==1 or cell.type==2:
                fillScalarValueCellLevel(self.scalarField,cell,cell.dict['ForceAngle'])       
    def finish(self):
        return
#......................................................................... 
#----------------------------Force X component-----------------------
from PySteppables import *
import CompuCell
import sys
from PlayerPython import *
import CompuCellSetup
from math import *
class ScalarFieldVisualizationSteppableFx(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)   
    def setScalarField(self,_field):
        self.scalarField=_field 
    def step(self,mcs):
        clearScalarValueCellLevel(self.scalarField)
        for cell in self.cellList:
            if cell.type==1 or cell.type==2:
                fillScalarValueCellLevel(self.scalarField,cell,cell.lambdaVecX)       
    def finish(self):
        return
#......................................................................... 
#----------------------------Force Y component-----------------------
from PySteppables import *
import CompuCell
import sys
from PlayerPython import *
import CompuCellSetup
from math import *
class ScalarFieldVisualizationSteppableFy(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)   
    def setScalarField(self,_field):
        self.scalarField=_field 
    def step(self,mcs):
        clearScalarValueCellLevel(self.scalarField)
        for cell in self.cellList:
            if cell.type==1 or cell.type==2:
                fillScalarValueCellLevel(self.scalarField,cell,cell.lambdaVecY)       
    def finish(self):
        return
#......................................................................... 
from PySteppables import *
import CompuCell
import sys
import random

from PySteppablesExamples import MitosisSteppableBase

G_beta = 0.04

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
        # 0 - parent child position will be randomized between mitosis event
        # negative integer - parent appears on the 'left' of the child
        # positive integer - parent appears on the 'right' of the child
        self.setParentChildPositionFlag(1)
        
    
    def step(self,mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        File_path = 'C:\Users\qq871\Desktop\Multi\ParameterScan\Alpha1'
        FileNumCell_path = os.path.join(File_path, "MatlabNumCell.txt")
        txt_Matlab_Num_Cell = open(FileNumCell_path, "a")
        
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        MCS_Barrier_removal = (MCS_Closed_Edge_cells_Zero_F + MCS_Closed_Edge_cells_NonZero_F) # Number of MCS with: barrier is removed and cell can apply force (cell polarity is linked to F).        
        MCS_polarizing_edge = MCS_Barrier_removal + 1
        numCell = 0
        G_count = 0
        
        if mcs > MCS_Closed_Edge_cells_Zero_F:
            for cell in self.cellList:
                if cell.type == 1 or cell.type == 2:
                    numCell +=1
                    if cell.volume >= 100:
                        G_count +=1
            
            txt_Matlab_Num_Cell.write(str(mcs) + " " + str(numCell) + "\n")

        totalG = G_count * G_beta
        if mcs > MCS_polarizing_edge:
            for cell in self.cellList:
                if (cell.type == 1 or cell.type == 2):
                    
                    if random.random()<=(G_beta/totalG) and cell.volume >= 100:
                        cell.dict['comXbeforePro'] = cell.xCOM
                        cell.dict['comYbeforePro'] = cell.yCOM
                        cells_to_divide.append(cell)
            
            for cell in cells_to_divide:
            # to change mitosis mode leave one of the below lines uncommented
                #if cell.dict['PolarityAngle'] == 180 or cell.dict['PolarityAngle'] == 0:
                    #self.divideCellOrientationVectorBased(cell,1,0,0)
                #else:
                    #self.divideCellOrientationVectorBased(cell,cell.dict['PolarityMagnitude']*math.cos(cell.dict['PolarityAngle']),cell.dict['PolarityMagnitude']*math.sin(cell.dict['PolarityAngle']),0)
                self.divideCellRandomOrientation(cell)
            # self.divideCellOrientationVectorBased(cell,1,0,0)                 # this is a valid option
            # self.divideCellAlongMajorAxis(cell)                               # this is a valid option
            # self.divideCellAlongMinorAxis(cell)                               # this is a valid option
        
        txt_Matlab_Num_Cell.close()
    def updateAttributes(self):
        
        self.cloneParent2Child()
        changeComX = self.parentCell.xCOM - self.parentCell.dict['comXbeforePro']
        changeComY = self.parentCell.yCOM - self.parentCell.dict['comYbeforePro']
        self.parentCell.dict['PreviousComX'] = self.parentCell.dict['PreviousComX'] + changeComX
        self.parentCell.dict['PreviousComY'] = self.parentCell.dict['PreviousComY'] + changeComY
        self.childCell.dict['PolarityX'] = 0   
        self.childCell.dict['PolarityY'] = 0 
        self.childCell.dict['PolarityMagnitude'] = 0
        self.childCell.dict['PolarityAngle'] = 0
        #self.childCell.lambdaVecX = 0
        #self.childCell.lambdaVecY = 0
        #self.childCell.dict['ForceMagnitude'] = 0
        #self.childCell.dict['ForceAngle'] = 0
        
        self.childCell.targetVolume = self.childCell.volume
        self.parentCell.targetVolume = self.parentCell.volume
        
        self.childCell.dict['PreviousComX'] = self.childCell.xCOM - self.parentCell.dict['velocityX']
        self.childCell.dict['PreviousComY'] = self.childCell.yCOM - self.parentCell.dict['velocityY']

        
        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.cloneAttributes(sourceCell=self.parentCell, targetCell=self.childCell, no_clone_key_dict_list = [attrib1, attrib2] )
    
    def finish(self):
        return
        
#......................................................................... 
from PySteppables import *
import CompuCell
import sys
import random

from PySteppablesExamples import SteppableBasePy      
class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        
        for cell in self.cellList:
            cell.targetVolume=100
            cell.lambdaVolume=70
            
            
#......................................................................... 
from PySteppables import *
import CompuCell
import sys
import random

from PySteppablesExamples import SteppableBasePy

class GrowthSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        MCS_Barrier_removal = (MCS_Closed_Edge_cells_Zero_F + MCS_Closed_Edge_cells_NonZero_F) # Number of MCS with: barrier is removed and cell can apply force (cell polarity is linked to F).        
        MCS_polarizing_edge = MCS_Barrier_removal + 1
        #if mcs > MCS_polarizing_edge:
        #    for cell in self.cellList:
        #        if cell.type == 1 or cell.type == 2:
        #            if cell.targetVolume < 100:
        #                cell.targetVolume +=1

