## Project: Mathematical Collective Cell Migration in Two-Dimension.

This project is my work for Summer Research Program 2019/202 offered by UQ. I developed this model under the supervision of Dr. Zoltan Neufeld and Dr. Hamid Khataee.

## What is this model?

Collective motion of cells is one of the principal modes of cell migration. Unlike single-cell modelling, moving collectively is more sophisticated. Biological experiments have shown that cellular motility arises from a few fundamental elements: polarity in response to the presence of a cell free region; adhesive interactions with adjacent cells; growth and proliferation. However, our understanding for how the interplay of the collective dynamics of multicellular systems attributes the emergence of cell movement is far from complete. In this paper, we investigate this influence using a computational model known as the Cellular Potts Model, incorporating the dynamics of polarisation and the mechanism of proliferation. The key finding presented by the model suggests that cell growth after proliferation impacts the collective migration pattern by inducing spontaneous polarity, which then forms swirling movements. Contractile tension, on the other hand, influences the magnitude of such pattern.

Here we model motile force as a Hill function dependent on polarity as follows:

################################################################################################
     velocity v = (dCOM_x/dt, dCOM_y/dt) 
     dp/dt = -Beta p + Gamma V    (in x and y directions)
    
     Centre of mass (COM_x,COM_y) COM = COM + v dt             (Euler's forward method in x and y directions)               
     polarity p = p + dp/dt dt
   
     Force F = Fmax (p/|p|)(|p|^n/[|p|^n + Alpha^n])   (in x and y directions, i.e. lambdaVecX and lambdaVecY in CompuCell3D)

     Force F (applied to a cell) is calculated for all cells. Then, F is pluged into energy function below.

                           Area elasticity       Contractility(actin)       Adhesion(E-cadherin)         Motile force
     Potts energy = SUM_AllCells_i[Ai - At]^2 + SUM_AllCells_i[Li - Lt]^2 - SUM_AllCells_i[Contact_i] +SUM_AllCells_i [F.s]
     
     Proliferation law G(A)= μ[(1/2) tanh⁡(270(√A  -(8.5√(At))/1000))+1]
     Growth A = (MCS - MCS_pro) + A_pro   for A < At
###############################################################################################
# This model takes in the approach of collecitve migration from Oelz et al. (2019) (link:https://link.aps.org/doi/10.1103/PhysRevE.100.032403)
 
# and proliferative mechanics from Baker et al. (2018) (link:https://doi.org/10.1016/j.jtbi.2018.12.025)
