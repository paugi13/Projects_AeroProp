# Projects_AeroProp

PR (220028) G06. Two separated folders for propulsion and aerodynamic calculations. Here you have how this repository works, what the most important scripts do and where to find them.
Every subdepartment has some notes folders where the most important theory knowledge serves as documentation. 

# Aerodynamics

This part of the repo is divided in three main blocks:

 - Airfoil design: Folder with the code to analyse one 4-digit NACA airfoil deeply and to compare characteristics. 
 - Airplane configuration: This folder aims to store the scripts that are used to trim the airplane once the previous studies have been correctly carried out.
 - Wing analysis: Necessary scripts using LLWING software to design the optimum wing for the plane. 

## Airfoil design



## Airplane configuration



## Wing analysis



# Propulsion

The propulsion part has been divided into 2 parts: One to solve the engine obtaining the most important parameters and the other one to analyze the rotary elements (compressor and turbine).

## Main solution
Provides a solution given an engine's main parameters. To find a suitable engine the mentioned book should be the better source of information. Additional info can be taken from Internet. The book is an aircraft engine database (available at BCT):

Élodie Roux. Turbofan and turbojet engines: database handbook. Eng. Blagnac:
Élodie Roux, 2007. ISBN: 9782952938013

Both codes for analysis at cruising altitudes and sea level are available. 

## Rotary elements

Scripts for High Pressure Compressor (HPC) and High Pressure Turbine (HPT) analysis are provided. 
Regarding the compressor polytropic efficiencies as well as temperature increases are calculated. Then rotor blade angles are calculated. The mean radius is calculated for the turbine and the difference between the inlet and outlet Mach in different stages. 

