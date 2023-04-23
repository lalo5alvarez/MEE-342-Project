# MEE-342-Project

**IGNORE ALL FILES BESIDES: 
 - 'input_shaft.m'
 - -countershaft.m'
 - output_shaft.m'

**Inputs: 

The inputs that the code accept are the total lengths of the shafts (all three), gear radius (two), fillet radii, surface finish, operating temperatures, and types of loading. 

 - Total Lengths of Shafts -> Up to the user to define! 
 - Gear Radius -> The first gear and its pinion are inputs, the code assings the second set!
 - Fillet Radii -> The code asks the user for the desired fillet radius within a constrained bound @ each geometry change.
 - Surface Finish -> User's input provides criteria for calculating Endurance Limit
 - Operating Temperature -> User's input provides criteria for calculating Endurance Limit
 - Type of Loading -> User's input provides criteria for calculating Endurance Limit
 - By altering the code, you can change material! (Sy, Sut)

**Note: The code shall be run IN ORDER!! 
1.) input_shaft.m
2.) countershaft.m
3.) output_shaft.m


**Outputs: 

 - Shear Diagrams
 - Moment Diagrams
 - Torque Diagrams
 - Optimized Diameters (Rounded to nearest .025 in.)
 - Locations of Diameter Steps

