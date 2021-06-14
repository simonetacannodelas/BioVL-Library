
# Minerva
Welcome to Minerva - the BioVL library that collects some of the models inside FermProc software implemented and solved in Python.
![logo_minerva](https://user-images.githubusercontent.com/51289658/121901186-775d8100-cd26-11eb-8afd-37f77fe883ea.png)

## Mechanistics models for Fermentation

Along their benefits, mechanistic models:

· Summary of knowledge about the process.

· Possibility to try different scenarios.

## Implementation of a fermentation mechanistic model 

In order to create a template, we follow the next steps:
 1. Create a class.
 
    1.1 Initialize variables, constants, and parameters.
  
    1.2 Create a process matrix in which it is enclosed the kinetic reaction.
   
    1.3 Create a process matrix in which it is enclosed the kinetic reaction.
    
    1.4 Solution of the model through a solver such as the "odeint" module of scipy.integrate for ODE.

## Models inside this repository

In this repository, you can find different models for the aerobic or anaerobic growth of microorganisms.
The models are implemented as an object and they are based on determistic principles.
In this library you can find:

 > Herbert-Monod in aerobic conditions
 
 > Herbert-Monod in anaerobic conditions
 
 > Aerobic growth of Saccharomyces cerevisiae in glucose (with metabolic switch controlled by the Cabtree Effect/Overflow and glucose inhibition)



