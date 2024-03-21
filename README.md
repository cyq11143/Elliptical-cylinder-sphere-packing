# Elliptical-cylinder-sphere-packing
## Introduction
This project provides the codes for Monte Carlo Simulation of ellpitical cylinder sphere packing. We have designed the **rand_function** module, the MC **compress** simulation module and the **output** module, along with a **main.cpp** for the overall execution. To run the simulation, simply enter the `g++ -std=c++0x rand_function.h rand_function.cpp compress.h compress.cpp output.h output.cpp main.cpp -o result` on  a Linux system. This command will generate an executable file named 'result'. Running `./result` will initiate the simulation.
## Relevant parameter setting
After initializing the program, you will need to set values for certain parameters in the interactive dialog box. Here, **r** signifies the ratio of the long diameter to the short diameter, while **a** corresponds to the short diameter. Additionally, **N_start** and **N_end** define the range within which the sphere number **N** will be iterated, and **N_loop** specifies the number of repetitions for each set of **(r, a, N)**.
## Other information
In the **main.cpp** file, we have defined numerous parameters, including the pressure **p**, the rate of temperature decay, and the initial height **L**. To analyze different systems, you will need to adjust these parameters accordingly.





 
