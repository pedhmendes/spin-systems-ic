# Spins Systems

This repo contains all thing about my 2020 IC research. I worked with Professor Heitor C. M. Fernandes, from the Complex Fluids Group IF-UFRGS. I studied Spins Systems for a while, my goal was to parallelize my codes. I pause this project for a while, but I'm still studying how to achieve this. I will talk about the Ising Model and what I learned about it.

## Ising 2D

If you know Portuguese I highly recommend this [wiki](https://fiscomp.if.ufrgs.br/index.php/Ising_2D) about the Ising Model. In english there is the usual [wiki](https://en.wikipedia.org/wiki/Ising_model). 

To compile my Ising Model code is necessary that ```mc.h``` is in the same directory. I used the RNG avaiable in this file. To run you may provide and temperature as argument. Follow the examples to compile and run.

```gcc ising_model_serial.c -lm```
```./a.out TEMP```

where ```TEMP``` is the temperature of the system. In the end the program returns the execution time. I did this in order to compare the execution times between the serial code and the future parallelized codes.


## OpenMP

Eventually I will talk about parallelization, the tool chose was OpenMP. More information about it can be found in this [link](https://www.openmp.org/). 
