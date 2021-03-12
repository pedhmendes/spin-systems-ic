# Spins Systems

This repo contains all thing about my 2020 IC research. I worked with Professor Heitor C. M. Fernandes, from the Complex Fluids Group IF-UFRGS. I studied Spins Systems for a while, my goal was to parallelize my codes. I pause this project for a while, but I'm still studying how to achieve this. I will talk about the Ising Model and what I learned about it.

## Ising 2D

In Portuguese I highly recommend this [wiki](https://fiscomp.if.ufrgs.br/index.php/Ising_2D) about the Ising Model. In english there is the usual [wiki](https://en.wikipedia.org/wiki/Ising_model). 

To compile my Ising Model code is necessary that ```mc.h``` is in the same directory. I used the RNG avaiable in this file. To run you may provide and temperature as argument. Follow the examples to compile and run.

```gcc ising_model_serial.c -lm```

```./a.out TEMP```

where ```TEMP``` is the temperature of the system. In the end the program returns the execution time. I did this in order to compare the execution times between the serial code and the future parallelized codes. There are some plots in a folder where I compared different matrixes sizes with different temperatures. There is also temporal series and correlation times for L = 16.

I did a presentation video for [SIC 2020](https://www.ufrgs.br/propesq1/sic2020/). You can find the youtube link in this [link](https://www.youtube.com/watch?v=nI9L4SJyBcA). The presentation file and the latex project can be found as well. The video and the presentation file are in portuguese.


## OpenMP

Eventually I will talk about parallelization, the tool chose was OpenMP. More information about it can be found in this [link](https://www.openmp.org/). 

I'm currently doing a course about OpenMp in Coursera, Fundamentals of Parallelism on Intel Architecture. More information about this course can be found [here](https://www.coursera.org/learn/parallelism-ia).
