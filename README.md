# README #

### AFiNeS: Active Filament Network Simulation ###

##### as detailed in : #####
#### A versatile framework for simulating the dynamics mechanical structure of cytoskeletal networks ####

### Authors / Contributors ###

* Simon Freedman (University of Chicago) 
* Shiladitya Banerjee (University College London) 
* Glen Hocky (University of Chicago)
* Aaron Dinner (University of Chicago)

#### created at the University of Chicago ####

### System Requirements ###
Minimally, this system requires gcc+11 and boost which you can load on midway via the commands
```
    > module load gcc
    > module load boost
```

On a Mac, you can install Boost via MacPorts, with the command
```
sudo port install boost
```

### QUICKSTART GUIDE ###

* If you don't already have a bin directory, create one with:
    ```
    > mkdir bin
    ```    
* If you don't already have an executable, run the command: 
    ```
    > make [clean] [tar] network 
    ```
    * [clean] will delete the old executable
    * [tar] will generate the file tars/amxbd.tar.gz
    * IF this doesn't work, then there's probably a dependency or linker issue. Try each of the following solutions in the order prescribed, and attempt to compile in between. 
        1. Make sure you have Boost installed 
        2. Find the folder with the "*boost*.dylib" or "*boost*.a" folders; when I installed Boost using MacPorts, it was `/opt/local/lib/` . Run the command:

                 export BOOST_ROOT=<my boost folder>

              With the Macports installation, <my boost folder>=/opt/local/lib . 

        3. Within BOOST_ROOT, identify if the library folders have a suffix, such as "-mt" or "-d" (e.g., the program-options library on my Mac is named "libboost_program_options-mt.dylib", instead of "libbbost_program_options.dylib"). If so, run the command:
            
                export BOOST_SUFFIX=<my boost suffix>
            
             With the Macports installation, <my boost suffix>=-mt

        4. Find the folder with the boost/*.h files; with MacPorts installation, it was `/opt/local/include/`. Add `-I <myincludefolder>` to the line that begins `INC :=` in the makefile.

* You should now have an executable file called bin/afines. NOTE: you only need to recreate this file if you edit the source
  code.

* Create an output directory for your simulation (e.g., "out") 

```
> mkdir out
```

* Run your simulation in the specified output output directory, e.g., 
    ``` 
    > bin/afines --dir out
    ```

* See below for other simulation configuration options that you can set from the command line or from a configuration
  file

* Once your simulation has completed, the following files will have been generated:
 * out/txt_stack/actins.txt //the trajectories of every actin bead
 * out/txt_stack/links.txt //the trajectories of every link 
 * out/txt_stack/amotors.txt //the trajectories of all active motors (e.g., myosin) at every time step
 * out/txt_stack/pmotors.txt //the trajectories of all passive motors (e.g., crosslinkers) at every time step
 * out/data/thermo.txt //the energies of actin filaments
 * out/data/output.txt //some metadata about the simulation

All files are tab delimited 

* txt_stack/actins.txt has the format

    * x y r idx
        * (x, y)  = bead position, 
        * r  = bead radius 
        * idx = index of filament that the bead is on

* txt_stack/links.txt has the format
     * x y dx dy idx
         * (x, y) = end position of link closer to the barbed end of filament 
         * (x + dx, y + dy) = end position of link farther from barbed end 
         * idx = index of filament the link is on

* txt_stack/amotors.txt and txt_stack/pmotors.txt have the format
    * x y dx dy fidx0 fidx1 lidx0 lidx1
        * (x, y) = position of head 0 of motor
        * (x + dx, y + dy) = position of head 1 of motor
        * fidx0 = index of filament that head 0 is attached to (-1 if not attached)
        * fidx1 = index of filament that head 1 is attached to (-1 if not attached)
        * lidx0 = index of link that head 0 is attached to (-1 if fidx0 = -1)
        * lidx1 = index of link that head 1 is attached to (-1 if fidx1 = -1)

* data/filament_e.txt is the energetics of each filament, and has the format 
    * KE PE TE idx
        * KE = total v^2 of filament 
        * PE = total potential energy of filament
        * TE = KE + PE
        * idx = filament index
The time associated with the positions/energies is on it's own line before 
each list of positions within the file. Thus the structure of actins.txt is:
```
    t = t1
    x1, y1, r1, idx1
    .
    .
    .
    xn, yn, rn, idxn
    t=t2
    x1, y1, r1, idx1,
    .
    .
    .
    t=tn
    .
    .
    .
    xn, yn, rn, idxn
```

* data/pe.txt is the total potential energy of all particles at a given time step and has the format
    * U(filament_stretch)  U(filament_bend) U(xlink_stretch) U(motor_stretch) where each U is total at that timestep
    * time isn't delineated in these files; rather line 1 is t=t1, line 2, is t=t2, etc. 

* data/config_full.cfg is the full set of configuration options used for the simulation. Thus if a simulation did not
  complete, you can restart with 
```
./bin/afines -c data/config_full.cfg --restart true
```

### Configurable settings ###

Currently the following options for a simulation can be set upon execution, either from the command line, or within a
configuration file:

For example, to run a 500 second of simulation of 10 rigid actin filaments, an active motor density of 0.5 and a crosslinker density
of 0.05 you would enter the command:
    ```
    > ./bin/afines --tf 500 --npolymer 10 --a_motor_density 0.5 --p_motor_density 0.05
    ```
(this would write to the default output directory)

## Simulation Parameters ##

All the available options for a binary can be displayed using the command
    ```
    > ./bin/afines -h
    ```
These are reproduced below, but the above command should be more accurate.

|variable name              |type   |default value  |units  |description                                                                                       |
|---------------------------|:-----:|:-------------:|:-----:|--------------------------------------------------------------------------------------------------|
|**ENVIRONMENT**            |       |               |       |                                                                                                  |
|bnd_cnd,bc                 |string |"PERIODIC"     |       |boundary conditions                                                                               |
|xrange                     |double |10             |um     |size of cell in horizontal direction                                                              |
|yrange                     |double |10             |um     |size of cell in vertical direction                                                                |
|dt                         |double |0.0001         |s      |length of individual timestep                                                                     |
|tinit                      |double |0              |s      |time that recording of simulation starts                                                          |
|tfinal                     |double |0.01           |s      |length of simulation                                                                              |
|nframes                    |int    |1000           |       |number of frames of actin/link/motor positions printed to file (equally distant in time)          |
|nmsgs                      |int    |10000          |       |number of timesteps between printing simulation progress to stdout                                |
|viscosity                  |double |0.001          |mg/um*s|Dynamic viscosity                                                                                 |
|temperature,temp           |double |0.004          |pN*um  |Temp in energy units                                                                              |
|dir                        |string |"."            |       |directory for output files                                                                        |
|myseed                     |int    |time(NULL)     |       |seed of random number generator                                                                   |
|restart                    |bool   |false          |       |whether to restart the simulation from the last recorded time step                                |
|restart_time               |double |-1             |s      |time to restart the simulation from                                                               |
|restart_strain             |double |0              |um     |starting strain for restarting simulation                                                         |
|steven_continuation_flag   |bool   |false          |       |whether to start with the strain from the last recorded time step                                 |
|continuation_fr            |int    |0              |       |last saved from from the simulation to be continued                                               |
|grid_factor                |double |2              |um^(-1)|number of grid boxes per micron                                                                   |
|quad_off_flag              |bool   |false          |       |flag to turn off neighbor list updating                                                           |
|quad_update_period         |int    |1              |       |number of timesteps between actin/link/motor position updates to update quadrants                 |
|circle_flag                |bool   |false          |       |flag to add a circular wall                                                                       |
|circle_radius              |double |INFINITY       |um     |radius of circular wall                                                                           |
|circle_spring_constant     |double |0              |pN/um  |spring constant of circular wall                                                                  |
|**ACTIN**                  |       |               |       |                                                                                                  |
|actin_in                   |string |""             |       |input actin positions file                                                                        |
|npolymer                   |int    |3              |       |number of polymers in the network                                                                 |
|nmonomer                   |int    |11             |       |number of beads per filament                                                                      |
|nmonomer_extra             |int    |0              |       |maximum number of monomers per filament                                                           |
|extra_bead_prob            |double |0.5            |       |probability of adding an extra bead from nmonomer...nmonomer_extra                                |
|actin_length               |double |0.5            |um     |Length of a single actin monomer                                                                  |
|actin_pos_str              |string |               |       |Starting positions of actin polymers, commas delimit coordinates; semicolons delimit positions    |
|link_length                |double |0              |       |Length of links connecting monomers                                                               |
|polymer_bending_modulus    |double |0.068          |pN*um^2|Bending modulus of a filament                                                                     |
|fracture_force             |double |1000000        |pN     |filament breaking poiafines                                                                       |
|link_stretching_stiffness  |double |1              |pN/um  |stiffness of link                                                                                 |
|rmax                       |double |0.25           |um     |cutoff distance for interactions between actin beads and filaments                                |
|kexv                       |double |1.0            |pN/um  |parameter of exv force calculation                                                                |
|kgrow                      |double |0              |s^(-1) |rate of filament growth                                                                           |
|lgrow                      |double |0              |um     |additional length of filament upon growth                                                         |
|l0min                      |double |0              |um     |minimum length a link can shrink to before disappearing                                           |
|l0max                      |double |2              |um     |maximum length a link can grow to before breaking into two links                                  |
|nlink_max                  |int    |25             |       |maximum number of links allowed on filament                                                       |
|**MOTORS**                 |       |               |       |                                                                                                  |
|a_motor_in                 |string |""             |       |input motor positions file                                                                        |
|a_motor_density            |double |0.05           |um^(-2)|density of active motors                                                                          |
|a_motor_pos_str            |string |""             |       |Starting positions of motors, commas delimit coordinates; semicolons delimit positions            |
|motor_intersect_flag       |boolean|false          |       |if true, then motors are placed at filament intersections                                         |
|a_linkage_prob             |double |1              |       |probability that filaments are linked by a motor if motor_intersect_flag = true                   |
|a_m_kon                    |double |1              |s^(-1) |active motor on rate                                                                              |
|a_m_koff                   |double |0.1            |s^(-1) |active motor off rate                                                                             |
|a_m_kend                   |double |0.1            |s^(-1) |active motor off rate at filament end                                                             |
|a_m_cut                    |double |0.063          |um     |cutoff distance for binding                                                                       |
|a_m_kon2                   |double |-1             |s^(-1) |active motor on rate when both heads bound                                                        |
|a_m_koff2                  |double |-1             |s^(-1) |active motor off rate when both heads bound                                                       |
|a_m_kend2                  |double |-1             |s^(-1) |active motor off rate at filament end when both heads bound                                       |
|a_motor_length             |double |0.4            |um     |length of motor                                                                                   |
|a_motor_stiffness          |double |1              |pN/um  |active motor spring stiffness                                                                     |
|a_m_bend                   |double |0.04           |pN     |bending force constant for active motors                                                          |
|a_m_ang                    |double |pi/2           |       |equilibrium angle of active motor-filament system                                                 |
|a_m_align                  |int    |0              |       |active motor alignment type                                                                       |
|a_m_kalign                 |double |0              |pN*um  |active motor alignment penalty                                                                    |
|a_motor_v                  |double |1              |um/s   |velocity along filaments towards barbed end when attached                                         |
|a_m_stall                  |double |0.5            |pN     |stall force of motors                                                                             |
|dead_head_flag             |boolean|false          |       |if true, then head [dead_head] of all motors remains stationary throughout sim                    |
|dead_head                  |int    |0              |       |can be 0 or 1; head that remains stationary if dead_head_flag=true                                |
|**CROSSLINKS**             |       |               |       |                                                                                                  |
|p_motor_in                 |string |""             |       |input crosslinker positions file                                                                  |
|p_motor_density            |double |0.05           |um^(-2)|number of passive motors                                                                          |
|p_motor_pos_str            |string |""             |       |Starting positions of xlinks, commas delimit coordinates; semicolons delimit positions            |
|link_intersect_flag        |boolean|false          |       |if true, then crosslinks are placed at filament intersections                                     |
|p_linkage_prob             |double |1              |       |probability that filaments are crosslinked if link_intersect_flag = true                          |
|p_m_kon                    |double |1              |s^(-1) |passive motor on rate                                                                             |
|p_m_koff                   |double |0.1            |s^(-1) |passive motor off rate                                                                            |
|p_m_kend                   |double |0.1            |s^(-1) |passive motor off rate at filament end                                                            |
|p_m_cut                    |double |0.063          |um     |cutoff distance for binding (um)                                                                  |
|p_m_kon2                   |double |-1             |s^(-1) |passive motor on rate when both heads bound                                                       |
|p_m_koff2                  |double |-1             |s^(-1) |passive motor off rate when both heads bound                                                      |
|p_m_kend2                  |double |-1             |s^(-1) |passive motor off rate at filament end when both heads bound                                      |
|p_motor_length             |double |0.150          |um     |length of xlink                                                                                   |
|p_motor_stiffness          |double |1              |pN/um  |xlink spring stiffness                                                                            |
|p_m_bend                   |double |0.04           |pN     |bending force constant of passive motors                                                          |
|p_m_ang                    |double |pi/2           |       |equilibrium angle of passive motor-filament system                                                |
|p_m_align                  |int    |0              |       |passive motor alignment type                                                                      |
|p_m_kalign                 |double |0              |pN*um  |passive motor alignment penality                                                                  |
|p_motor_v                  |double |0              |um/s   |passive motor velocity                                                                            |
|p_m_stall                  |double |0              |pN     |stall force                                                                                       |
|p_dead_head_flag           |boolean|false          |       |if true, then head [p_dead_head] of all xlinks remains stationary throughout sim                  |
|p_dead_head                |int    |0              |       |can be 0 or 1; head that remains stationary if p_dead_head_flag=true                              |
|static_cl_flag             |boolean|false          |       |should be set to true if xlinks start off attached to filaments and never detach                  |
|**SHEAR**                  |       |               |       |                                                                                                  |
|n_bw_shear                 |int    |10^9           |s      |number of timesteps between subsequent differential strains                                       |
|d_strain_freq              |double |1              |Hz     |frequency of differential oscillatory strain                                                      |
|d_strain_pct               |double |0              |s      |differential strain (e.g., 0.5 means a strain of 0.5*xrange)                                      |
|time_of_dstrain            |double |10000          |s      |time when differential strain begins                                                              |
|time_of_dstrain2           |double |10000          |s      |time when second differential strain begins                                                       |
|diff_strain_flag           |boolean|false          |       |flag to use if differential strain should be linear (in one direction)                            |
|osc_strain_flag            |boolean|false          |       |flag to use if differential strain should be oscillatory (like Gardel, Science 2004)              |
|stress_flag                |boolean|false          |       |flag to turn on constant stress                                                                   |
|stress                     |double |0              |pN/um  |value of constant stress                                                                          |
|stress_rate                |double |0              |s/mg   |decay rate to specified value of stress                                                           |
|stress2                    |double |0              |pN/um  |second value of constant stress                                                                   |
|stress_rate2               |double |0              |s/mg   |second decay rate to the specified value of stress                                                |
|shear_motor_flag           |boolean|false          |       |flag to turn on shearing for motors                                                               |

### Configuration file Example ###
Below is an example of a configuration file named example.cfg. 
To run a simulation using this configuration, enter the command
     ```   
    >./bin/afines -c example.cfg
    ```
#### example.cfg ####
```
npolymer=500
nmonomer=1
dt=0.001
nframes=2000
tfinal=100
actin_length=0.5
a_motor_density=0.05
link_intersect_flag=true
actin_pos_str=0,0,0:1,2,3.141
```

### Contribution guidelines ###

* Email: dinner@uchicago.edu

### Who do I talk to? ###

* Simon Freedman (simonfreedman@uchicago.edu)
* Aaron Dinner (dinner@uchicago.edu)
