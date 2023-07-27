# ChordLengthDistributions
Allows calculating various chord length distributions (CLDs) in 3D  in C++ for binary 8bit tif-image sequences.

The program should build on any Linux distribution by executing the script *build_chordlengthdistributions.sh*.
<br>

The newly create binary needs to be executed with at least one mandatory arguments (*-i* followed by a **full absolute path** to the data. There is now safety check for faulty inputs.)
<br>

Default settings execute the code in **random** mode assuming a label of **255 for material phase, 0 for void phase and 128 for interface voxels**.
It is recommended to at least also increase the amount of sampled locations via **-n** to an amount that provides statistics with the desired precision (typically 10⁵ to 10⁶). The default of 10000 is for debugging purposes only.
<br>

Default settings create/append a *logfile.txt* in the *root* directory, i.e., one level above the input directory. The logfile tracks the statistical moments of the calculated CLDs.
Two additional files are created in the input directory. These are *chord_densityestimate_material.csv* and *chord_densityestimate_void.csv* which contain the CLDs for the two phases in the input image with chordlength in the first and frequency in the second column. Distributions are provided from a kernel density estimation using a Epanechnikov kernel with IQR bandwidth estimator and a stepsize of 1.

### Available Modes

The mode of calculation can be set with the -mode argument followed by a string defining the mode:

- *random* samples chords at a random location in a random correction until at least the amount of chords set by -n is reached in both phases.
- *average* samples all directions at random locations and records the mean chordlength which requires substantially more runtime.
- *minchord* is similar to *average* but records the smallest chord which may be usefull for describing necks.
- *directionality* is a quick port to another CLD program. The output is a csv-file (*directionality.csv*) in the root directory that contains an individual CLD for all the directions sampled. 

### Input Data

Input data should be binary and a 8bit tif sequence. However interface voxels may be defined with another label color. This is useful for segmentations acquired via watershed segmentation where the interface is often undecided between foreground and background. Chords are calculated simply by expanding a chord along the trajectory around the sample location via integer sized steps. The expansions terminates when the other phase is reached. If an interface voxel is reached the expansions also terminates but half a step is still added. This can be significant in small/fine structures.

| argument | type | explanation |
|----------|------|-------------|
| -i       | string | full path to input directory which needs to contain an 8bit tif-image sequence. Repeating the argument queues multiple directories. |
| -n       | int    | Default is 10000. Defines the amount of random locations sampled. |
| -mode    | string | Sets the mode of calculation. |
| -vxl     | float  | Default 1.0. Allows setting the voxel size for scaling the output to physical length scales. |
| -angle   | double | Default 15. Angular spacing in degrees between chord directions. Defines the granularity of chord directions sampled. |
| -cm   | uint8  | grayscale value of the material phase to be evaluated. Default is 255. |
| -cv   | uint8  | grayscale value of the void phase to be evaluated. Default is 0. |
| -ci   | uint8  | grayscale value of otional interface phase. Default is 128. |
| -n_cpu   | int    | Default is 128. Allows limiting the amount of threads for not clogging up all computational resources. |


