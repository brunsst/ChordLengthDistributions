# ChordLengthDistributions
Allows calculating various chord length distributions (CLDs) in 3D  in C++ for labeled 8bit tif-image sequences.

The program should build on any Linux distribution by executing the script *build_chordlengthdistributions.sh*.
<br>

The newly create binary needs to be executed with at least one mandatory arguments (*-i* followed by a **full absolute path** to the data. There is now safety check for faulty inputs.)
<br>

Default settings execute the code in **random** mode assuming a label of **255 for material phase, 0 for void phase and 128 for interface voxels**.
It is recommended to at least also increase the amount of sampled locations via **-n** to an amount that provides statistics with the desired precision (typically 10⁵ to 10⁶). The default of 10000 is for debugging purposes only.
<br>

Default settings create/append a *logfile.txt* in the *root* directory, i.e., one level above the input directory. The logfile tracks the statistical moments of the calculated CLDs.
Two additional files are created in the input directory. These are *chord_densityestimate_material.csv* and *chord_densityestimate_void.csv*


| argument | type | explanation |
|----------|------|-------------|
| -i       | string | full path to input directory which needs to contain an 8bit tif-image sequence. |
| -o       | string | full path to desired output directory. |
| -color   | uint8  | grayscale value of the label to be evaluated. Default is 255. |
| -r_min   | int    | Default is 1. Allows discarding incscribed spheres with small radii. |
| -vxl     | float  | Default 1.0. Allows setting the voxel size for scaling the output to physical length scales. |
| -n_cpu   | int    | Default is 128. Allows limiting the amount of threads for not clogging up all computational resources. |


