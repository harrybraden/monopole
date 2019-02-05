The output of the calculation is a 5 dimensional dataset, bucketed in values of
K. Because we can think of the energy density as a scalar value and the other three
dimensions as dimensions of space, we can therefore consider the data as three
dimensional volumetric data over time.

In order to initially visualise the data, an interactive viewer was created [1]
which allows the data to be dragged around. This viewer uses the energy density
as opacity, and hides all volumes below this value in order to look inside the
volume.

One obvious way to visualise the volumetric data is to define a threshold above
which to consider as solid, and use the Marching Cubes algorithm[2] to construct
a mesh of that threshold's contour.

These meshes can be visualised with the many mesh viewers, or even 3D printed.



[1] https://peterbraden.co.uk/monopole/browse.html
[2] Lorensen, W. E.; Cline, Harvey E. (1987). "Marching cubes: A high resolution 3d surface construction algorithm". ACM Computer Graphics. 21 (4): 163–169. doi:10.1145/37402.37422

