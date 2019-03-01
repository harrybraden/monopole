# Monopole Visualisation

This repository contains code to generate and visualise charge 2 monopoles.

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

- [1] https://peterbraden.co.uk/monopole/browse.html
- [2] Lorensen, W. E.; Cline, Harvey E. (1987). "Marching cubes: A high resolution 3d surface construction algorithm". ACM Computer Graphics. 21 (4): 163–169. doi:10.1145/37402.37422


## Contents

- `README.md`: This documentation.
- `Makefile`: to build the visualisations. See below.
- `visualise/`: The interactive visualiser.
- `generate-image-data.py`: Generates the image data for the interactive web visualiser
- `data.py`: Load and manipulate the data
- `contours-image.py`: Generate the 'Tomogram' image
- `generatemesh.py`: Generate a 3d obj file from the data

## Running the visualisations

All the visualisations can be run using the Makefile.

To make all:
```sh
make
```

To make a 3d model of a particular set of k and threshold parameters:
```
make models/k=0.55.th=0.57.obj
```

The interactive viewer is a web page in the `visualise` folder and can be run
locally with:
```
make run-webapp
```
Then navigate to [http://localhost:8080/browse.html](http://localhost:8080/browse.html)

