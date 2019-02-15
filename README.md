# Monopole Visualisation

This repository contains code to generate and visualise charge 2 monopoles.



## Contents

`Makefile`: to build the visualisations. See below.
`visualise/`: The interactive dashboard. 

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

