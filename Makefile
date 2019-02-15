
all: installed models/k=0.55.th=0.57.obj models/k=0.01.th=0.57.obj contours.png
.PHONY: all


installed:
	pip install numpy
	pip install --upgrade PyMCubes
	python -m pip install numpy scipy matplotlib
	pip install pillow
	pip install -U scikit-image
	touch installed

models/%.obj:
	python generatemesh.py $(shell echo $(@) | sed -r 's/models\/k=(.*)\.th.*/\1/') $(shell echo $(@) | sed -r 's/models\/k=.*\.th=(.*)\.obj/\1/') > $(@)


contours.png: contours-image.py
	python contours-image.py

run-webapp:
	cd visualise && python -m SimpleHTTPServer 8080
.PHONY: run-webapp
