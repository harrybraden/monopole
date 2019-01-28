import generatemesh
from PIL import Image



def draw_contours():
    print "# Load data" 
    volume = generatemesh.load_data(0.55) 
    print "# Trace contours"

if __name__ == "__main__":
    draw_contours()
