# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

BOOTSTRAP_MPI = 100
SIZE_COEFFICIENT = 9.75 # coefficient to swap between coordinates in numbers of pixels and micrometers
VOLUME_OFFSET = 1 # what does this mean ???
PERIPHERAL_FRACTION_THRESHOLD = 10
NUM_CONTOURS = 100 # number of isolines for periferal distance map
MAX_CELL_RADIUS = 300
#MAX_CELL_RADIUS = 20
DET_TOLERANCE = 0.00000001
PIXEL_PER_VOXEL = 16   # pixels per voxel width
VOXELS_PER_IMAGE = 512 # voxels per image width
RIPLEY_K_SIMULATION_NUMBER=20
STRIPE_NUM=3
SLICE_THICKNESS=0.3
Z_LINE_SPACING=20
PIXELS_IN_SLICE= SLICE_THICKNESS * SIZE_COEFFICIENT
VOLUME_COEFFICIENT= ((1 / SIZE_COEFFICIENT)**2) * 0.3 #conversion from pixel volume to micrometer volume