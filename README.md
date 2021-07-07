# abaqus-knits

currently 2 objects in `beam_knit_obj.py` (one day that'll get renamed hopefully):
1) beam knit (notched + perf)
2) 3d knit (notched + perf)

each object is initialized w/ 2 named tuples that represent the (1) yarn properties and (2) the stitch size properties along with r the radius of the yarn. The radius is a little fake compared to real life measurments b/c abaqus is...picky (eg for stitch size 5 I used to have r=0.255 but for the notched 3d knit having an r>0.245 causes parts of it to not mesh...for some reason)

then each knit has some class functions to create the `.cae` file for a perf/notched knit and create an `.inp` file from that. There are then global functions to run the `.inp` files
