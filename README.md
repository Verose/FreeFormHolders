# FreeFromHolders
This project loads a 3D OBJ model and creates the geometry of a holder.  
Our implementation follows the lines of AutoConnect's free form holders expansion algorithm:  
https://koyama.xyz/project/AutoConnect/autoconnect.pdf

Instructions:
* A 3D model geometry is loaded
* Select a source point on the model
* The cut is calculated using geodesic distances from the source point, and is marked red
* This parameter is configurable using the 'Max Distance &' box
* Pressing '1' finalizes the cut and creates the holder
* Pressing '2' clears the mesh in order to choose a new source point
* Pressing '3' loads a new 3D model
* Pressing '4' saves the holder to and OBJ file


