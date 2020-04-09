# NetworkExtraction
Compiled version of the Max Ball Medial Axis algorithm
This is a executable program for extracting a Medial Axis Pore Network from binary images of porous media.
 Usage: PoreNetworkMedial.exe Berea500.nii 2.25 
 The input parameters are a NIFITI type volume and the voxel resolution. The choice for NIFITI is to include Magnetic Resonance Imaging related images. 
Here I have also included python scripts to load and save NIFITI type images using nibabel

The output are .txt files:
_validate.txt 
_skel.txt
_nodes_skel.txt
_all_nodes.txt
_edges_graph.txt

The validate file contains:

500 500 500 // the shape of the volume x,y,z

Berea500.nii // name of the file

2.25 // resolution that was input

Total volume: 27094864 // total volume of the porous space in voxels

:: Total sphere volume:27094864 // the total volume by summing the volume assigned to each sphere in voxels

Total internal surface: 10974798 // estimated internal surface of the pore space in voxels

:: Total sphere surface: 10974798 // total surface by summing the surface voxels assigned to each sphere in voxels

Porosity %: 0.216759 // sample´s porosity (considering the shape in the first line)


_skel.txt  // contains the edges of the medial axis, the file is organized as:

pore1_id pore2_id thoat_diameter


_nodes_skel.txt // contains the node information of the spheres selected of the medial axis, the file is organized as:

nodeId x_index y_index z_index radius isBorder border_type  volume contact_surface


x_index,y_index,z_index : sphere's center position

radius: sphere's radius

isBorder: is sphere is near the borders of the volume, this is set to 1, 0 otherwise

border_type: integer value to identify the face that the sphere is connected to

	border_type = 1;//"inlet_x";
	
	border_type = 2;//"inlet_y";
	
	border_type = 3;//"inlet_z";
	
	border_type = 4;//"outlet_x";
  
	border_type = 5;//"outlet_y";
	
	border_type = 6;//"outlet_z";
	
volume: Sphere's assigned volume in voxels.Each voxel is assigned to one sphere only.

contact_surface: Sphere's contact area, an estimate of the voxels assigned to that sphere's volume that also belong to the pore space internal surface 


_all_nodes.txt // is the collection of all maximal spheres that completely fill the volume, the file is organized the same way as _nodes_skel.txt

_edges_graph.txt // the connections between spheres, this file is organized as _skel.txt




