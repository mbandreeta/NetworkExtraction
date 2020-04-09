import numpy as np;
import networkx as nx;


def getNiftiObject(filepath):
	import nibabel as nib;
	img = nib.load(filepath);
	header = img.header
	pixdim = header.get('pixdim');
	unit = pixdim[1]
	img_array = img.get_data();
	s = img_array.shape;
	if len(s)>3:
		img_array = img_array[:,:,:,0];
	return [img_array,unit];
	

	
def saveNiftiObject(data,unit,filepath):
	import nibabel as nib;
	print(unit);
	header_ = nib.Nifti1Header();
	header_['pixdim'] = np.array([ unit ,unit,unit,unit ,  1. ,  1. ,  1. ,  1. ], dtype=np.float32);
	nifti_img = nib.Nifti1Image(data, affine=None, header=header_);
	nib.save(nifti_img, filepath)



def saveGraph(G,filepath):
    import gzip
    import json
    from networkx.readwrite import json_graph
    data = json_graph.node_link_data(G)
    json_str = json.dumps(data, default=lambda x: float(x))

    
    json_bytes = json_str.encode('utf-8') 
    with gzip.GzipFile(filepath, 'w') as fout:   # 4. gzip
        fout.write(json_bytes)        

def loadGraph(filepath):
    import gzip
    import json
    from networkx.readwrite import json_graph
    with gzip.GzipFile(filepath, 'rb') as file:
        json_str = file.read().decode('utf-8').replace('\0', '')
    data = json.loads(json_str[json_str.find('{')::], encoding='utf-8')
    G = json_graph.node_link_graph(data)
    return G,data;
    
    
def getPrefix():
    from tkinter.filedialog import askopenfilename
    from tkinter import Tk
    import re
    Tk().withdraw()
    filepath = askopenfilename();
    dir_names = filepath.split("/");
    filetype = dir_names[-1].split("_")[-1];
    name = dir_names[-1].replace('_'+filetype,'');
    filepath = filepath.replace(dir_names[-1],name)
    return filepath;

def getPathPrefix():
    from tkinter.filedialog import askopenfilename
    from tkinter import Tk
    import re
    Tk().withdraw()
    filepath = askopenfilename();
    dir_names = filepath.split("/");
    prefix = dir_names[-1].split("_")[0];
    filepath = filepath.replace(dir_names[-1],'')
    return filepath,prefix;

def getBorderType(v):
    if v == -1:
        return "not_border"
    if v == 1:
        return "inlet_x"
    if v == 2:
        return "inlet_y"
    if v == 3:
        return "inlet_z"
    if v == 4:
        return "outlet_x"        
    if v == 5:
        return "outlet_y"
    if v == 6:
        return "outlet_z"

def distanceBetweenNodes(nodes,node01,node02):
    x01 = nodes[node01]['metadata']['node_coordinates']['x']*1.0;
    y01 = nodes[node01]['metadata']['node_coordinates']['y']*1.0;
    z01 = nodes[node01]['metadata']['node_coordinates']['z']*1.0;        
    x02 = nodes[node02]['metadata']['node_coordinates']['x']*1.0;
    y02 = nodes[node02]['metadata']['node_coordinates']['y']*1.0;
    z02 = nodes[node02]['metadata']['node_coordinates']['z']*1.0;
    return np.sqrt(((x01-x02)**2)+((y01-y02)**2)+((z01-z02)**2));

def load_data_skel(return_tot=False):
    import re;
    skel = np.loadtxt("_skel.txt",dtype=np.float32);
    nodes = np.loadtxt("_nodes_skel.txt",dtype=np.float32);
    file = open("_validate.txt", "r") 
    aux = (file.readline())
    filename = (file.readline())
    unit = float(file.readline())
    file.close();
    image_shape = re.findall(r"[-+]?\d*\.\d+|\d+", aux)
    nlist=list();
    ndict={};
    border = list();
    for i in range(nodes.shape[0]):
        ndict[nodes[i][0]]= [nodes[i][4],nodes[i][7]];
        if nodes[i][6]==0:
            border.append(nodes[i][0]);
        border_type = getBorderType(nodes[i][6]);
        nlist.append((nodes[i][0],{
                            'metadata': {'node_coordinates': {'x': nodes[i][1], 'y': nodes[i][2], 'z': nodes[i][3]},
                                            'node_border': {'border_type': border_type},
                                            'node_squared_radius': nodes[i][4]**2,
                                            'node_volume': nodes[i][7],
                                            'node_isPore': False
                                            }
                                            }));
    Gskel = nx.Graph();
    Gskel.graph['unit']=unit;
    Gskel.graph['shape'] = [int(image_shape[0]), int(image_shape[1]), int(image_shape[2])]
    Gskel.add_nodes_from(nlist); 
    nodes = Gskel.node_dict_factory(Gskel.nodes(data=True));
    
    for i in range(skel.shape[0]):
        _diameter = skel[i][2];
        _distance = distanceBetweenNodes(nodes,int(skel[i][0]),int(skel[i][1]));
        w = (_diameter/2.0)**4/_distance;
        Gskel.add_edge(skel[i][0],skel[i][1],diameter=_diameter, distance = _distance, weight = 1.0/w);
    c = max(nx.nx.connected_components(Gskel), key=len);
    Gc = Gskel.subgraph(c).copy();
    filename = filename.rstrip('\n');
    filename = filename.split('.');
    if return_tot:
        return Gc,Gskel,filename[0],unit; 
    else:
        return Gc,filename[0],unit;    
    


    
def load_full_data(filepath):
    import re;
    import pandas as pd
    edge_data = pd.read_csv("_all_edges_graph.txt", sep=" ", header=None)
    print("read edges");
    nodesdf = pd.read_csv("_all_nodes.txt", sep=" ", header=None)
    print("read nodes");
    file = open("D:/Rochas/_validate.txt", "r") 
    aux = (file.readline())
    filename = (file.readline())
    unit = float(file.readline())
    file.close();
    filename = filename.rstrip('\n');
    filename = filename.split('.');
    image_shape = re.findall(r"[-+]?\d*\.\d+|\d+", aux)
    nlist=list();
    ndict={};
    border = list();
    nodes = nodesdf.to_numpy()
    for i in range(nodes.shape[0]):
        ndict[nodes[i][0]]= [nodes[i][4],nodes[i][7]];
        if nodes[i][6]==0:
            border.append(nodes[i][0]);
        border_type = getBorderType(nodes[i][6]);
        nlist.append((nodes[i][0],{
                            'metadata': {'node_coordinates': {'x': nodes[i][1], 'y': nodes[i][2], 'z': nodes[i][3]},
                                            'node_border': {'border_type': border_type},
                                            'node_squared_radius': nodes[i][4]**2,
                                            'node_volume': nodes[i][7],
                                            'node_surface': nodes[i][8],
                                            'node_isPore': False
                                            }
                                            }));
    print('finished nodes')
    G = nx.Graph();
    G.graph['unit']=unit;
    G.graph['shape'] = [int(image_shape[0]), int(image_shape[1]), int(image_shape[2])]
    G.add_nodes_from(nlist); 
    nodes = G.node_dict_factory(G.nodes(data=True));
    edge_data = edge_data.to_numpy();
    for i in tqdm.tqdm(range(edge_data.shape[0])):
        _diameter = edge_data[i][2];
        _distance = distanceBetweenNodes(nodes,int(edge_data[i][0]),int(edge_data[i][1]));
        G.add_edge(edge_data[i][0],edge_data[i][1],diameter=_diameter, distance = _distance);
    Gc = max(nx.connected_component_subgraphs(G), key=len);

    return Gc,filename[0],unit;
    
