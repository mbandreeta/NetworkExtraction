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

def load_data(return_tot=False):
    import re;
    skel = np.loadtxt("_skel.txt",dtype=np.float32);
    nodes = np.loadtxt("_nodes_skel.txt",dtype=np.float32);
    file = open("_validate.txt", "r") 
    aux = (file.readline())
    filename = (file.readline())
    unit = float(file.readline())
    file.close();
    image_shape = re.findall(r"[-+]?\d*\.\d+|\d+", aux)
    G = nx.Graph();
    G.graph['unit']=unit;
    G.graph['shape'] = [int(image_shape[0]), int(image_shape[1]), int(image_shape[2])]
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
    

def load_data_vol():
    import re;
    skel = np.loadtxt("D:/Rochas/_skel.txt",dtype=np.float32);
    nodes = np.loadtxt("D:/Rochas/_nodes_skel.txt",dtype=np.float32);
    file = open("D:/Rochas/_validate.txt", "r") 
    aux = (file.readline())
    filename = (file.readline())
    unit = float(file.readline())
    file.close();
    image_shape = re.findall(r"[-+]?\d*\.\d+|\d+", aux)
    G = nx.Graph();
    G.graph['unit']=unit;
    G.graph['shape'] = [int(image_shape[0]), int(image_shape[1]), int(image_shape[2])]
    nlist=list();
    ndict={};
    border = list();
    nvol = local_volume_surface();
    
    for i in range(nodes.shape[0]):
        ndict[nodes[i][0]]= [nodes[i][4],nodes[i][7]];
        if nodes[i][6]==0:
            border.append(nodes[i][0]);
        border_type = getBorderType(nodes[i][6]);
        volume = nodes[i][7]
        if nodes[i][0]  in nvol.keys():
            volume = nodes[i][7]+nvol[nodes[i][0]]
        nlist.append((nodes[i][0],{
                            'metadata': {'node_coordinates': {'x': nodes[i][1], 'y': nodes[i][2], 'z': nodes[i][3]},
                                            'node_border': {'border_type': border_type},
                                            'node_squared_radius': nodes[i][4]**2,
                                            'node_volume': volume,
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
        Gskel.add_edge(skel[i][0],skel[i][1],diameter=_diameter, distance = _distance);
    Gc = max(nx.connected_component_subgraphs(Gskel), key=len);
    filename = filename.rstrip('\n');
    filename = filename.split('.');
    return Gc,filename[0],unit;    
    
def local_volume_surface():
    import pandas as pd
    edgesdf = pd.read_csv("D:/Rochas/_edges_graph.txt", sep=" ", header=None)
    print("read edges");
    nodesdf = pd.read_csv("D:/Rochas/_nodes.txt", sep=" ", header=None)
    print("read nodes");
    skel = np.loadtxt("D:/Rochas/_skel.txt",dtype=np.float32);
    Gskel = nx.Graph()
    Gtot = nx.Graph()
    for i in range(skel.shape[0]):
        Gskel.add_edge(skel[i][0],skel[i][1]);
    Gc = max(nx.connected_component_subgraphs(Gskel), key=len);
    nlist=list();
    volume={};
    surface={};
    border = list();
    nodes = nodesdf.to_numpy()
    edges = edgesdf.to_numpy()
    visited = {}
    for i in range(nodes.shape[0]):
        volume[nodes[i][0]]= nodes[i][7]
        surface[nodes[i][0]]= nodes[i][8]
        visited[nodes[i][0]] = False;
    print('finished nodes')

    for i in (range(edges.shape[0])):
        Gtot.add_edge(edges[i][0],edges[i][1]);

    nvol = {}
    for node in Gc.nodes():
        nvol[node] = 0
        visited[node] = True;
        for n in nx.neighbors(Gtot,node):
            if(visited[n]==False):
                nvol[node]+=volume[n];
                visited[n] =  True;
    return nvol;
            
            
            
    
    
def load_full_data(filepath):
    import re;
    import pandas as pd
    skeldf = pd.read_csv("_all_edges_graph.txt", sep=" ", header=None)
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
                                            'node_isPore': False
                                            }
                                            }));
    print('finished nodes')
    Gskel = nx.Graph();
    Gskel.graph['unit']=unit;
    Gskel.graph['shape'] = [int(image_shape[0]), int(image_shape[1]), int(image_shape[2])]
    Gskel.add_nodes_from(nlist); 
    nodes = Gskel.node_dict_factory(Gskel.nodes(data=True));
    skel = skeldf.to_numpy();
    for i in tqdm.tqdm(range(skel.shape[0])):
        _diameter = skel[i][2];
        _distance = distanceBetweenNodes(nodes,int(skel[i][0]),int(skel[i][1]));
        Gskel.add_edge(skel[i][0],skel[i][1],diameter=_diameter, distance = _distance);
    Gc = max(nx.connected_component_subgraphs(Gskel), key=len);

    return Gc,filename[0],unit;
    
def load_data_mb(prefix):
    import re;
    link1 = prefix+"_link1.dat"
    link2 = prefix+"_link2.dat"
    node1 = prefix+"_node1.dat"
    node2 = prefix+"_node2.dat"
    i=0;
    nodeDict = {}
    with open(node1, "r") as file:
        line = file.readline().rstrip('\n');
        values=[l for l in line.split(" ") if l!=''];
        Npores = int(values[0]);
        image_shape = [float(values[1]),float(values[2]),float(values[3])];
        print(Npores)
        for i in range(1,Npores+1):
            line = file.readline().rstrip('\n');
            values=[l for l in line.split(" ") if l!=''];
            nodeDict[int(values[0])]={
                                'metadata': {'node_coordinates': {'x': float(values[1]), 'y': float(values[2]), 'z': float(values[3])},
                                                'node_border': {'border_type': None},
                                                'node_squared_radius': 0,
                                                'node_volume': 0,
                                                'node_shape_factor': 0,
                                                'node_clay_volume': 0
                                                }
                                    }    
    with open(node2, "r") as file:
        line = file.readline().rstrip('\n');
        while line:
            values=[l for l in line.split(" ") if l!=''];            
            x,y,z = getXYZ(nodeDict[int(values[0])])
            nodeDict[int(values[0])]['metadata']['node_volume']= float(values[1]);
            nodeDict[int(values[0])]['metadata']['node_squared_radius']= float(values[2])**2;
            nodeDict[int(values[0])]['metadata']['node_border']['border_type']= checkBorder(image_shape,float(values[2]),x,y,z);
            nodeDict[int(values[0])]['metadata']['node_shape_factor']= float(values[3]);
            nodeDict[int(values[0])]['metadata']['node_clay_volume']= float(values[4]);    
            line = file.readline().rstrip('\n');

    nlist=[(k,v) for k,v in nodeDict.items()];
    Gskel = nx.Graph();
    Gskel.graph['unit']=1e6;
    Gskel.graph['shape'] = [(image_shape[0]), (image_shape[1]), (image_shape[2])]
    Gskel.add_nodes_from(nlist); 
    nodes = Gskel.node_dict_factory(Gskel.nodes(data=True));
    tDict={}
    Nthroats=0;
    i=0;
    with open(link1, "r") as file:
        line = file.readline().rstrip('\n');
        while line:
            values=[l for l in line.split(" ") if l!=''];
            if i ==0:
                Nthroats = float(values[0]);
            else:
                if(int(values[1])>0 and int(values[2])>0 ):
                    tDict[int(values[0])]= { 'p1':int(values[1]),
                                                "p2":int(values[2]),
                                                "radius":float(values[3]),
                                                "shape_factor":float(values[4]),
                                                "total_length":float(values[5]),
                                                "length_pore1": 0,
                                                "length_pore2": 0,
                                                "length_throat": 0,
                                                "throat_volume": 0,
                                                "throat_clay_volume": 0
                                                }    
            i+=1;
            line = file.readline().rstrip('\n');    
    with open(link2, "r") as file:
        line = file.readline().rstrip('\n');
        while line:
            values=[l for l in line.split(" ") if l!=''];
            if int(values[0]) in tDict.keys():
                tDict[int(values[0])]["length_pore1"] = float(values[3])
                tDict[int(values[0])]["length_pore2"] = float(values[4])
                tDict[int(values[0])]["length_throat"] = float(values[5])
                tDict[int(values[0])]["throat_volume"] = float(values[6])
                tDict[int(values[0])]["throat_clay_volume"] = float(values[7])
            line = file.readline().rstrip('\n');
            
    for k,v in tDict.items():
        Gskel.add_edge(v["p1"],v["p2"],diameter=v["radius"]*2, distance = v["length_throat"],
        length_pore1=v["length_pore1"], 
        length_pore2=v["length_pore2"], 
        length_throat=v["length_throat"], 
        throat_volume=v["throat_volume"], 
        throat_clay_volume=v["throat_clay_volume"], 
        
        );
    return Gskel;    