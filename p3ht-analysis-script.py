import gsd.hoomd
import numpy as np
from morphct.chromophores import amber_dict, get_chromo_ids_smiles
from morphct.system import System
from morphct.mobility_kmc import snap_molecule_indices

gsdfile = "crystalline_AA.gsd"
system = System(gsdfile, "output_p3ht", frame=0, scale=1, conversion_dict=amber_dict)
front_carbs = np.array(list(range(0, 11)))
end_carbs = np.array(list(range(154, 165)))
middle_list_carbs = []
middle_chromo_carbs = np.array(list(range(11,22)))
for x in range(13):
    middle_list_carbs.append(middle_chromo_carbs)
    middle_chromo_carbs = middle_chromo_carbs + 11
middle_list_carbs.append(front_carbs)
middle_list_carbs.append(end_carbs)
front_hydros = np.array(list(range(165000, 165015)))
end_hydros = np.array(list(range(165197, 165212)))
middle_list_hydros = []
middle_chromo_hydros = np.array(list(range(165015,165029)))
for x in range(13):
    middle_list_hydros.append(middle_chromo_hydros)
    middle_chromo_hydros = middle_chromo_hydros + 14
middle_list_hydros.append(front_hydros)
middle_list_hydros.append(end_hydros)
master_carbs_list = []
sublist = middle_list_carbs
for i in range(1000):    
    for j in sublist:
        master_carbs_list.append(j)
    sublist = [x + 165 for x in sublist]
    
master_hydros_list = []
sublist = middle_list_hydros
for i in range(1000):    
    for j in sublist:
        master_hydros_list.append(j)
    sublist = [x + 212 for x in sublist]    
master_list = []
for i in range(15000):
    chromo = np.append(master_carbs_list[i],master_hydros_list[i])
    master_list.append(chromo)
system.add_chromophores(master_list, "donor")
system.compute_energies(dcut=8)
system.set_energies()
lifetimes = [1e-11,1e-10,1e-9]
temp = 300
system.run_kmc(lifetimes, temp, n_holes=500, verbose=1)

