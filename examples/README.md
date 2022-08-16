# Examples

These example notebooks demonstrate the workflow employed in morphct.
To create the conda environment:
```bash
conda env create -f environment-ex.yml
conda activate morphct-ex
cd ../
pip install .
```

To run the examples:
```bash
cd examples
jupyter lab
```
# Outline of the examples

1. **transfer-integral-workflow.ipynb** : Calculate an organic molecule's HOMO/LUMO and Transfer Integral (electronic coupling) for various orientations with morphCT.
2. **chromophore-picking.ipynb** : How to use VMD to select atoms from a morphology and create an array of the chosen atom's indices to create an array of chromophore indices
3. **new-acceptor-molecule-workflow.ipynb** : This workflow seeks to outline the path of least resistance for investigation of a new acceptor molecule with morphCT. To predict the mobililty of a morphology, we need 2 things: (1) a gsd of your morphology and (2) and list of atom indices that belong to each chromophore. This workflow automates the creation of (2). I WOULD START BY RUNNING THIS TO SEE IF MORPHCT IS WORKING. 
4. **workflow-itic.ipynb** : A notebook for breaking ITIC into multiple chromophores. This notebook also contains an example for creating a csv containing chromophore atom ids for using morphct-flow.  
5. **workflow-p3ht.ipynb** : How to use smarts matching to delineate chromophores. 

# You can find some further (uncurated) workflows that use the morphct api here:
https://github.com/JimmyRushing/thesis/tree/main/notebooks

(1) **fused-experiment.ipynb** : I used smiles strings to check out the homo-lumo gap of fused thiophenes. 
(2) **voronoi_figures.ipynb** : I created a plot that visualized the voronoi analysis that takes place in the creation of neighbor lists in morphCT
(3) other scratch workbooks for generating figures from morphCT data
