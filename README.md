# FRASE-bot-RDKit

FRASE-bot-RDKit is a collection of cheminformatics protocols designed for data processing on FRASE-bot screening, developed in Python with RDKit, a powerful cheminformatics toolkit. This version eliminates the need for a Pipeline Pilot commercial license and provides flexibility for users to leverage the capabilities of RDKit in their cheminformatics workflows.

## Features

- **Python and RDKit:** FRASE-bot-RDKit made use of the capabilities of RDKit to perform a variety of cheminformatics tasks in a Python environment.

- **Elimination of Pipeline Pilot Requirement:** Unlike the original FRASE-bot version, FRASE-bot-RDKit operates independently of Pipeline Pilot, allowing users to run protocols without the need for a commercial license.

- **Two Protocols:** The workflow is divided into two distinct protocols, one for FRASE-bot screening processing, and the other for subsequent neural network for final screening.
- **Efficiency:** Current scripts including bash and python scripts are specifically designed for optimal performance in a Linux environment. This ensures efficient execution of calculations.

## System Requirements

To use FRASE-bot-RDKit, ensure that you have the following prerequisites:

- **Python:** Install Python on your system. The code is compatible with Python 3.x.

- **RDKit:** Install the RDKit cheminformatics toolkit. Detailed instructions can be found on the [RDKit website](https://www.rdkit.org/).
  
## Installation 
We strongly recommend installing RDKit on Linux using conda. Alternatively, you can run the scripts on Jupyter Notebook.

1. **To install RDKit cheminformatics toolkit, it's recommended to create a new conda environment for RDKit.**
  ```
  conda create --name rdkit-tools python=3.8
  ```
2. **Activate the rdkit enviroment:**
  ```
  conda activate rdkit-tools
  ```
3. **Install RDKit:**
  ```
  conda install -c conda-forge rdkit
  ```

## Workflow:
1. FRASE database is stored in the FRASE_database_PDB/ folder, utilized for FRASE-bot screening. The folder contains a total of 10,599 FRASEs.
2. Execute the main script:
   ```
   submit FRASE-screening.sh.
   ```
   This bash script, `FRASE-bot-RDKit-batch.sh`, use python script `FRASE-bot-RDKit.py` to perform batch FRASE screening based on 10,599 FARASEs. Note that certain High-Performance Computing (HPC) systems may have limitations on job submission numbers. Adjust the submission array accordingly in `FRASE-bot-RDKit-batch.sh`: `##SBATCH --array 0-10598`. In addition, modyfy the path of FRASE_database_PDB/ directory in `FRASE-bot-RDKit.py`: `pdb_files = glob.glob(f'/your_path/FRASE_database_PDB{folder_id}/{frase_inp}')`
3. After execution, generate new files for subsequent machine learning using the following bash script.
   ```
   submit generate_IF_target.sh
   ```
4. Employ the `FRASE-bot-RDKit-NN.py` script in the final step to generate screening outcomes. 
   ```
   submit FRASE-bot-RDKit-NN.sh
   ```
   The information on the interaction fingerprint of 10,598 FRASEs is stored in the file `IF_FRASEdb.json` and is utilized as training data. Save this file to your working directory and specify its actual path in `FRASE-bot-RDKit-NN.py`: `FRASEdb_file = '/your_path/IF_FRASEdb.json'`
