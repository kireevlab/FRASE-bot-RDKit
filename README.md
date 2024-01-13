# FRASE-bot-RDKit

FRASE-bot-RDKit is a collection of cheminformatics protocols designed for data processing on FRASE-bot screening, developed in Python with RDKit, a powerful cheminformatics toolkit. This version eliminates the need for a Pipeline Pilot commercial license and provides flexibility for users to leverage the capabilities of RDKit in their cheminformatics workflows.

## Features

- **Python and RDKit:** FRASE-bot-RDKit harnesses the capabilities of RDKit to perform a variety of cheminformatics tasks in a Python environment.

- **Elimination of Pipeline Pilot Requirement:** Unlike the original FRASE-bot version, FRASE-bot-RDKit operates independently of Pipeline Pilot, allowing users to run protocols without the need for a commercial license.

- **Two Protocols:** The workflow is divided into two distinct protocols, each focusing on specific cheminformatics tasks.

## System Requirements

To use FRASE-bot-RDKit, ensure that you have the following prerequisites:

- **Python:** Install Python on your system. The code is compatible with Python 3.x.

- **RDKit:** Install the RDKit cheminformatics toolkit. Detailed instructions can be found on the [RDKit website](https://www.rdkit.org/).
  
## Installation 
We encourage to install RDKit on Linux using conda. Current scripts were specially designed in Linux enviroment to run calculations for improving computational efficiency.  

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
   This bash script, FRASE-screening.sh, use python script FRASE-screening.py to perform batch FRASE screening based on 10,599 FARASEs. Note that certain High-Performance Computing (HPC) systems may have limitations on job submission numbers. Adjust the submission array accordingly in FRASE-screening.sh: `##SBATCH --array 0-10598`. In addition, modyfy the path of FRASE_database_PDB/ directory in FRASE-screening.py: `pdb_files = glob.glob(f'/your_path/FRASE_database_PDB{folder_id}/{frase_inp}')`
3. After execution, generate new files for subsequent machine learning using the following bash script.
   ```
   submit generate_IF_target.sh
   ```
4. Employ the `neural-network.py` script in the final step to generate screening outcomes. 
   ```
   submit neural-network.sh
   ```
   The information on the interaction fingerprint of 10,598 FRASEs is stored in the file `IF_FRASEdb.json` and is utilized as training data in the neural-network.py script. Save this file to your working directory and specify its actual path in neural-network.py: `FRASEdb_file = '/your_path/IF_FRASEdb.json'`
