# Requirement to Install TeloTracker

## Environment telotracker: is the main environment

## Environment telotracker-clair3: contains clair3 (needs python 3.9) and some other basic sequence mapping packages in the defaul telotracker environment


# Once the telotracker environment is created, activate the environment and run the following commands:

## 1. To fully install quast:

quast-download-gridss
quast-download-silva

### Change the dataset web address for Busco
### Works as of 11/19/2024

OLD_BUSCO_BACTERIA_DB="https://busco-archive.ezlab.org/v3/datasets/bacteria_odb9.tar.gz"
OLD_BUSCO_FUNGI_DB="https://busco-archive.ezlab.org/v3/datasets/fungi_odb9.tar.gz"
OLD_BUSCO_EUKARYOTA_DB="https://busco-archive.ezlab.org/v3/datasets/eukaryota_odb9.tar.gz"

NEW_BUSCO_BACTERIA_DB="https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.2024-01-08.tar.gz"
NEW_BUSCO_FUNGI_DB="https://busco-data.ezlab.org/v5/data/lineages/fungi_odb10.2024-01-08.tar.gz"
NEW_BUSCO_EUKARYOTA_DB="https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2024-01-08.tar.gz"

BUSCO_INSTALL_PATH="$CONDA_PREFIX/lib/python3.*/site-packages/quast_libs/run_busco.py"

sed -i "s|$OLD_BUSCO_BACTERIA_DB|$NEW_BUSCO_BACTERIA_DB|g" $BUSCO_INSTALL_PATH
sed -i "s|$OLD_BUSCO_FUNGI_DB|$NEW_BUSCO_FUNGI_DB|g" $BUSCO_INSTALL_PATH
sed -i "s|$OLD_BUSCO_EUKARYOTA_DB|$NEW_BUSCO_EUKARYOTA_DB|g" $BUSCO_INSTALL_PATH

# sed -i 's|https://busco-archive.ezlab.org/v3/datasets/bacteria_odb9.tar.gz|https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.2024-01-08.tar.gz|g' $BUSCO_INSTALL_PATH
# sed -i 's|https://busco-archive.ezlab.org/v3/datasets/fungi_odb9.tar.gz|https://busco-data.ezlab.org/v5/data/lineages/fungi_odb10.2024-01-08.tar.gz|g' $BUSCO_INSTALL_PATH
# sed -i 's|https://busco-archive.ezlab.org/v3/datasets/eukaryota_odb9.tar.gz|https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2024-01-08.tar.gz|g' $BUSCO_INSTALL_PATH

quast-download-busco


### Quast install notification

"""
Executing transaction: - The default QUAST package does not include:
* GRIDSS (needed for structural variants detection)
* SILVA 16S rRNA database (needed for reference genome detection in metagenomic datasets)
* BUSCO tools and databases (needed for searching BUSCO genes) -- works in Linux only!

To be able to use those, please run
    quast-download-gridss
    quast-download-silva
    quast-download-busco
                          
"""

## 2. To install basecalling software:

### Dowload the current version of dorado
### Current version as of 11/19/2024: Dorado 0.8.3 @ https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.8.3-linux-x64.tar.gz

wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.8.3-linux-x64.tar.gz
gunzip https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.8.3-linux-x64.tar.gz
tar -xvf https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.8.3-linux-x64.tar

## To Do ##
## Also need pod5 -- ** Need to add pip and pip pod5 in telotracker environment
