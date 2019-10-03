# KnowEnG's Network Prepper
 This is the Knowledge Engine for Genomics (KnowEnG), an NIH BD2K Center of Excellence, Network Prepper tool.
This tool **prepares** the data of a user-supplied network for subsequent processing by KnowEnG Analytics Platform.

## Detailed assumptions and cleanup steps

The original network input file (`raw_edgefile_full_path`) is assumed to be a tab separated file with at least three columns 
where each row contains information about an edge in the network.  
The first column is treated as the source nodes of the network, the second column is treated as the target nodes.
The third column is required and assumed to have the weights of each edge.  The weights are assumed to be non-null and greater than zero.

The Network Prepper tool begins by checking that the above assumptions are met.  If not, it will FAIL and return an ERROR.

It then used the redis entity database (`redis_credential`) of the KnowEnG Knowledge Network to map the source and target node identifiers to 
stable Ensembl identifiers using the species identifier (`taxonid`) and (`source_hint`). 

The tool then checks that a sufficient percentage of entities were successfully mapped (at least `network_threshold` of the original nodes and the original edges).

Finally, the tool performs some network cleanup.  It makes each edge bi-directional if `make_symmetric` and removes any duplicates. 
If there are fewer clean edges than clean nodes, the tool FAILs and returns an ERROR.

The outputs of the tools are:

#### A) *.full_mapped_edges.tsv
- This file is just the original network file with four columns prepended:
  1. source_mapped_id
  2. source_mapped_alias
  3. target_mapped_id
  4. target_mapped_alias
Columns will often contain the mapped value, or if not mapped, the reason mapping failed ("unmapped-none" if no match was found, or "unmapped-many" if multiple matches were found).

#### B) *.full_mapped_nodes.tsv
- This file lists all of the original node names as well as their mapped values:
  1. original input identifier
  2. Knowledge Network mapped identifier
  3. type of entity, 'Gene' or 'Property'
  4. official symbol or alias
  5. gene or property description
  6. gene biotype (e.g. protein coding)
Columns will often contain the mapped value, or if not mapped, the reason mapping failed ("unmapped-none" if no match was found, or "unmapped-many" if multiple matches were found).

#### C) *.clean.edge
- This file contains only the cleaned edge list after unmapped edges are dropped, it is made symmetric (if required), and only the largest weighted edge of duplicates is kept.  This input is ready for other KnowEnG pipelines.

#### D) *.clean.node_map
- This file is in the sample formate as output B but contains only the nodes that were successfully mapped.

#### E) *.metadata 
This yaml file contains information about the prepared network. Its keys include summarizations about the network size, information about the meaning of its edges, and some commands and configurations used in its construction.


* * * 
## How to run this pipeline with your data
* * * 

### 1. Clone the Network_Prepper Repo
```
 git clone https://github.com/KnowEnG/Network_Prepper.git
```
 
### 2. Install the following (Ubuntu or Linux)
```
 apt-get install -y python3-pip
 apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran
 pip3 install numpy
 pip3 install pandas
 pip3 install scipy==0.19.1
 pip3 install scikit-learn==0.19.2
 apt-get install -y libfreetype6-dev libxft-dev
 pip3 install xmlrunner
 pip3 install pyyaml
 pip3 install knpackage
 pip3 install redis
```

### 3. Change directory to Network_Prepper

```
cd Network_Prepper 
```

### 4. Create your run directory

 ```
 mkdir run_directory
 ```

### 5. Change directory to the run_directory

 ```
 cd run_directory
 ```

### 6. Create your results directory

 ```
 mkdir results_directory
 ```
 
### 7. Create run_paramters file  (YAML Format)
 ``` 
Look for examples of run_parameters in ./Network_Prepper/data/run_files/TEST_0_small_success.yml
 ```
### 8. Modify run_paramters file  (YAML Format)
```
set the raw_edgefile_full_path to point to your data
```

### 9 Run the Network Prepper tool:

  * Update PYTHONPATH enviroment variable
   ``` 
   export PYTHONPATH='../src':$PYTHONPATH    
   ```
   
  * Run (these relative paths assume you are in the test directory with setup as described above)
   ```
  python3 ../src/network_prep.py -run_directory ./run_dir -run_file TEST_user_job.yml
   ```

* * * 
## Description of "run_parameters" file
* * * 

| **Key**                    | **Value**                            | **Comments**                                      |
| -------------------------- | ------------------------------------ | ------------------------------------------------- |
| raw_edgefile_full_path     | directory+network_data_name          | Path and file name of user network                |
| make_symmetric             | boolean                              | True or False to make network symmetricadsheet    |
| network_threshold          | float (0.6)                          | Proportion of nodes and edges that must map       |
| results_directory          | directory                            | Directory to save the output files                |
| redis_credential           | host, password and port              | Credential to access gene names lookup            |
| taxonid                    | 9606                                 | Taxon id of the genes                             |
| source_hint                | ' '                                  | Hint for lookup ensembl names                     |

