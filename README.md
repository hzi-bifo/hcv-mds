# Neutralization plots for HCV

## Requirements
* Clone this repository:
  ```
  git clone https://github.com/hzi-bifo/hcv-mds 
  ```
* ```Python 2``` with the following packages:
    * **numpy**
    * **scipy**
    * **matplotlib**
    * **pandas**
    * **sklearn**
    * **gapkmean**
   
  If these are not installed, you can install them with ``` pip ```. 
  ```
   cd ./scripts
   pip2 install -r requirements.txt 
   ```
* ```R version 3.4 or higher``` with the following packages:
    * **plotrix**
    * **smacof**
    * **rgl**
    * **polynom**
    * **Hmisc**
  
    If these are not installed, you can install them by running:
    ```
    R
    install.packages(c("plotrix", "smacof", "rgl", "polynom", "Hmisc"))
    ```
    
## Usage
* Put the file containing neutralization matrix (residual infectivity at maximum concentration) in the ```input``` folder. The file should be a comma separated **csv** file with the viruses as columns and sera as rows. Use **NA** for missing values. Here is an example:
  
  Sera / Viruses | Gt1b(Con1) |  Gt2a(Jc1) | Gt2a(2a-3) | Gt2b(J8)
  :---: | :----------: | :----------: | :----------: | :--------:
  **37**| 80.2895218854| 57.5271175607| 74.9637533368|80.3419066704
  **39**| 65.5662622736| 72.1845541645| 70.7020142884|71.0087697487
  **41**| 84.0995706752| 50.7927096018| 41.204264273 |46.658115028
  **46**| 83.2643944991| 76.0598932927| 74.5105488987|NA
  **47**| 70.9418542308| 69.7642736513| 75.3173168277|49.7798284265

* Go in the ```scripts``` folder and run the ```run_r.py``` file.
  ```
  cd scripts
  python run_r.py data.csv 
  ```
  There are different flags for selecting the weighting scheme, normalization method and number of dimensions for the multi-dimensional scaling.
  ```
  python run_r.py inputfile.csv [-h ] [-dim ] [-w ] [-norm ] [-loo ]
  
  -h, --help          :  show this help message and exit
  
  -norm, --normalize  :  Normalization method 
                         0 - No Normalization (Default) 
                         1 - take column minimum 
                         2 - take global minimum 
                       >=3 - minimum column basis 
  
  -w, --weight WEIGHT :  Weighing method 
                         0 - No weights (Default)
                         1 - 1/dij as weights 
  
  -dim, --dim DIM     :  Number of dimensions 
                         2 - Two dimensions (Default) 
                         3 - Three dimensions
                         4 or higher
  
  -loo, --loo LOO     :  Leave-one-out tests 
                         0 - No loo test (Default)) 
                         1 - Compute loo error
  ```
 * After running this, the Neutralization map will be displayed and you'll be prompted to enter the value of K for K-means clustering. You can select the K by looking at the gap-statistics plot which will also be displayed along with the neutralization map. Once you enter the value of K, the final map with clustered viruses will also be displayed. When you close all the plots, program will automatically exit.
 
## Results
* The log file for the entire process can be found in ```log``` folder with the name ```log.txt```.
* The final plots and the summary of the result will be in the ```results``` folder.
  *inputfile_nmds_plot.eps* : eps file containing the neutralization plot of the input data
  *inputfile_nmds_plot.png* : png file containing the neutralization plot of the input data
  *inputfile_summary.txt*   : txt file containing the stress values, loo error, cluster information etc.
  *inputfile_cluser.png*    : plot of final clusters of the viruses
* The normalization files, LOO-error files, clustering analysis files etc. will be in the ```processing``` folder.




