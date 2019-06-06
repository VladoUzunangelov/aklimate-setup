### These are the steps for setting up an environment to run AKLIMATE.

1. Start up an instance of Ubuntu 18.04 LTS in m1.huge flavor.

2. Update the system via apt
```
sudo apt update
sudo apt upgrade -y
sudo apt dist-upgrade -y
sudo apt autoremove -y
```

3. install R via apt
```
sudo apt install -y r-base
```

4. install X11 via apt
```
sudo apt-get install -y xorg xvfb xauth xfonts-base
```

5. install the INTEL MKL library
    1. There is a handy script that makes this very easy. Get it from GitHub.
    ```
    git clone https://github.com/eddelbuettel/mkl4deb.git
    ```
    2. run the cloned script to install.
    ```
    sudo bash ./mkl4deb/script.sh
    ```
    3. confirmed installation in R
    ```
    sessionInfo()
    ```
    4. If you need, you can switch between linear algebra implementations.
    ```
    sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu
    ```

6. install other various system package dependencies for R packages in the subsequent steps
    ```
    sudo apt install -y libcurl4-openssl-dev libssl-dev libxml2-dev libmariadbclient-dev
    ```

7. in R, install many packages
```
install.packages(c("foreach","doParallel","ranger","plyr","abind","ROCR","caret","proxy","purrr","pracma","fastmatch","devtools","mlr","e1071","igraph","circlize","RColorBrewer","RMySQL"))
```

8. in R, install bioconductor following these [instructions](https://bioconductor.org/install/).
Following legacy instructions because we have R v3.4.4 (R<3.5.0)
```
source("https://bioconductor.org/biocLite.R")
```

9. in R, install ComplexHeatmap package
```
BiocInstaller::biocLite(c("ComplexHeatmap","FDb.InfiniumMethylation.hg19"))
```

10. Install Vlado's forked Similarity R package.
    1. Similarity requires RcppParallel, which needs the Intel TBB library.
    ```
    sudo apt install -y libtbb-dev
    ```
    2. Set the system environment variable:
    ```
    export tbb_os=linux
    ```
    3. Here are [instructions](https://cran.r-project.org/web/packages/githubinstall/vignettes/githubinstall.html) for using devtools to install R package hosted on GitHub.
    4. The package name is `VladoUzunangelov/similarity`.
    ```
    devtools::install_github("VladoUzunangelov/similarity")
    ```
