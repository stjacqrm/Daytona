# run_daytona.sh and Daytona-Combine
This wrapper is used to tidy and analyse Clear Labs data using a custom cleaner and the [`Daytona-combine pipeline`](https://github.com/BPHL-Molecular/Daytona_combine). 


## Prerequisites

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install Python3.7+ and BioPython

May need to use update-alternatives to set correct version of python, ie:

```
$ sudo update-alternatives --install /usr/bin/python3 python3 /home/${USER}/miniconda3/bin/python3.8 4
```

BioPython:

```
$ python3 -m pip install biopython
```

3. Install Bash Needs (may need to install individually)

```
$ sudo apt-get update && sudo apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    pkg-config
```

```
$ sudo apt-get install libtool-bin
```

4. Install Go

```
$ export VERSION=1.11 OS=linux ARCH=amd64 && wget https://dl.google.com/go/go1.22.2.linux-amd64.tar.gz && sudo tar -C /usr/local -xzvf go1.22.2.linux-amd64.tar.gz && rm go1.22.2.linux-amd64.tar.gz
```

```
$ echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc
```

5. Install [`Singularity`](https://singularity-tutorial.github.io/01-installation/)

```
$ wget https://github.com/sylabs/singularity/releases/download/v4.1.2/singularity-ce-4.1.2.tar.gz
$ tar -xzf singularity-ce-4.1.2.tar.gz
$ cd ./singularity-ce-4.1.2
$ ./mconfig
$ make -C ./builddir
$ sudo make -C ./builddir install
```

6. Install the Code

```
$ git clone https://github.com/VADGS/Daytona.git
```

Set to path in bashrc:

```
export PATH="/home/${USER}/Applications/daytona:$PATH"
```

## How to run
### 1) If the data is "raw" from Clear Labs:
1. Download the tar file from Clear Labs
2. Move the tar file into the desired directory
3. Use the following command to invoke the pipeline:

```
$ run_daytona.sh -t tar
```

### 2) If the data has been previously downloaded and extracted/analzed from Clear Labs:
1. Use the following command to invoke the pipeline:

```
$ run_daytona.sh -r path/to/directory -o desired/output/name
```

## Results
All results can be found in the directory /output.

