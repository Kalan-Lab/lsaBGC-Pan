# This code is adapted from BiG-SCAPE's Dockerfile - if you are using lsaBGC-Pan - you should definitely check out
# BiG-SCAPE/CORASON incase it suites your needs better (https://www.nature.com/articles/s41589-019-0400-9) - eg.
# you are interested in clustering diverse gene-clusters based on protein domain simlarity or are interested
# in investigating the variability of contexts for a single reference gene! 

FROM continuumio/miniconda3
LABEL maintainer="Rauf Salamzade - Kalan Lab, UW-Madison"

WORKDIR /usr/src
SHELL ["/bin/bash", "-c"]

# Clone lsaBGC-Pan v1.1.0 with multiprocessing in place of concurrent.futures from raufs 
# github repo and create conda environment
RUN apt-get update && apt-get install -y git wget gcc && \
    conda install -n base conda-libmamba-solver && conda config --set solver libmamba && \
    conda create -p /usr/src/lsaBGC_conda_env/ -c conda-forge -c bioconda lsabgc panaroo && \
    conda clean --all -y && \
    echo "source activate /usr/src/lsaBGC_conda_env/" > ~/.bashrc && source ~/.bashrc && \
    apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

# Set env variables and download minimal databases (MIBiG v3.1 + PGAP HMMs)
ENV PATH /usr/src/lsaBGC_conda_env/bin:$PATH
ENV ZOL_EXEC_PATH /usr/src/lsaBGC_conda_env/bin/
ENV ZOL_DATA_PATH /usr/src/lsaBGC_conda_env/share/zol-1.4.12/db/
RUN setup_annotation_dbs.py -ld && chmod 777 /home
USER 1000:1000
RUN mkdir /home/input /home/output
WORKDIR /home
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

ENTRYPOINT ["lsaBGC-Pan"]
CMD ["--help"]
