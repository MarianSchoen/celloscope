FROM rocker/shiny:4.1.0

RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y apt-utils software-properties-common locales

RUN DEBIAN_FRONTEND=noninteractive
RUN apt-get -qq install -y \
    cmake \
    libxml2-dev \
    fftw-dev \
    libfftw3-dev \
    fftw3 \
    libcurl4-openssl-dev \
    libssl-dev \
    tree \
    gcc \
    g++ \
    git \
    libgeos-dev

RUN apt-get clean && apt-get autoremove -y

WORKDIR /usr/local/r-job/

COPY R/install_scripts/install_cran_packages.R install_cran_packages.R
RUN Rscript install_cran_packages.R

COPY R/install_scripts/install_bioconductor_packages.R install_bioconductor_packages.R
RUN Rscript install_bioconductor_packages.R

COPY R/install_scripts/install_seurat.R install_seurat.R
RUN Rscript install_seurat.R 


# # add those R-package repositories from our gitlab.spang-lab.de website, with a specific tag:
# # the branch argument takes a tag aswell:
RUN git clone --branch v-0.2.1 https://jupyter:ywsQosELJUQUVk4g7-My@gitlab.spang-lab.de/jsimeth/dtdpipeline.git
RUN git clone --branch v1.0 https://jupyter:fztP7WY-hkDKyeSrsfey@gitlab.spang-lab.de/jsimeth/diggeR.git
COPY R/install_scripts/install_gitlab-spang-lab_packages.R install_gitlab-spang-lab_packages.R
RUN Rscript install_gitlab-spang-lab_packages.R

COPY . .
COPY . /srv/shiny-server/

RUN chown -R shiny /srv/shiny-server/data 

USER shiny

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
