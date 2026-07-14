####################################################################################################################################
# Adapted from following:
# - Rocker RStudio container using new versioned paradigm: https://github.com/rocker-org/rocker-versioned2/blob/master/dockerfiles/Dockerfile_rstudio_4.0.3
# - license: GPLV2
####################################################################################################################################

FROM rocker/rstudio:4.1.2

ENV S6_VERSION=v1.21.7.0
ENV RSTUDIO_VERSION=latest
ENV PATH=/usr/lib/rstudio-server/bin:$PATH

# install key dependencies of rstudio preview edition
# install radian via python and pip3
RUN apt-get update \
    && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
    libpq5 openssh-client openssh-server software-properties-common curl

RUN /rocker_scripts/install_pandoc.sh

# no login timeouts for local use
RUN echo "auth-timeout-minutes=0" >> /etc/rstudio/rserver.conf

# install key dependencies of certain packages that could be installed later
RUN apt-get update \
    && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends libxml2-dev libsodium-dev libssl-dev \
    imagemagick libmagick++-dev libgit2-dev libssh2-1-dev zlib1g-dev librsvg2-dev \
    libudunits2-dev libfontconfig1-dev libfreetype6-dev 

# install R packages needed for VSCode interaction and package management
RUN install2.r languageserver renv remotes

# install radian via python and pip3
RUN apt-get update \
    && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends python3-setuptools python3-pip

RUN pip3 install radian

# RUN echo "RENV_PATHS_CACHE=~/renv/cache" >> /usr/local/lib/R/etc/Renviron

#  adds addtional system dependencies needed for project
RUN apt-get update \
    && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
    glpk-utils libglpk-dev glpk-doc \
    libgtk2.0-dev libcairo2-dev xvfb xauth xfonts-base libxt-dev \
    libudunits2-dev libgdal-dev libgeos-dev libproj-dev \
    libzmq3-dev

# optional Copies the imagemagick config included here to expand the resource limits on 
# image magick the defaults are low and big heatmaps may been higher.
# COPY policy.xml /etc/ImageMagick-6/policy.xml

EXPOSE 8787

CMD ["/init"]
