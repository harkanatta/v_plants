FROM rocker/tidyverse:latest

RUN apt-get update && apt-get install -y --no-install-recommends \
    libudunits2-0 libudunits2-dev \
    libgdal-dev gdal-bin \
    libgeos-dev libproj-dev \
    && ldconfig \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN echo '.libPaths(c("/home/rstudio/R/library", .libPaths()))' >> \
    /usr/local/lib/R/etc/Rprofile.site

EXPOSE 8787
CMD ["/init"]