# Base Image
FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Install R Packages
RUN R -e "devtools::install_github('chiba-ai-med/wPGSA', \
    upgrade='always', force=TRUE, INSTALL_opts = '--install-tests');\
    tools::testInstalledPackage('wPGSA')"