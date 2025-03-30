# Base Image
FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Install R Packages
RUN R -e "BiocManager::install(version = '3.20'); BiocManager::valid(); BiocManager::install('BiocFileCache', update = TRUE, ask = FALSE); devtools::install_github('chiba-ai-med/wPGSAR', \
    upgrade='always', force=TRUE, INSTALL_opts = '--install-tests');\
    tools::testInstalledPackage('wPGSAR')"