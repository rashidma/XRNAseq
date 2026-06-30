# Custom R analysis image for the differential-expression + GSEA branch.
#
# No existing public image was confidently verified to contain this exact
# DESeq2 + fgsea + apeglm + annotation-package combination, so this image
# must be built and hosted by you (e.g. push to ghcr.io/<you>/xrnaseq-r-analysis:1.0.0)
# before modules/deseq2.nf and modules/gsea_fgsea.nf can run.
#
# Verify the base tag is still current before building:
#   https://hub.docker.com/r/bioconductor/bioconductor_docker/tags
FROM bioconductor/bioconductor_docker:3.23-r-4.6.0

RUN R -e "BiocManager::install(c('DESeq2', 'fgsea', 'apeglm', 'AnnotationDbi', 'org.Hs.eg.db'), update = FALSE, ask = FALSE)" \
    && R -e "install.packages(c('optparse', 'ggplot2', 'dplyr', 'readr', 'tibble'), repos = 'https://cloud.r-project.org')"

COPY bin/run_deseq2.R /usr/local/bin/run_deseq2.R
COPY bin/run_fgsea.R  /usr/local/bin/run_fgsea.R
RUN chmod +x /usr/local/bin/run_deseq2.R /usr/local/bin/run_fgsea.R
