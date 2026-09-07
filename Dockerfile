FROM rocker/tidyverse:4.3.2

# --- Системные зависимости (дополнительные) ---
RUN apt-get update && apt-get install -y \
    libglpk-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libnode-dev \
    npm \
    pandoc-citeproc \
    dos2unix \
    && rm -rf /var/lib/apt/lists/*

# --- CRAN пакеты ---
RUN R -e "install.packages(c('optparse','data.table','stringr','plotly','htmlwidgets','ggrepel'), repos='https://cloud.r-project.org')"

# --- Bioconductor пакеты ---
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('decoupleR','dorothea','limma'), update=FALSE, ask=FALSE)"

# --- Папка для скриптов ---
RUN mkdir -p /usr/local/my-scripts
COPY scripts/ /usr/local/my-scripts/

# --- Делаем их исполняемыми ---
RUN chmod +x /usr/local/my-scripts/*.R && dos2unix /usr/local/my-scripts/*.R

# --- Алиасы ---
RUN for f in /usr/local/my-scripts/*.R; do ln -s "$f" "/usr/local/bin/$(basename ${f%.R})"; done

# --- Рабочая директория ---
WORKDIR /usr/local/my-scripts

ENTRYPOINT ["/bin/bash"]