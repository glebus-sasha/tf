FROM rocker/r-ver:4.3.2

# --- Системные зависимости ---
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libglpk-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libx11-dev \
    libxt-dev \
    libcairo2-dev \
    libgit2-dev \
    libicu-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpq-dev \
    libsqlite3-dev \
    libssh2-1-dev \
    zlib1g-dev \
    pandoc \
    pandoc-citeproc \
    dos2unix \
    && rm -rf /var/lib/apt/lists/*

# --- CRAN пакеты (устанавливаем в правильном порядке) ---
RUN R -e "install.packages(c('optparse','data.table','dplyr','tidyr','tibble','readr','stringr'), repos='https://cloud.r-project.org')"

RUN R -e "install.packages(c('ggplot2','ggrepel'), repos='https://cloud.r-project.org')"

RUN R -e "install.packages(c('htmlwidgets','plotly'), repos='https://cloud.r-project.org')"

# --- Bioconductor пакеты ---
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('decoupleR','dorothea','limma'), update=FALSE, ask=FALSE)"

# --- Папка для твоих скриптов, отдельно от системных бинарников ---
RUN mkdir -p /usr/local/my-scripts
COPY scripts/ /usr/local/my-scripts/

# --- Делаем их исполняемыми и убираем возможные CRLF (Windows) ---
RUN chmod +x /usr/local/my-scripts/*.R && dos2unix /usr/local/my-scripts/*.R

# --- Алиасы для удобного вызова без .R ---
RUN for f in /usr/local/my-scripts/*.R; do ln -s "$f" "/usr/local/bin/$(basename ${f%.R})"; done

# --- Рабочая директория ---
WORKDIR /usr/local/my-scripts

# --- По умолчанию интерактивный shell ---
ENTRYPOINT ["/bin/bash"]