MS2Search is a Shiny app for the annotation from MS/MS (MS2) data by custom rules of the fragmentation pattern.

## Features

- Load MS2 spectra from MGF files (via the `Spectra` package)
- Rule-based lipid identification: rules are JSON objects (neutral losses, fragment ions, chain compositions)
- Interactive spectra view (neutral losses, chains, validations) and feature tables
- Custom Rule based Annotation file for MetaboScape >= 2026b
- Optional: spectral similarity and clustering with [DreaMS](https://github.com/pluskal-lab/DreaMS) embeddings (separate Docker service)
- Optional: import of feature tables from MetaboScape via its REST API 

## Installation

### Option 1: R package

```r
# install.packages("remotes")
remotes::install_github("cKeuters/MS2Search")
```

Start the app:

```r
MS2Search::run_app()
```

Without DreaMS, the cluster view is not available, everything else works.

### Option 2: Docker (with DreaMS)

The folder [`inst/deploy`](inst/deploy) contains a `Dockerfile` and a `docker-compose.yml` that start two services:

| Service     | Description                                                        | Port |
|-------------|--------------------------------------------------------------------|------|
| `ms2search` | The Shiny app (R, based on `rocker/verse`, packages via `renv`)    | 80   |
| `dreams`    | Python service (FastAPI/uvicorn) that computes DreaMS embeddings   | 8000 (internal) |

The app reaches the DreaMS service via `DREAMS_URL=http://dreams:8000`. The DreaMS model checkpoint is downloaded from Zenodo during the image build, so no download is needed at runtime.

```bash
git clone https://github.com/cKeuters/MS2Search.git
cd MS2Search
docker compose -f inst/deploy/docker-compose.yml up --build -d
```

If the DreaMS service is not reachable, the app keeps running and only the DreaMS features are disabled.

