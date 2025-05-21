# plastid


## Directory Organization

* [`evaluation/`](./evaluation/): `cutoff` evaluations

* [`Brassicaceae/Atha_1001.md`](Brassicaceae/Atha_1001.md)
  * *Arabidopsis thaliana* 1001 Genomes Project
  * 1134 - 705 genomes

* [`Gramineae/Osat_50.md`](Gramineae/Osat_50.md)
  * *Oryza sativa* 50 accessions
  * 50 genomes

* [`Leguminosae/Mtru_hapmap.md`](Leguminosae/Mtru_384)
  * *Medicago truncatula* Hapmap Project
  * 384 - 156 genomes

* [`Solanaceae/Slyc_360.md`](Solanaceae/Slyc_360.md)
  * *Solanum lycopersicum* 360 accessions
  * 360 - 250 genomes

* [`Rosaceae/peach.md`](Rosaceae/peach.md)
  * Peaches (*Prunus* spp.)
  * 58 genomes

## Some other dependencies

* `mdtable2csv`

```shell script
mkdir -p ${HOME}/bin
curl -fsSL $(
    curl -fsSL https://api.github.com/repos/515hikaru/mdtable2csv/releases/latest |
        jq -r '.assets[] | select(.name == "mdtable2csv_linux_x86_64.tar.gz").browser_download_url'
    ) |
    tar xvz mdtable2csv
mv mdtable2csv ${HOME}/bin

```

* `slivar`

```shell script
mkdir -p ${HOME}/bin
curl -fsSL https://github.com/brentp/slivar/releases/download/v0.2.1/slivar -O
chmod +x slivar
mv slivar ${HOME}/bin

```
