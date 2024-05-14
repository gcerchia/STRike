# STRike plugin user guide

This is a plugin for the Ensembl Variant Effect Predictor (VEP) that generates HGVS notation at genomic and transcript level for Short Tandem Repeats. The plugin requires the installation of Ensembl VEP and Perl (version used: v5.34.0).

## Installation

First of all, Ensembl VEP must be installed:

```bash
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl
```

For more [detail about Ensembl VEP installation](https://github.com/Ensembl/ensembl-vep).

Then, download STRike and place it under VEP plugins folder:

```bash
git clone https://github.com/gcerchia/STRike
cp STRike.pm ~/.vep/Plugins/
```

If your VEP cache directory is not set to default, then:
```bash
cp STRike.pm /ensembl-vep/{VERSION}/{CACHE_DIR}/Plugins/
```

## Usage

For more detailed information about the usage of VEP plugins, please check [this documentation](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html).

In order to annnotate a VCF file, the following command can be executed:

```bash
vep -i input.vcf --plugin STRike --vcf --force --assembly GRCh38 --cache --port 3337
```

If your VEP configurations are different than default, then tou need to specify the cache directory:

```bash
vep -i input.vcf --plugin STRike --vcf --force --assembly GRCh38 --port 3337 --cache --dir /ensembl-vep/{VERSION}/{CACHE_DIR}/Plugins/
```

For further details about VEP plugin parameters, please check [this documentation](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_plugin).
