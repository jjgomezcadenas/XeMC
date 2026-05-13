# Running Recipe

Complete commands for the three-stage background pipeline.
Run from `lxe_mc/` with the `sources` branch checked out.

## Prerequisites

```bash
cd /Users/jjgomezcadenas/Projects/XeMC/lxe_mc
conda activate itacatf   # for python scripts
```

## Directory layout

```
results/
  flux_tables/
    cryostat/{barrel,top,bottom}/{Bi214,Tl208}/
    pmt/{top,bottom,bottom_lxe,barrel}/{Bi214,Tl208}/
  backgrounds/
    cryostat/{barrel,top,bottom}/{Bi214,Tl208}/
    pmt/{top,bottom,bottom_lxe,barrel}/{Bi214,Tl208}/
```


## Stage 1: Flux table generation

Use `julia -t 8` for 8-thread parallelism. Adjust `--n` for statistics
(1M recommended for cryostat, 100k sufficient for PMTs since survival = 50%).

### Cryostat sources

```bash
# Barrel
julia -t 8 --project=LXeMC LXeMC/scripts/generate_flux_tables.jl \
  --source cryostat_barrel --isotope Bi214 --n 1000000 --seed 42 \
  --outdir LXeMC/results/flux_tables/cryostat/barrel/Bi214

julia -t 8 --project=LXeMC LXeMC/scripts/generate_flux_tables.jl \
  --source cryostat_barrel --isotope Tl208 --n 1000000 --seed 42 \
  --outdir LXeMC/results/flux_tables/cryostat/barrel/Tl208

# Top
julia -t 8 --project=LXeMC LXeMC/scripts/generate_flux_tables.jl \
  --source cryostat_top --isotope Bi214 --n 1000000 --seed 42 \
  --outdir LXeMC/results/flux_tables/cryostat/top/Bi214

julia -t 8 --project=LXeMC LXeMC/scripts/generate_flux_tables.jl \
  --source cryostat_top --isotope Tl208 --n 1000000 --seed 42 \
  --outdir LXeMC/results/flux_tables/cryostat/top/Tl208

# Bottom
julia -t 8 --project=LXeMC LXeMC/scripts/generate_flux_tables.jl \
  --source cryostat_bottom --isotope Bi214 --n 1000000 --seed 42 \
  --outdir LXeMC/results/flux_tables/cryostat/bottom/Bi214

julia -t 8 --project=LXeMC LXeMC/scripts/generate_flux_tables.jl \
  --source cryostat_bottom --isotope Tl208 --n 1000000 --seed 42 \
  --outdir LXeMC/results/flux_tables/cryostat/bottom/Tl208
```

### PMT sources

```bash
# PMT top
julia -t 8 --project=LXeMC LXeMC/scripts/generate_flux_tables.jl \
  --source pmt_top --isotope Bi214 --n 100000 --seed 42 \
  --outdir LXeMC/results/flux_tables/pmt/top/Bi214

julia -t 8 --project=LXeMC LXeMC/scripts/generate_flux_tables.jl \
  --source pmt_top --isotope Tl208 --n 100000 --seed 42 \
  --outdir LXeMC/results/flux_tables/pmt/top/Tl208

# PMT bottom
julia -t 8 --project=LXeMC LXeMC/scripts/generate_flux_tables.jl \
  --source pmt_bottom --isotope Bi214 --n 100000 --seed 42 \
  --outdir LXeMC/results/flux_tables/pmt/bottom/Bi214

julia -t 8 --project=LXeMC LXeMC/scripts/generate_flux_tables.jl \
  --source pmt_bottom --isotope Tl208 --n 100000 --seed 42 \
  --outdir LXeMC/results/flux_tables/pmt/bottom/Tl208

# PMT bottom (with passive LXe propagation)
julia -t 8 --project=LXeMC LXeMC/scripts/generate_flux_tables.jl \
  --source pmt_bottom_lxe --isotope Bi214 --n 100000 --seed 42 \
  --outdir LXeMC/results/flux_tables/pmt/bottom_lxe/Bi214

julia -t 8 --project=LXeMC LXeMC/scripts/generate_flux_tables.jl \
  --source pmt_bottom_lxe --isotope Tl208 --n 100000 --seed 42 \
  --outdir LXeMC/results/flux_tables/pmt/bottom_lxe/Tl208

# PMT barrel
julia -t 8 --project=LXeMC LXeMC/scripts/generate_flux_tables.jl \
  --source pmt_barrel --isotope Bi214 --n 100000 --seed 42 \
  --outdir LXeMC/results/flux_tables/pmt/barrel/Bi214

julia -t 8 --project=LXeMC LXeMC/scripts/generate_flux_tables.jl \
  --source pmt_barrel --isotope Tl208 --n 100000 --seed 42 \
  --outdir LXeMC/results/flux_tables/pmt/barrel/Tl208
```

## Stage 2: Background processing (fast kernel + FV stack)

Adjust `--n` for statistics (500k-1M recommended).

### Cryostat sources

```bash
# Barrel
julia -t 8 --project=LXeMC LXeMC/scripts/run_source_backgrounds.jl \
  --source cryostat_barrel --isotope Bi214 --n 500000 --seed 42 \
  --indir LXeMC/results/flux_tables/cryostat/barrel/Bi214 \
  --outdir LXeMC/results/backgrounds/cryostat/barrel/Bi214

julia -t 8 --project=LXeMC LXeMC/scripts/run_source_backgrounds.jl \
  --source cryostat_barrel --isotope Tl208 --n 500000 --seed 42 \
  --indir LXeMC/results/flux_tables/cryostat/barrel/Tl208 \
  --outdir LXeMC/results/backgrounds/cryostat/barrel/Tl208

# Top
julia -t 8 --project=LXeMC LXeMC/scripts/run_source_backgrounds.jl \
  --source cryostat_top --isotope Bi214 --n 500000 --seed 42 \
  --indir LXeMC/results/flux_tables/cryostat/top/Bi214 \
  --outdir LXeMC/results/backgrounds/cryostat/top/Bi214

julia -t 8 --project=LXeMC LXeMC/scripts/run_source_backgrounds.jl \
  --source cryostat_top --isotope Tl208 --n 500000 --seed 42 \
  --indir LXeMC/results/flux_tables/cryostat/top/Tl208 \
  --outdir LXeMC/results/backgrounds/cryostat/top/Tl208

# Bottom
julia -t 8 --project=LXeMC LXeMC/scripts/run_source_backgrounds.jl \
  --source cryostat_bottom --isotope Bi214 --n 500000 --seed 42 \
  --indir LXeMC/results/flux_tables/cryostat/bottom/Bi214 \
  --outdir LXeMC/results/backgrounds/cryostat/bottom/Bi214

julia -t 8 --project=LXeMC LXeMC/scripts/run_source_backgrounds.jl \
  --source cryostat_bottom --isotope Tl208 --n 500000 --seed 42 \
  --indir LXeMC/results/flux_tables/cryostat/bottom/Tl208 \
  --outdir LXeMC/results/backgrounds/cryostat/bottom/Tl208
```

### PMT sources

```bash
# PMT top
julia -t 8 --project=LXeMC LXeMC/scripts/run_source_backgrounds.jl \
  --source pmt_top --isotope Bi214 --n 500000 --seed 42 \
  --indir LXeMC/results/flux_tables/pmt/top/Bi214 \
  --outdir LXeMC/results/backgrounds/pmt/top/Bi214

julia -t 8 --project=LXeMC LXeMC/scripts/run_source_backgrounds.jl \
  --source pmt_top --isotope Tl208 --n 500000 --seed 42 \
  --indir LXeMC/results/flux_tables/pmt/top/Tl208 \
  --outdir LXeMC/results/backgrounds/pmt/top/Tl208

# PMT bottom
julia -t 8 --project=LXeMC LXeMC/scripts/run_source_backgrounds.jl \
  --source pmt_bottom --isotope Bi214 --n 500000 --seed 42 \
  --indir LXeMC/results/flux_tables/pmt/bottom/Bi214 \
  --outdir LXeMC/results/backgrounds/pmt/bottom/Bi214

julia -t 8 --project=LXeMC LXeMC/scripts/run_source_backgrounds.jl \
  --source pmt_bottom --isotope Tl208 --n 500000 --seed 42 \
  --indir LXeMC/results/flux_tables/pmt/bottom/Tl208 \
  --outdir LXeMC/results/backgrounds/pmt/bottom/Tl208

# PMT bottom (with passive LXe propagation)
julia -t 8 --project=LXeMC LXeMC/scripts/run_source_backgrounds.jl \
  --source pmt_bottom_lxe --isotope Bi214 --n 500000 --seed 42 \
  --indir LXeMC/results/flux_tables/pmt/bottom_lxe/Bi214 \
  --outdir LXeMC/results/backgrounds/pmt/bottom_lxe/Bi214

julia -t 8 --project=LXeMC LXeMC/scripts/run_source_backgrounds.jl \
  --source pmt_bottom_lxe --isotope Tl208 --n 500000 --seed 42 \
  --indir LXeMC/results/flux_tables/pmt/bottom_lxe/Tl208 \
  --outdir LXeMC/results/backgrounds/pmt/bottom_lxe/Tl208

# PMT barrel
julia -t 8 --project=LXeMC LXeMC/scripts/run_source_backgrounds.jl \
  --source pmt_barrel --isotope Bi214 --n 500000 --seed 42 \
  --indir LXeMC/results/flux_tables/pmt/barrel/Bi214 \
  --outdir LXeMC/results/backgrounds/pmt/barrel/Bi214

julia -t 8 --project=LXeMC LXeMC/scripts/run_source_backgrounds.jl \
  --source pmt_barrel --isotope Tl208 --n 500000 --seed 42 \
  --indir LXeMC/results/flux_tables/pmt/barrel/Tl208 \
  --outdir LXeMC/results/backgrounds/pmt/barrel/Tl208
```

## Stage 3: Analysis (python)

### Plot flux tables

```bash
# Cryostat
python py/plot_fluxes.py LXeMC/results/flux_tables/cryostat/barrel/Bi214
python py/plot_fluxes.py LXeMC/results/flux_tables/cryostat/barrel/Tl208
python py/plot_fluxes.py LXeMC/results/flux_tables/cryostat/top/Bi214
python py/plot_fluxes.py LXeMC/results/flux_tables/cryostat/top/Tl208
python py/plot_fluxes.py LXeMC/results/flux_tables/cryostat/bottom/Bi214
python py/plot_fluxes.py LXeMC/results/flux_tables/cryostat/bottom/Tl208

# PMT
python py/plot_fluxes.py LXeMC/results/flux_tables/pmt/top/Bi214
python py/plot_fluxes.py LXeMC/results/flux_tables/pmt/top/Tl208
python py/plot_fluxes.py LXeMC/results/flux_tables/pmt/bottom/Bi214
python py/plot_fluxes.py LXeMC/results/flux_tables/pmt/bottom/Tl208
python py/plot_fluxes.py LXeMC/results/flux_tables/pmt/bottom_lxe/Bi214
python py/plot_fluxes.py LXeMC/results/flux_tables/pmt/bottom_lxe/Tl208
python py/plot_fluxes.py LXeMC/results/flux_tables/pmt/barrel/Bi214
python py/plot_fluxes.py LXeMC/results/flux_tables/pmt/barrel/Tl208
```

### Plot candidates

```bash
# Cryostat
python py/plot_candidates.py LXeMC/results/backgrounds/cryostat/barrel/Bi214
python py/plot_candidates.py LXeMC/results/backgrounds/cryostat/barrel/Tl208
python py/plot_candidates.py LXeMC/results/backgrounds/cryostat/top/Bi214
python py/plot_candidates.py LXeMC/results/backgrounds/cryostat/top/Tl208
python py/plot_candidates.py LXeMC/results/backgrounds/cryostat/bottom/Bi214
python py/plot_candidates.py LXeMC/results/backgrounds/cryostat/bottom/Tl208

# PMT
python py/plot_candidates.py LXeMC/results/backgrounds/pmt/top/Bi214
python py/plot_candidates.py LXeMC/results/backgrounds/pmt/top/Tl208
python py/plot_candidates.py LXeMC/results/backgrounds/pmt/bottom/Bi214
python py/plot_candidates.py LXeMC/results/backgrounds/pmt/bottom/Tl208
python py/plot_candidates.py LXeMC/results/backgrounds/pmt/bottom_lxe/Bi214
python py/plot_candidates.py LXeMC/results/backgrounds/pmt/bottom_lxe/Tl208
python py/plot_candidates.py LXeMC/results/backgrounds/pmt/barrel/Bi214
python py/plot_candidates.py LXeMC/results/backgrounds/pmt/barrel/Tl208
```

### Compute background rates (smearing + ROI)

Default: sigma=25 keV, ROI=50 keV around Q_bb=2.458 MeV.

```bash
# Cryostat
python py/compute_backgrounds.py \
  --source cryostat_barrel --isotope Bi214 --sigma 25 --roi 50 \
  --bgdir LXeMC/results/backgrounds/cryostat/barrel/Bi214 \
  --fluxdir LXeMC/results/flux_tables/cryostat/barrel/Bi214

python py/compute_backgrounds.py \
  --source cryostat_barrel --isotope Tl208 --sigma 25 --roi 50 \
  --bgdir LXeMC/results/backgrounds/cryostat/barrel/Tl208 \
  --fluxdir LXeMC/results/flux_tables/cryostat/barrel/Tl208

python py/compute_backgrounds.py \
  --source cryostat_top --isotope Bi214 --sigma 25 --roi 50 \
  --bgdir LXeMC/results/backgrounds/cryostat/top/Bi214 \
  --fluxdir LXeMC/results/flux_tables/cryostat/top/Bi214

python py/compute_backgrounds.py \
  --source cryostat_top --isotope Tl208 --sigma 25 --roi 50 \
  --bgdir LXeMC/results/backgrounds/cryostat/top/Tl208 \
  --fluxdir LXeMC/results/flux_tables/cryostat/top/Tl208

python py/compute_backgrounds.py \
  --source cryostat_bottom --isotope Bi214 --sigma 25 --roi 50 \
  --bgdir LXeMC/results/backgrounds/cryostat/bottom/Bi214 \
  --fluxdir LXeMC/results/flux_tables/cryostat/bottom/Bi214

python py/compute_backgrounds.py \
  --source cryostat_bottom --isotope Tl208 --sigma 25 --roi 50 \
  --bgdir LXeMC/results/backgrounds/cryostat/bottom/Tl208 \
  --fluxdir LXeMC/results/flux_tables/cryostat/bottom/Tl208

# PMT
python py/compute_backgrounds.py \
  --source pmt_top --isotope Bi214 --sigma 25 --roi 50 \
  --bgdir LXeMC/results/backgrounds/pmt/top/Bi214 \
  --fluxdir LXeMC/results/flux_tables/pmt/top/Bi214

python py/compute_backgrounds.py \
  --source pmt_top --isotope Tl208 --sigma 25 --roi 50 \
  --bgdir LXeMC/results/backgrounds/pmt/top/Tl208 \
  --fluxdir LXeMC/results/flux_tables/pmt/top/Tl208

python py/compute_backgrounds.py \
  --source pmt_bottom --isotope Bi214 --sigma 25 --roi 50 \
  --bgdir LXeMC/results/backgrounds/pmt/bottom/Bi214 \
  --fluxdir LXeMC/results/flux_tables/pmt/bottom/Bi214

python py/compute_backgrounds.py \
  --source pmt_bottom --isotope Tl208 --sigma 25 --roi 50 \
  --bgdir LXeMC/results/backgrounds/pmt/bottom/Tl208 \
  --fluxdir LXeMC/results/flux_tables/pmt/bottom/Tl208

python py/compute_backgrounds.py \
  --source pmt_bottom_lxe --isotope Bi214 --sigma 25 --roi 50 \
  --bgdir LXeMC/results/backgrounds/pmt/bottom_lxe/Bi214 \
  --fluxdir LXeMC/results/flux_tables/pmt/bottom_lxe/Bi214

python py/compute_backgrounds.py \
  --source pmt_bottom_lxe --isotope Tl208 --sigma 25 --roi 50 \
  --bgdir LXeMC/results/backgrounds/pmt/bottom_lxe/Tl208 \
  --fluxdir LXeMC/results/flux_tables/pmt/bottom_lxe/Tl208

python py/compute_backgrounds.py \
  --source pmt_barrel --isotope Bi214 --sigma 25 --roi 50 \
  --bgdir LXeMC/results/backgrounds/pmt/barrel/Bi214 \
  --fluxdir LXeMC/results/flux_tables/pmt/barrel/Bi214

python py/compute_backgrounds.py \
  --source pmt_barrel --isotope Tl208 --sigma 25 --roi 50 \
  --bgdir LXeMC/results/backgrounds/pmt/barrel/Tl208 \
  --fluxdir LXeMC/results/flux_tables/pmt/barrel/Tl208
```

## Source summary

| Source | Isotopes | Components | Type |
|:-------|:---------|:-----------|:-----|
| cryostat_barrel | Bi214, Tl208 | OCV, MLI, ICV | Dense (KN + compound) |
| cryostat_top | Bi214, Tl208 | OCV, ICV | Dense (KN + compound) |
| cryostat_bottom | Bi214, Tl208 | OCV, ICV | Dense (KN + compound) |
| pmt_top | Bi214, Tl208 | PMTs, bases, structure | Transparent (merged PDisk) |
| pmt_bottom | Bi214, Tl208 | PMTs, bases, structure, R8778_dome | Transparent (merged PDisk) |
| pmt_bottom_lxe | Bi214, Tl208 | Same as pmt_bottom + passive LXe to cathode | Compound (transparent + LXe slab) |
| pmt_barrel | Bi214, Tl208 | cables, R8520 skin PMTs, R8778 lower-ring | Transparent (merged PCylShell) |

Total: 7 sources x 2 isotopes = 14 pipeline runs.
