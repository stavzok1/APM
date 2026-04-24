# Reproducing your Python + tools environment on a Linux cluster

This guide assumes the cluster runs **Linux** (like your WSL/Ubuntu box) and you already clone the repo with Git and sync `data/` separately (see `docs/GIT_GUIDE.md`).

You care about three layers:

1. **OS packages** (`apt`, or the site’s **environment modules** / **Spack** — use what your cluster documents).
2. **APM’s Python env** — usually a **venv** at `APM/.venv` with everything `pip install`’d.
3. **Separate tool envs** — e.g. **VEP** (`vep_env`), R, Java, etc., which are *not* the same as `.venv`.

---

## 1. On your laptop (once): capture what you actually use

### A. APM project virtualenv (`.venv`)

From the repo root, with the env **activated**:

```bash
cd /home/stavz/masters/gdc/APM
source .venv/bin/activate
python --version          # note this (e.g. 3.12) — match on cluster if you can
pip install -U pip wheel
pip freeze > requirements-lock.txt
deactivate
```

- **`pip freeze`** writes **exact versions** of every installed package. That file is what you use on the cluster for a **reproducible** install.
- Put **`requirements-lock.txt`** in Git (or rsync it) so the cluster always gets the same pins.

Optional (human-readable top-level deps only — harder to reproduce exactly):

```bash
pip list --format=freeze > requirements-lock.txt   # same as pip freeze on modern pip
```

### B. Separate `vep_env` (or any second environment)

**If it is a Conda/Mamba env:**

```bash
conda activate vep_env
conda env export --no-builds > vep_env.yml
# or, from any env:
conda list -n vep_env --explicit > vep_env_explicit.txt
```

`--no-builds` keeps the YAML smaller and more portable across machines.

**If it is another `venv`:**

```bash
source /path/to/vep_env/bin/activate
pip freeze > vep_requirements-lock.txt
deactivate
```

Keep those files **next to the repo** or **in `docs/`** if they contain no secrets (they usually do not).

### C. System / “tools” you rely on (optional but useful)

On Ubuntu/WSL, list manually what you installed for bioinformatics / building wheels, e.g.:

```bash
# examples — adjust to what you actually use
dpkg -l | grep -E 'zlib|bz2|lzma|ssl|curl|git|build-essential' 
```

On the cluster, you often **do not** use `apt` yourself; instead you run:

```bash
module avail
module load python/3.12  # example name — use your site’s module names
```

Document for yourself: **Python version**, **GCC**, **zlib**, **git**, **Java** (for Picard/GATK-style tools), **Perl** (for VEP), etc.

---

## 2. On the cluster: OS + Python version

### Match Python major.minor when possible

If the laptop used **Python 3.12**, prefer **3.12** on the cluster so wheels and behaviour match.

Check:

```bash
python3 --version
which python3
```

If the cluster only has 3.10, you can still try installing from the same `requirements-lock.txt`, but **some pins may not have wheels** for 3.10; you may need a relaxed file or a container.

### OS packages (when you *are* allowed `sudo apt`)

Typical build deps for `pip` compiling extensions (install **before** `pip install -r`):

```bash
sudo apt update
sudo apt install -y \
  build-essential \
  python3-dev \
  python3-venv \
  git \
  curl \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  libssl-dev \
  libffi-dev
```

Your site may use **`module load`** instead; follow HPC docs.

---

## 3. On the cluster: recreate APM’s `.venv`

From the **repository root** (after `git clone` + `cd APM`):

```bash
cd /path/to/APM

# Create a fresh venv (do not copy .venv/ from laptop — Linux paths + compiled bits differ)
python3 -m venv .venv
source .venv/bin/activate

python -m pip install -U pip wheel
pip install -r requirements-lock.txt

# Quick sanity check
python -c "import pandas, numpy; print('ok', pandas.__version__)"
python -c "import pyranges" 2>/dev/null || echo "pyranges missing (optional for some tests)"
deactivate
```

**Why not `rsync` the whole `.venv` from home?**  
Cross-machine venvs often break (absolute shebangs, CPU-specific wheels, different glibc). **Recreate with `pip install -r`** on the cluster.

---

## 4. On the cluster: recreate `vep_env` (Conda example)

If you exported `vep_env.yml`:

```bash
# install Miniforge/Mambaforge once per user if you do not have conda
conda env create -f vep_env.yml
conda activate vep_env
# verify
vep --help || which vep
```

If the cluster provides **centrally installed VEP**, prefer that and only mirror **plugin / cache paths** in your job script instead of maintaining your own conda VEP.

---

## 5. Running the pipeline on the cluster

Always activate the APM venv first:

```bash
cd /path/to/APM
source .venv/bin/activate
export APM_WORKING_DIR=/path/to/your/data   # if data is not under <repo>/data
python -m pipeline.main   # or your entrypoint / Slurm script
```

Use the same **`APM_*`** variables as on your laptop if paths differ (see `pipeline/config.py` / comments in `PathConfig`).

---

## 6. Optional: faster / stricter installs

- **`uv`** (`pip install uv` then `uv pip sync requirements-lock.txt`) — often faster; still respects pinned versions if the file is a freeze.
- **Containers** (Singularity/Apptainer) — best when the site cannot give you the right Python or system libs; you build an image once with all deps.

---

## 7. Checklist

| Step | Laptop | Cluster |
|------|--------|---------|
| Lock APM Python deps | `pip freeze > requirements-lock.txt` | `pip install -r requirements-lock.txt` in new `.venv` |
| Lock VEP / other env | `conda env export` or second `pip freeze` | `conda env create` or second venv |
| Match Python | note `python --version` | same or adjust lock file |
| System tools | note `apt` / modules you use | `apt` or `module load` equivalents |
| Data | — | `rsync` or shared filesystem + `APM_WORKING_DIR` |

---

## 8. If `pip install` fails on the cluster

- **Missing compiler / headers** — install `build-essential`, `python3-dev`, `zlib1g-dev`, etc. (see §2).
- **Old OpenSSL / glibc** — try an older wheel, loosen one package in a **copy** of the lock file, or use a **container** / newer `module load python`.
- **Private Git URLs in freeze** — replace with PyPI versions or vendor wheels.

After you generate **`requirements-lock.txt`** on your laptop, add and commit it so the cluster always installs the same stack as your dev machine.
