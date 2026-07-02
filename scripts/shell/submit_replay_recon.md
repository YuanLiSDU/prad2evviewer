# submit_replay_recon.sh

Submit one PRad-II replay pipeline job to the JLab ifarm Slurm system.

This script is the batch-job version of `replay_recon.sh`. It keeps the same
interactive style, but instead of running the full replay pipeline directly on
the login node, it generates one Slurm job for one run and submits it with
`sbatch`.

Pipeline:

```text
cache check -> optional jcache staging -> sbatch
    -> replay recon -> filter -> live charge -> quick check
```

The generated job keeps all products directly under
`<OUTPUT_BASE>/prad_<RUN>/`.  `prad2ana_replay_recon` does grouped `hadd`
merging internally by default, with 80 split ROOT files per merged recon
ROOT.  The same CPU count is used for replay, filter, quick check, and Slurm
`--cpus-per-task`.

---

## Usage

Run this script from the ifarm login node:

```bash
cd /path/to/prad2evviewer
chmod +x scripts/shell/submit_replay_recon.sh
bash scripts/shell/submit_replay_recon.sh
```

You can also run it from another working directory, but then make sure the
prompted `prad2evviewer directory` points to the correct installation.

> Important: run with `bash`. Do not use `source`.

Before first use, edit the hard-coded default near the top of the script:

```bash
DEFAULT_PRAD2_SOFT="/w/hallb-scshelf2102/prad/your_path/prad2evviewer"
```

If `PRAD2_SOFT` is already set in the environment, for example after sourcing a
local PRad-II setup script, the environment value overrides this hard-coded
default.

---

## Batch Submit Multiple Runs

For submitting multiple one-run jobs, use:

```bash
chmod +x scripts/shell/submit_replay_recon_m.sh
bash scripts/shell/submit_replay_recon_m.sh
```

This wrapper still submits one Slurm job per run. It simply calls
`submit_replay_recon.sh` repeatedly with the same replay and Slurm parameters.

There are two run input modes.

Range mode:

```text
Enter start run number, or press Enter to input a run list: 024650
Enter end run number [024650]: 024655
```

This submits all runs from `024650` through `024655`, inclusive.

List mode:

```text
Enter start run number, or press Enter to input a run list:
Enter run list separated by spaces or commas: 024650 024652,024660
```

Before submitting, the wrapper checks the cache directory for every requested
run. Runs already in cache are submitted immediately. Runs missing from cache
are skipped first; after all ready jobs have been submitted, the wrapper asks
for an email address and submits `jcache get` requests for the missing runs.

At the end, it prints the run list that needs to wait for staging and be
submitted again later:

```text
Wait for staging, then re-run this script for these run(s):
024651 024653
```

---

## What This Script Does

The submit script itself runs only the lightweight preparation steps:

1. Ask for the run number and replay/job parameters.
2. Check the EVIO input directory in `/cache`.
3. If the cache directory is empty, submit a `jcache get` tape-staging request
   and exit.
4. Create the output directory.
5. Generate a Slurm batch script:
   ```text
   <OUTPUT_BASE>/prad_<RUN>/replay_recon_<RUN>.sbatch
   ```
6. Submit that batch script with `sbatch`.

The Slurm job then runs the full replay/recon pipeline on a compute node.

---

## Interactive Parameters

Press Enter to accept the default value shown in brackets.

| Prompt | Description | Default |
|--------|-------------|---------|
| `Enter run number` | 6-digit PRad-II run number, e.g. `024650` | required |
| `Enter prad2evviewer directory` | Source/installation directory | hard-coded `DEFAULT_PRAD2_SOFT`, or environment `PRAD2_SOFT` |
| `Enter executable directory` | Directory containing `prad2ana_*` executables | `<PRAD2_SOFT>/build/bin` |
| `Enter EVIO cache base directory` | Base directory for raw EVIO data | `/cache/clas12/rg-o/data` |
| `Enter output base directory` | Base directory for all output files | `./` |
| `Enter number of parallel jobs (-j)` | Shared CPU count for `replay_recon -j`, `replay_filter -t`, `quick_check -j`, and Slurm `--cpus-per-task` | `50` |
| `Enter GEM zero suppression (-z)` | GEM zero-suppression sigma threshold passed to replay | `5` |
| `Enter max number of files to process (-f)` | Maximum number of EVIO sub-files to process | `10000` |
| `Enter replay merge group size (-m, 0 disables)` | Number of split recon ROOT files per merged output | `80` |
| `Cut JSON [default]` | Cut config for `prad2ana_replay_filter` | `<PRAD2_SOFT>/analysis/cuts/prad2_default.json` |
| `Enter ROOT setup script` | Optional ROOT setup script to source inside the job | `none` |
| `Enter Slurm account` | Slurm account | `hallb` |
| `Enter Slurm partition (production or priority)` | Slurm partition | `production` |
| `Enter Slurm time limit` | Job wall time | `12:00:00` |
| `Enter Slurm mem-per-cpu MB` | Memory per CPU in MB | `1500` |

All relative paths entered at the prompts are resolved relative to the directory
where `submit_replay_recon.sh` is launched.

---

## Environment Variable Overrides

You can override defaults before running the script:

```bash
PRAD2_SOFT=/w/halla-scshelf2102/your/path/prad2evviewer \
OUTPUT_BASE=/volatile/halla/your/output/recon \
REPLAY_CORES=32 \
bash scripts/shell/submit_replay_recon.sh
```

For example, if your setup script exports `PRAD2_SOFT`, you can do:

```bash
source scripts/shell/prad2_setup.sh
bash scripts/shell/submit_replay_recon.sh
```

| Variable | Description |
|----------|-------------|
| `PRAD2_SOFT` | Root directory of the prad2evviewer source/installation |
| `PRAD2_BIN` | Executable directory |
| `CACHE_BASE` | Base directory for EVIO input data |
| `OUTPUT_BASE` | Base directory for output files |
| `REPLAY_CORES` | Shared CPU count for replay, filter, quick check, and Slurm `--cpus-per-task` |
| `REPLAY_ZERO_SUPPRESS` | GEM zero-suppression sigma threshold |
| `REPLAY_MAX_FILES` | Maximum number of EVIO sub-files |
| `REPLAY_MERGE_FILES` | Number of split recon ROOT files per merged output; `0` disables merging |
| `DEFAULT_CUTS` | Default replay filter cut JSON |
| `ROOT_SETUP` | ROOT setup script sourced inside the Slurm job |
| `SLURM_ACCOUNT` | Slurm account |
| `SLURM_PARTITION` | Slurm partition |
| `SLURM_TIME` | Slurm wall time |
| `SLURM_MEM_PER_CPU` | Slurm memory per CPU in MB |

If `PRAD2_BIN` or `DEFAULT_CUTS` are not set explicitly, they are updated
automatically when you change `PRAD2_SOFT` at the prompt.

---

## Input Files

For run `024650`, the input directory is:

```text
/cache/clas12/rg-o/data/prad_024650
```

The script checks for regular files in this directory before submitting the
Slurm job.

If the directory exists but is empty, the script assumes the data may still be
on tape and asks for an email address. It then runs:

```bash
jcache get /mss/clas12/rg-o/data/prad_<RUN>/* -e <your@email>
```

After you receive the staging notification, run `submit_replay_recon.sh` again.

---

## Output Files

All files for one run are written to:

```text
<OUTPUT_BASE>/prad_<RUN>/
```

Expected outputs:

| File | Description |
|------|-------------|
| `replay_recon_<RUN>.sbatch` | Generated Slurm batch script |
| `replay_recon_<RUN>-<JOBID>.out` | Slurm stdout |
| `replay_recon_<RUN>-<JOBID>.err` | Slurm stderr |
| `prad_<RUN>.evio.XXXXX_recon.root` | Per-split recon ROOT output from replay |
| `prad_<RUN>_recon_000.root`, `prad_<RUN>_recon_001.root`, ... | Merged recon ROOT files, grouped by `REPLAY_MERGE_FILES` |
| `prad_<RUN>_filter_000.root`, `prad_<RUN>_filter_001.root`, ... | Filtered ROOT files, one per recon input |
| `prad_<RUN>_epics.root` | Run-level slow-control ROOT containing only `scalers`, `epics`, and `runinfo` |
| `prad_<RUN>_filter_report.json` | Run-level replay filter report |
| `prad_<RUN>_live_charge.json` | Live charge result from all filtered ROOT files |
| `prad_<RUN>_quick_check.root` | Quick-check histograms from all recon ROOT files |

If `REPLAY_MERGE_FILES=0`, the downstream filter and quick-check steps use the
per-split `prad_<RUN>.evio.XXXXX_recon.root` files directly.

---

## Slurm Job Details

The generated batch script uses:

```bash
#SBATCH --job-name=recon<RUN>
#SBATCH --account=<SLURM_ACCOUNT>
#SBATCH --partition=<SLURM_PARTITION>
#SBATCH --output=<OUT_DIR>/replay_recon_<RUN>-%j.out
#SBATCH --error=<OUT_DIR>/replay_recon_<RUN>-%j.err
#SBATCH --mail-user=<username>@jlab.org
#SBATCH --time=<SLURM_TIME>
#SBATCH --mem-per-cpu=<SLURM_MEM_PER_CPU>
#SBATCH --cpus-per-task=<REPLAY_CORES>
```

Inside the job, the script changes into the run output directory before running
the replay pipeline. This helps keep logs and automatically named outputs close
to the final ROOT/JSON files.

---

## Pipeline Inside the Job

The submitted job runs these steps:

```text
1. Optionally source ROOT_SETUP
2. Check root
3. Check the cache input directory again
4. Run prad2ana_replay_recon with -j REPLAY_CORES and -m REPLAY_MERGE_FILES
5. Select merged recon ROOTs, or per-split recon ROOTs when -m 0
6. Run prad2ana_replay_filter on all selected recon ROOTs with -t REPLAY_CORES
7. Run prad2ana_live_charge on all filtered ROOTs together
8. Run prad2ana_quick_check on all selected recon ROOTs with -j REPLAY_CORES
```

`live_charge` does not fall back to unfiltered recon files. It requires the
filtered ROOT files produced by `prad2ana_replay_filter`.

---

## Example Session

```text
$ bash scripts/shell/submit_replay_recon.sh
Submit one PRad-II replay/recon Slurm job

Enter run number (e.g. 024650): 024650
Enter prad2evviewer directory [/w/halla-scshelf2102/user/prad2evviewer]:
Enter executable directory [/w/halla-scshelf2102/user/prad2evviewer/build/bin]:
Enter EVIO cache base directory [/cache/clas12/rg-o/data]:
Enter output base directory [./]: /volatile/halla/user/recon
Enter number of parallel jobs (-j) [50]:
Enter GEM zero suppression (-z) [5]:
Enter max number of files to process (-f) [10000]:
Enter replay merge group size (-m, 0 disables) [80]:
Enter cut JSON file for replay_filter (...):
Cut JSON [default]:
Enter ROOT setup script, or 'none' to rely on submitted environment [none]:
Enter Slurm account [hallb]:
Enter Slurm partition (production or priority) [production]:
Enter Slurm time limit [12:00:00]:
Enter Slurm mem-per-cpu MB [1500]:

Checking input directory: /cache/clas12/rg-o/data/prad_024650

Prepared Slurm script: /volatile/halla/user/recon/prad_024650/replay_recon_024650.sbatch
Submitting job...
Submitted batch job 12345678
```

Check job status with:

```bash
squeue -u $USER
```

After the job starts, follow the Slurm stdout log:

```bash
tail -f /volatile/halla/user/recon/prad_024650/replay_recon_024650-12345678.out
```

---

## Notes and Caveats

- This script submits exactly one run as one Slurm job.
- It is intended for JLab ifarm and assumes `sbatch`, `jcache`, ROOT, and the
  PRad-II replay executables are available there.
- The script uses bash both for submission and for the generated batch job.
- If ROOT is not available by default on compute nodes, provide a valid
  `ROOT_SETUP`, for example:
  ```bash
  /path/to/root/bin/thisroot.sh
  ```
- The generated `.sbatch` file is kept in the output directory so the exact job
  configuration can be inspected or resubmitted manually.
