# replay_recon.sh

A one-shot shell script that automates the full replay pipeline for a single PRad-II run on JLab ifarm:
replay → merge → quick check → filter → live charge.

---

## Usage

```bash
cd /path/to/working/dir
chmod +x /path/to/prad2evviewer/scripts/shell/replay_recon.sh
bash /path/to/prad2evviewer/scripts/shell/replay_recon.sh
```

> **Important**: Run with `bash` directly. Do **not** use `source` — ifarm's default shell is `tcsh` and the bash-specific syntax (`${VAR:-default}`) will cause errors.

---

## Prerequisites

- ROOT environment is set up (`root` command available):
  ```bash
  source /path/to/root/bin/thisroot.sh
  ```
- `hadd` and `jcache` are available (standard on JLab ifarm)
- `prad2evviewer` has been compiled (`build/bin/` contains the executables)

---

## Interactive Parameters

The script prompts for the following inputs at startup. Press Enter to accept the default value shown in brackets.

| Prompt | Description | Default |
|--------|-------------|---------|
| `Enter run number` | 6-digit run number, e.g. `024650` | *(required)* |
| `Enter output base directory` | Root directory for all output files | `./` |
| `Enter number of parallel jobs (-j)` | Number of replay threads | `50` |
| `Enter GEM zero suppression (-z)` | zlib compression level | `5` |
| `Enter max number of files to process (-f)` | Maximum number of EVIO sub-files to process | `10000` |
| `Cut JSON [default]` | Path to the cut config file for `replay_filter`; type `default` to use the built-in template | `prad2evviewer/analysis/cuts/prad2_default.json` |

All defaults can also be overridden via environment variables before running:

```bash
PRAD2_SOFT=/data/soft/prad2evviewer OUTPUT_BASE=/data/recon bash replay_recon.sh
```

| Variable | Description |
|----------|-------------|
| `PRAD2_SOFT` | Root directory of the prad2evviewer source/installation |
| `PRAD2_BIN` | Executable directory (default: `$PRAD2_SOFT/build/bin`) |
| `CACHE_BASE` | Base directory for EVIO input data |
| `OUTPUT_BASE` | Base directory for all output files |
| `REPLAY_CORES` | Number of parallel replay threads |
| `REPLAY_ZERO_SUPPRESS` | zlib compression level |
| `REPLAY_MAX_FILES` | Maximum number of EVIO sub-files to process |
| `DEFAULT_CUTS` | Path to the default cut JSON file |

---

## Input Files

| Content | Path |
|---------|------|
| EVIO raw data | `/cache/clas12/rg-o/data/prad_<RUN>/` |
| Cut configuration | User-specified, or `prad2evviewer/analysis/cuts/prad2_default.json` |

---

## Output Files

All output files are written to `<OUTPUT_BASE>/prad_<RUN>/`:

| File | Description |
|------|-------------|
| `prad_<RUN>.00**_recon.root` | Per-sub-file ROOT output from replay (deleted after merging) |
| `prad_<RUN>_recon.root` | Merged recon ROOT file |
| `prad_<RUN>_quick_check.root` | Quick-check histograms (auto-named by `prad2ana_quick_check`) |
| `prad_<RUN>_filtered.root` | Physics events passing the slow-control cuts |
| `prad_<RUN>_filter_report.json` | Per-checkpoint filter verdict from `replay_filter` |
| `prad_<RUN>_live_charge.json` | Live charge result in nC from `live_charge` |

---

## Pipeline Steps

```
0. Check ROOT environment
1. Collect parameters interactively
2. Count files in /cache input directory
   └─ If empty → submit jcache tape-staging request and exit
3. Create output directory
4. prad2ana_replay_recon    (multi-threaded replay)
5. hadd                     (merge per-sub-file ROOT outputs)
6. Delete per-sub-file ROOT outputs
7. prad2ana_quick_check     (fast quality-check histograms)
8. prad2ana_replay_filter   (apply slow-control cuts)
9. prad2ana_live_charge     (compute live charge)
```

---

## Data Not Yet in Cache (jcache Tape Staging)

If `/cache/clas12/rg-o/data/prad_<RUN>/` is empty, the data is likely still on tape.
The script will prompt for your email address and submit a staging request:

```bash
jcache get /mss/clas12/rg-o/data/prad_<RUN>/* -e <your@email>
```

You will receive an email notification when the files have been moved to cache.
Re-run the script after receiving the notification.

---

## Example Session

```
$ bash replay_recon.sh (or just "./replay_recon.sh")
ROOT is available: 6.32.02
Enter run number (e.g. 024650): 024650
Enter output base directory [./]: /home/liyuan/PRad2Analysis/data/recon
Enter number of parallel jobs (-j) [50]:
Enter GEM zero suppression (-z) [5]:
Enter max number of files to process (-f) [10000]:
Enter cut JSON file for replay_filter (...):
Cut JSON [default]:

Checking input directory: /cache/clas12/rg-o/data/prad_024650
Found 42 file(s) in /cache/clas12/rg-o/data/prad_024650:
...
Output directory: /home/liyuan/PRad2Analysis/data/recon/prad_024650
Starting replay...
...
Merging 42 ROOT file(s) into: .../prad_024650/prad_024650_recon.root
...
Running quick check...
Running replay filter...
Running live charge calculation...
Live charge JSON: .../prad_024650/prad_024650_live_charge.json
```
