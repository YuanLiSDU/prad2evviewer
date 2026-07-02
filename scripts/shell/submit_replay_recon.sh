#!/bin/bash
# submit_replay_recon.sh
# Submit one PRad-II replay pipeline job to JLab ifarm Slurm.
#
# The generated job runs:
#   1. prad2ana_replay_recon, with optional grouped hadd merging via -m
#   2. prad2ana_replay_filter over all recon ROOTs
#   3. prad2ana_live_charge over all filtered ROOTs together
#   4. prad2ana_quick_check over all recon ROOTs
#
# All outputs are written directly under <output_base>/prad_<RUN>/.
# The requested CPU count is shared by replay_recon -j, replay_filter -t,
# quick_check -j, and Slurm --cpus-per-task.

set -euo pipefail

# ---------------------------------------------------------------------------
# User defaults
# ---------------------------------------------------------------------------
# Edit this path for your ifarm installation. If PRAD2_SOFT is already set in
# the environment, for example by sourcing prad2_setup.sh, that value overrides
# this hard-coded default.
DEFAULT_PRAD2_SOFT="/w/hallb-scshelf2102/prad/your_path/prad2evviewer"

prompt_default() {
    local prompt="$1"
    local default="$2"
    local value
    read -rp "${prompt} [${default}]: " value
    if [[ -n "${value}" ]]; then
        printf '%s\n' "${value}"
    else
        printf '%s\n' "${default}"
    fi
}

to_abs_path() {
    local path="$1"
    local base="$2"
    if [[ "${path}" == /* ]]; then
        printf '%s\n' "${path}"
    else
        printf '%s\n' "${base}/${path#./}"
    fi
}

shell_quote() {
    printf "'%s'" "$(printf '%s' "$1" | sed "s/'/'\\\\''/g")"
}

SUBMIT_DIR="$(pwd)"
USER_NAME="$(whoami)"

PRAD2_BIN_WAS_SET=1
DEFAULT_CUTS_WAS_SET=1
[[ -z "${PRAD2_BIN:-}" ]] && PRAD2_BIN_WAS_SET=0
[[ -z "${DEFAULT_CUTS:-}" ]] && DEFAULT_CUTS_WAS_SET=0

PRAD2_SOFT="${PRAD2_SOFT:-${DEFAULT_PRAD2_SOFT}}"
PRAD2_BIN="${PRAD2_BIN:-${PRAD2_SOFT}/build/bin}"
CACHE_BASE="${CACHE_BASE:-/cache/clas12/rg-o/data}"
OUTPUT_BASE="${OUTPUT_BASE:-./}"
REPLAY_CORES="${REPLAY_CORES:-15}"
REPLAY_ZERO_SUPPRESS="${REPLAY_ZERO_SUPPRESS:-5}"
REPLAY_MAX_FILES="${REPLAY_MAX_FILES:-10000}"
REPLAY_MERGE_FILES="${REPLAY_MERGE_FILES:-62}"
DEFAULT_CUTS="${DEFAULT_CUTS:-${PRAD2_SOFT}/analysis/cuts/prad2_default.json}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-hallb}"
SLURM_PARTITION="${SLURM_PARTITION:-production}"
SLURM_TIME="${SLURM_TIME:-12:00:00}"
SLURM_MEM_PER_CPU="${SLURM_MEM_PER_CPU:-1500}"
ROOT_SETUP="${ROOT_SETUP:-}"

echo "Submit one PRad-II replay/recon Slurm job"
echo ""

read -rp "Enter run number (e.g. 024650): " RUN_NUMBER
if [[ -z "${RUN_NUMBER}" ]]; then
    echo "ERROR: run number cannot be empty."
    exit 1
fi

PRAD2_SOFT="$(prompt_default "Enter prad2evviewer directory" "${PRAD2_SOFT}")"
if [[ "${PRAD2_BIN_WAS_SET}" -eq 0 ]]; then
    PRAD2_BIN="${PRAD2_SOFT}/build/bin"
fi
if [[ "${DEFAULT_CUTS_WAS_SET}" -eq 0 ]]; then
    DEFAULT_CUTS="${PRAD2_SOFT}/analysis/cuts/prad2_default.json"
fi
PRAD2_BIN="$(prompt_default "Enter executable directory" "${PRAD2_BIN}")"
CACHE_BASE="$(prompt_default "Enter EVIO cache base directory" "${CACHE_BASE}")"
OUTPUT_BASE="$(prompt_default "Enter output base directory" "${OUTPUT_BASE}")"
REPLAY_CORES="$(prompt_default "Enter number of parallel jobs (-j)" "${REPLAY_CORES}")"
REPLAY_ZERO_SUPPRESS="$(prompt_default "Enter GEM zero suppression (-z)" "${REPLAY_ZERO_SUPPRESS}")"
REPLAY_MAX_FILES="$(prompt_default "Enter max number of files to process (-f)" "${REPLAY_MAX_FILES}")"
REPLAY_MERGE_FILES="$(prompt_default "Enter replay merge group size (-m, 0 disables)" "${REPLAY_MERGE_FILES}")"

echo "Enter cut JSON file for replay_filter (path to cuts.json, or 'default' to use ${DEFAULT_CUTS}):"
read -rp "Cut JSON [default]: " _INPUT
if [[ -z "${_INPUT}" || "${_INPUT}" == "default" ]]; then
    CUT_JSON="${DEFAULT_CUTS}"
else
    CUT_JSON="${_INPUT}"
fi

ROOT_SETUP="$(prompt_default "Enter ROOT setup script, or 'none' to rely on submitted environment" "${ROOT_SETUP:-none}")"
if [[ "${ROOT_SETUP}" == "none" ]]; then
    ROOT_SETUP=""
fi

SLURM_ACCOUNT="$(prompt_default "Enter Slurm account" "${SLURM_ACCOUNT}")"
SLURM_PARTITION="$(prompt_default "Enter Slurm partition (production or priority)" "${SLURM_PARTITION}")"
SLURM_TIME="$(prompt_default "Enter Slurm time limit" "${SLURM_TIME}")"
SLURM_MEM_PER_CPU="$(prompt_default "Enter Slurm mem-per-cpu MB" "${SLURM_MEM_PER_CPU}")"

PRAD2_SOFT="$(to_abs_path "${PRAD2_SOFT}" "${SUBMIT_DIR}")"
PRAD2_BIN="$(to_abs_path "${PRAD2_BIN}" "${SUBMIT_DIR}")"
CUT_JSON="$(to_abs_path "${CUT_JSON}" "${SUBMIT_DIR}")"
if [[ -n "${ROOT_SETUP}" ]]; then
    ROOT_SETUP="$(to_abs_path "${ROOT_SETUP}" "${SUBMIT_DIR}")"
fi
OUTPUT_BASE="$(to_abs_path "${OUTPUT_BASE}" "${SUBMIT_DIR}")"
OUT_DIR="${OUTPUT_BASE}/prad_${RUN_NUMBER}"
RUN_DIR="${CACHE_BASE}/prad_${RUN_NUMBER}"
REPORT="${OUT_DIR}/prad_${RUN_NUMBER}_filter_report.json"
SLOW_ROOT="${OUT_DIR}/prad_${RUN_NUMBER}_epics.root"
QC_OUTPUT="${OUT_DIR}/prad_${RUN_NUMBER}_quick_check.root"
LC_JSON="${OUT_DIR}/prad_${RUN_NUMBER}_live_charge.json"
SBATCH_SCRIPT="${OUT_DIR}/replay_recon_${RUN_NUMBER}.sbatch"

echo ""
echo "Checking input directory: ${RUN_DIR}"
if [[ ! -d "${RUN_DIR}" ]]; then
    echo "ERROR: directory not found: ${RUN_DIR}"
    exit 1
fi

FILE_COUNT=$(find "${RUN_DIR}" -maxdepth 1 -type f | wc -l)
if [[ "${FILE_COUNT}" -eq 0 ]]; then
    echo "No files found in ${RUN_DIR}."
    echo "The data may still be on tape (MSS)."
    read -rp "Enter your email address for jcache notification: " USER_EMAIL
    if [[ -z "${USER_EMAIL}" ]]; then
        echo "ERROR: email cannot be empty."
        exit 1
    fi
    MSS_DIR="/mss/clas12/rg-o/data/prad_${RUN_NUMBER}"
    echo "Submitting: jcache get ${MSS_DIR}/* -e ${USER_EMAIL}"
    jcache get "${MSS_DIR}"/* -e "${USER_EMAIL}"
    echo "jcache request submitted. Re-run this submit script after staging finishes."
    exit 0
fi

mkdir -p "${OUT_DIR}"

for exe in prad2ana_replay_recon; do
    if [[ ! -x "${PRAD2_BIN}/${exe}" ]]; then
        echo "ERROR: executable not found: ${PRAD2_BIN}/${exe}"
        exit 1
    fi
done

cat > "${SBATCH_SCRIPT}" <<EOF
#!/bin/bash
#SBATCH --job-name=recon${RUN_NUMBER}
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --output=${OUT_DIR}/replay_recon_${RUN_NUMBER}-%j.out
#SBATCH --error=${OUT_DIR}/replay_recon_${RUN_NUMBER}-%j.err
#SBATCH --mail-user=${USER_NAME}@jlab.org
#SBATCH --time=${SLURM_TIME}
#SBATCH --mem-per-cpu=${SLURM_MEM_PER_CPU}
#SBATCH --cpus-per-task=${REPLAY_CORES}

set -euo pipefail

RUN_NUMBER=$(shell_quote "${RUN_NUMBER}")
RUN_DIR=$(shell_quote "${RUN_DIR}")
OUT_DIR=$(shell_quote "${OUT_DIR}")
PRAD2_BIN=$(shell_quote "${PRAD2_BIN}")
CUT_JSON=$(shell_quote "${CUT_JSON}")
REPLAY_CORES=$(shell_quote "${REPLAY_CORES}")
REPLAY_ZERO_SUPPRESS=$(shell_quote "${REPLAY_ZERO_SUPPRESS}")
REPLAY_MAX_FILES=$(shell_quote "${REPLAY_MAX_FILES}")
REPLAY_MERGE_FILES=$(shell_quote "${REPLAY_MERGE_FILES}")
REPORT=$(shell_quote "${REPORT}")
SLOW_ROOT=$(shell_quote "${SLOW_ROOT}")
QC_OUTPUT=$(shell_quote "${QC_OUTPUT}")
LC_JSON=$(shell_quote "${LC_JSON}")
ROOT_SETUP=$(shell_quote "${ROOT_SETUP}")

echo "Host: \$(hostname)"
echo "Work dir: \$(pwd)"
echo "Run: \${RUN_NUMBER}"
echo "Input: \${RUN_DIR}"
echo "Output: \${OUT_DIR}"
date

if [[ -n "\${ROOT_SETUP}" ]]; then
    echo "Sourcing ROOT setup: \${ROOT_SETUP}"
    source "\${ROOT_SETUP}"
fi

if ! command -v root >/dev/null 2>&1; then
    echo "ERROR: root is not found in PATH. Set ROOT_SETUP or submit from a ROOT-ready environment."
    exit 1
fi
echo "ROOT is available: \$(root-config --version 2>/dev/null || root --version 2>&1 | head -1)"

if [[ ! -d "\${RUN_DIR}" ]]; then
    echo "ERROR: input directory not found: \${RUN_DIR}"
    exit 1
fi
FILE_COUNT=\$(find "\${RUN_DIR}" -maxdepth 1 -type f | wc -l)
if [[ "\${FILE_COUNT}" -eq 0 ]]; then
    echo "ERROR: no input files found in \${RUN_DIR}. Stage them with jcache first."
    exit 1
fi

mkdir -p "\${OUT_DIR}"
cd "\${OUT_DIR}"
echo "Job work dir: \$(pwd)"

REPLAY_CMD="\${PRAD2_BIN}/prad2ana_replay_recon"
QUICK_CHECK_CMD="\${PRAD2_BIN}/prad2ana_quick_check"
FILTER_CMD="\${PRAD2_BIN}/prad2ana_replay_filter"
LIVE_CHARGE_CMD="\${PRAD2_BIN}/prad2ana_live_charge"

if [[ ! -x "\${REPLAY_CMD}" ]]; then
    echo "ERROR: executable not found: \${REPLAY_CMD}"
    exit 1
fi

echo ""
echo "Starting replay..."
echo "\${REPLAY_CMD} \${RUN_DIR} -o \${OUT_DIR} -j \${REPLAY_CORES} -z \${REPLAY_ZERO_SUPPRESS} -f \${REPLAY_MAX_FILES} -m \${REPLAY_MERGE_FILES}"
"\${REPLAY_CMD}" "\${RUN_DIR}" -o "\${OUT_DIR}" -j "\${REPLAY_CORES}" -z "\${REPLAY_ZERO_SUPPRESS}" -f "\${REPLAY_MAX_FILES}" -m "\${REPLAY_MERGE_FILES}"

mapfile -t RECON_INPUTS < <(find "\${OUT_DIR}" -maxdepth 1 -type f -name "prad_\${RUN_NUMBER}_recon_*.root" | sort)
if [[ "\${#RECON_INPUTS[@]}" -eq 0 ]]; then
    mapfile -t RECON_INPUTS < <(find "\${OUT_DIR}" -maxdepth 1 -type f -name "prad_\${RUN_NUMBER}.*_recon.root" | sort)
fi
if [[ "\${#RECON_INPUTS[@]}" -eq 0 ]]; then
    echo "ERROR: no reconstructed ROOT files found in \${OUT_DIR}."
    exit 1
fi
echo "Downstream input ROOT file(s): \${#RECON_INPUTS[@]}"

LC_INPUTS=()
if [[ ! -x "\${FILTER_CMD}" ]]; then
    echo "WARNING: executable not found: \${FILTER_CMD}, skipping replay filter."
elif [[ ! -f "\${CUT_JSON}" ]]; then
    echo "WARNING: cut JSON not found: \${CUT_JSON}, skipping replay filter."
else
    echo ""
    echo "Running replay filter..."
    if [[ "\${#RECON_INPUTS[@]}" -gt 1 ]]; then
        FILTER_OUT="\${OUT_DIR}"
    else
        _BASE="\$(basename "\${RECON_INPUTS[0]}")"
        if [[ "\${_BASE}" =~ ^(.+)_recon_([^/]+)\.root$ ]]; then
            FILTER_OUT="\${OUT_DIR}/\${BASH_REMATCH[1]}_filter_\${BASH_REMATCH[2]}.root"
        elif [[ "\${_BASE}" =~ ^(.+)\.evio\.([0-9]+)_recon\.root$ ]]; then
            FILTER_OUT="\${OUT_DIR}/\${BASH_REMATCH[1]}_filter_\${BASH_REMATCH[2]}.root"
        elif [[ "\${_BASE}" =~ ^(.+)\.([0-9]+)_recon\.root$ ]]; then
            FILTER_OUT="\${OUT_DIR}/\${BASH_REMATCH[1]}_filter_\${BASH_REMATCH[2]}.root"
        else
            FILTER_OUT="\${OUT_DIR}/\${_BASE%.root}_filter.root"
        fi
    fi
    "\${FILTER_CMD}" "\${RECON_INPUTS[@]}" -o "\${FILTER_OUT}" -c "\${CUT_JSON}" -j "\${REPORT}" -t "\${REPLAY_CORES}"
    mapfile -t FILTER_ROOTS < <(find "\${OUT_DIR}" -maxdepth 1 -type f -name "prad_\${RUN_NUMBER}_filter*.root" | sort)
    if [[ "\${#FILTER_ROOTS[@]}" -gt 0 ]]; then
        LC_INPUTS=("\${FILTER_ROOTS[@]}")
    else
        echo "ERROR: no filtered ROOT files found in \${OUT_DIR}; live_charge requires prad_\${RUN_NUMBER}_filter*.root inputs."
        exit 1
    fi
fi

if [[ ! -x "\${LIVE_CHARGE_CMD}" ]]; then
    echo "WARNING: executable not found: \${LIVE_CHARGE_CMD}, skipping live charge."
else
    if [[ "\${#LC_INPUTS[@]}" -eq 0 ]]; then
        echo "ERROR: live_charge requires filtered ROOT files, but none were selected."
        exit 1
    fi
    echo ""
    echo "Running live charge calculation..."
    "\${LIVE_CHARGE_CMD}" "\${LC_INPUTS[@]}" -j "\${LC_JSON}"
fi

if [[ ! -x "\${QUICK_CHECK_CMD}" ]]; then
    echo "WARNING: executable not found: \${QUICK_CHECK_CMD}, skipping quick check."
else
    echo ""
    echo "Running quick check..."
    "\${QUICK_CHECK_CMD}" "\${RECON_INPUTS[@]}" -o "\${QC_OUTPUT}" -j "\${REPLAY_CORES}"
fi

echo ""
echo "Replay/recon job finished."
echo "Recon inputs: \${#RECON_INPUTS[@]}"
echo "Output dir: \${OUT_DIR}"
echo "Filter report: \${REPORT}"
echo "Slow ROOT: \${SLOW_ROOT}"
echo "Live charge JSON: \${LC_JSON}"
echo "Quick check ROOT: \${QC_OUTPUT}"
date
EOF

chmod +x "${SBATCH_SCRIPT}"

echo ""
echo "Prepared Slurm script: ${SBATCH_SCRIPT}"
echo "Submitting job..."
sbatch "${SBATCH_SCRIPT}"
