#!/bin/bash
# submit_replay_recon_m.sh
# Submit multiple single-run PRad replay pipeline Slurm jobs.
#
# Prompts once for shared settings, then feeds each cached run to
# submit_replay_recon.sh.  The CPU count is entered once and reused by
# replay_recon -j, replay_filter -t, quick_check -j, and Slurm
# --cpus-per-task in each generated job.  Each run writes all products
# directly under <output_base>/prad_<RUN>/.

set -euo pipefail

# ---------------------------------------------------------------------------
# User defaults
# ---------------------------------------------------------------------------
# Keep this in sync with submit_replay_recon.sh, or export PRAD2_SOFT before
# running this script to override the hard-coded default.
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

join_by_space() {
    local first=1
    local item
    for item in "$@"; do
        if [[ "${first}" -eq 1 ]]; then
            printf '%s' "${item}"
            first=0
        else
            printf ' %s' "${item}"
        fi
    done
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUBMIT_ONE="${SCRIPT_DIR}/submit_replay_recon.sh"
SUBMIT_DIR="$(pwd)"

if [[ ! -f "${SUBMIT_ONE}" ]]; then
    echo "ERROR: single-run submit script not found: ${SUBMIT_ONE}"
    exit 1
fi

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

echo "Submit multiple PRad replay/recon Slurm jobs"
echo ""

RUNS=()
read -rp "Enter start run number, or press Enter to input a run list: " START_RUN
if [[ -z "${START_RUN}" ]]; then
    read -rp "Enter run list separated by spaces or commas: " RUN_LIST
    if [[ -z "${RUN_LIST}" ]]; then
        echo "ERROR: run list cannot be empty."
        exit 1
    fi
    RUN_LIST="${RUN_LIST//,/ }"
    for run in ${RUN_LIST}; do
        if [[ ! "${run}" =~ ^[0-9]+$ ]]; then
            echo "ERROR: invalid run number: ${run}"
            exit 1
        fi
        RUNS+=("${run}")
    done
else
    if [[ ! "${START_RUN}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: invalid start run number: ${START_RUN}"
        exit 1
    fi
    read -rp "Enter end run number [${START_RUN}]: " END_RUN
    END_RUN="${END_RUN:-${START_RUN}}"
    if [[ ! "${END_RUN}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: invalid end run number: ${END_RUN}"
        exit 1
    fi
    if (( 10#${END_RUN} < 10#${START_RUN} )); then
        echo "ERROR: end run number must be greater than or equal to start run number."
        exit 1
    fi

    WIDTH="${#START_RUN}"
    if (( ${#END_RUN} > WIDTH )); then
        WIDTH="${#END_RUN}"
    fi
    for ((run=10#${START_RUN}; run<=10#${END_RUN}; run++)); do
        RUNS+=("$(printf "%0${WIDTH}d" "${run}")")
    done
fi

while true; do
    read -rp "Enter replay mode (prad2, x17, or prad1) [prad2]: " REPLAY_MODE_INPUT
    REPLAY_MODE_INPUT="${REPLAY_MODE_INPUT,,}"
    case "${REPLAY_MODE_INPUT}" in
        ""|prad2|-prad2)
            REPLAY_MODE_INPUT_FOR_ONE="prad2"
            REPLAY_MODE_NAME="PRad2"
            break
            ;;
        x17|-x17)
            REPLAY_MODE_INPUT_FOR_ONE="x17"
            REPLAY_MODE_NAME="X17"
            break
            ;;
        prad1|-prad1)
            REPLAY_MODE_INPUT_FOR_ONE="prad1"
            REPLAY_MODE_NAME="PRad1"
            break
            ;;
        *)
            echo "ERROR: enter prad2, x17, or prad1."
            ;;
    esac
done
echo "Replay mode: ${REPLAY_MODE_NAME}"

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
read -rp "Cut JSON [default]: " CUT_INPUT
if [[ -z "${CUT_INPUT}" || "${CUT_INPUT}" == "default" ]]; then
    CUT_JSON="${DEFAULT_CUTS}"
    CUT_INPUT_FOR_ONE="default"
else
    CUT_JSON="${CUT_INPUT}"
    CUT_INPUT_FOR_ONE="${CUT_INPUT}"
fi

ROOT_SETUP="$(prompt_default "Enter ROOT setup script, or 'none' to rely on submitted environment" "${ROOT_SETUP:-none}")"
if [[ "${ROOT_SETUP}" == "none" ]]; then
    ROOT_SETUP=""
    ROOT_SETUP_FOR_ONE="none"
else
    ROOT_SETUP_FOR_ONE="${ROOT_SETUP}"
fi

SLURM_ACCOUNT="$(prompt_default "Enter Slurm account" "${SLURM_ACCOUNT}")"
SLURM_PARTITION="$(prompt_default "Enter Slurm partition (production or priority)" "${SLURM_PARTITION}")"
SLURM_TIME="$(prompt_default "Enter Slurm time limit" "${SLURM_TIME}")"
SLURM_MEM_PER_CPU="$(prompt_default "Enter Slurm mem-per-cpu MB" "${SLURM_MEM_PER_CPU}")"

PRAD2_SOFT="$(to_abs_path "${PRAD2_SOFT}" "${SUBMIT_DIR}")"
PRAD2_BIN="$(to_abs_path "${PRAD2_BIN}" "${SUBMIT_DIR}")"
CUT_JSON="$(to_abs_path "${CUT_JSON}" "${SUBMIT_DIR}")"
if [[ "${CUT_INPUT_FOR_ONE}" != "default" ]]; then
    CUT_INPUT_FOR_ONE="${CUT_JSON}"
fi
if [[ -n "${ROOT_SETUP}" ]]; then
    ROOT_SETUP="$(to_abs_path "${ROOT_SETUP}" "${SUBMIT_DIR}")"
    ROOT_SETUP_FOR_ONE="${ROOT_SETUP}"
fi
OUTPUT_BASE="$(to_abs_path "${OUTPUT_BASE}" "${SUBMIT_DIR}")"

echo ""
echo "Runs requested: $(join_by_space "${RUNS[@]}")"
echo "Checking cache under: ${CACHE_BASE}"

CACHED_RUNS=()
MISSING_RUNS=()
SUBMIT_FAILED_RUNS=()

for run in "${RUNS[@]}"; do
    RUN_DIR="${CACHE_BASE}/prad_${run}"
    if [[ -d "${RUN_DIR}" ]] && [[ "$(find "${RUN_DIR}" -maxdepth 1 -type f | wc -l)" -gt 0 ]]; then
        CACHED_RUNS+=("${run}")
    else
        MISSING_RUNS+=("${run}")
    fi
done

echo "Ready to submit: ${#CACHED_RUNS[@]} run(s)"
if [[ "${#CACHED_RUNS[@]}" -gt 0 ]]; then
    echo "  $(join_by_space "${CACHED_RUNS[@]}")"
fi
echo "Need jcache staging: ${#MISSING_RUNS[@]} run(s)"
if [[ "${#MISSING_RUNS[@]}" -gt 0 ]]; then
    echo "  $(join_by_space "${MISSING_RUNS[@]}")"
fi

echo ""
for run in "${CACHED_RUNS[@]}"; do
    echo "Submitting replay/recon job for run ${run}..."
    if ! printf '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n' \
            "${run}" \
            "${REPLAY_MODE_INPUT_FOR_ONE}" \
            "${PRAD2_SOFT}" \
            "${PRAD2_BIN}" \
            "${CACHE_BASE}" \
            "${OUTPUT_BASE}" \
            "${REPLAY_CORES}" \
            "${REPLAY_ZERO_SUPPRESS}" \
            "${REPLAY_MAX_FILES}" \
            "${REPLAY_MERGE_FILES}" \
            "${CUT_INPUT_FOR_ONE}" \
            "${ROOT_SETUP_FOR_ONE}" \
            "${SLURM_ACCOUNT}" \
            "${SLURM_PARTITION}" \
            "${SLURM_TIME}" \
            "${SLURM_MEM_PER_CPU}" | bash "${SUBMIT_ONE}"; then
        echo "WARNING: submit failed for run ${run}."
        SUBMIT_FAILED_RUNS+=("${run}")
    fi
    echo ""
done

if [[ "${#MISSING_RUNS[@]}" -gt 0 ]]; then
    echo "Some requested runs are not ready in cache."
    read -rp "Enter your email address for jcache notification: " USER_EMAIL
    if [[ -z "${USER_EMAIL}" ]]; then
        echo "ERROR: email cannot be empty."
        exit 1
    fi

    echo ""
    echo "Submitting jcache requests..."
    JCACHE_FAILED_RUNS=()
    for run in "${MISSING_RUNS[@]}"; do
        MSS_DIR="/mss/clas12/rg-o/data/prad_${run}"
        echo "jcache get ${MSS_DIR}/* -e ${USER_EMAIL}"
        if ! jcache get "${MSS_DIR}"/* -e "${USER_EMAIL}"; then
            echo "WARNING: jcache request failed for run ${run}."
            JCACHE_FAILED_RUNS+=("${run}")
        fi
    done

    echo ""
    echo "jcache requests finished."
    echo "Wait for staging, then re-run this script for these run(s):"
    echo "$(join_by_space "${MISSING_RUNS[@]}")"
    if [[ "${#JCACHE_FAILED_RUNS[@]}" -gt 0 ]]; then
        echo "jcache failed for these run(s); check manually:"
        echo "$(join_by_space "${JCACHE_FAILED_RUNS[@]}")"
    fi
else
    echo "All requested runs were submitted; no jcache staging is needed."
fi

if [[ "${#SUBMIT_FAILED_RUNS[@]}" -gt 0 ]]; then
    echo ""
    echo "Submit failed for these cached run(s); check the messages above:"
    echo "$(join_by_space "${SUBMIT_FAILED_RUNS[@]}")"
fi
