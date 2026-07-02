#!/bin/bash
# replay_recon.sh — run the single-run PRad-II replay pipeline on JLab ifarm
#
# Usage: chmod +x scripts/shell/replay_recon.sh 
#       ./replay_recon.sh
#   Prompts for run number, then:
#     1. Counts EVIO files in /cache/clas12/rg-o/data/prad_<RUN>/
#     2. Runs prad2ana_replay_recon; by default, recon ROOTs are merged in
#        groups of 62 into prad_<RUN>_recon_<NNN>.root
#     3. Runs prad2ana_replay_filter on all recon ROOTs and writes the
#        corresponding *_filter.root files directly under the run output dir
#     4. Runs prad2ana_live_charge over all filtered ROOTs together
#     5. Runs prad2ana_quick_check over all recon ROOTs
#
# The one "parallel jobs" prompt is reused for replay_recon -j,
# replay_filter -t, and quick_check -j. All products are written under
# <output_base>/prad_<RUN>/; no filter subdirectory is created.

# ---------------------------------------------------------------------------
# Configurable defaults (override via environment variables)
# ---------------------------------------------------------------------------
PRAD2_SOFT=${PRAD2_SOFT:-./prad2evviewer}
PRAD2_BIN="${PRAD2_BIN:-${PRAD2_SOFT}/build/bin}"
CACHE_BASE="${CACHE_BASE:-/cache/clas12/rg-o/data}"

# Default output base directory (can be overridden by user input)
OUTPUT_BASE="${OUTPUT_BASE:-./}"
REPLAY_CORES="${REPLAY_CORES:-15}"
REPLAY_ZERO_SUPPRESS="${REPLAY_ZERO_SUPPRESS:-5}"
REPLAY_MAX_FILES="${REPLAY_MAX_FILES:-10000}"
REPLAY_MERGE_FILES="${REPLAY_MERGE_FILES:-62}"
DEFAULT_CUTS="${DEFAULT_CUTS:-${PRAD2_SOFT}/analysis/cuts/prad2_default.json}"

# ---------------------------------------------------------------------------
# 0. Check ROOT environment
# ---------------------------------------------------------------------------
if ! command -v root &>/dev/null; then
    echo "ERROR: 'root' is not found in PATH. Please set up the ROOT environment first."
    echo "       e.g.: source /path/to/root/bin/thisroot.sh"
    exit 1
fi
if ! root -l -q &>/dev/null; then
    echo "ERROR: 'root -l' failed. The ROOT installation may be incomplete or misconfigured."
    exit 1
fi
echo "ROOT is available: $(root-config --version 2>/dev/null || root --version 2>&1 | head -1)"

# ---------------------------------------------------------------------------
# 1. Ask for parameters
# ---------------------------------------------------------------------------
read -rp "Enter run number (e.g. 024650): " RUN_NUMBER
if [[ -z "${RUN_NUMBER}" ]]; then
    echo "ERROR: run number cannot be empty."
    exit 1
fi

read -rp "Enter output base directory [${OUTPUT_BASE}]: " _INPUT
[[ -n "${_INPUT}" ]] && OUTPUT_BASE="${_INPUT}"

read -rp "Enter number of parallel jobs (-j) [${REPLAY_CORES}]: " _INPUT
[[ -n "${_INPUT}" ]] && REPLAY_CORES="${_INPUT}"

read -rp "Enter GEM zero suppression (-z) [${REPLAY_ZERO_SUPPRESS}]: " _INPUT
[[ -n "${_INPUT}" ]] && REPLAY_ZERO_SUPPRESS="${_INPUT}"

read -rp "Enter max number of files to process (-f) [${REPLAY_MAX_FILES}]: " _INPUT
[[ -n "${_INPUT}" ]] && REPLAY_MAX_FILES="${_INPUT}"

read -rp "Enter replay merge group size (-m, 0 disables) [${REPLAY_MERGE_FILES}]: " _INPUT
[[ -n "${_INPUT}" ]] && REPLAY_MERGE_FILES="${_INPUT}"

echo "Enter cut JSON file for replay_filter (path to cuts.json, or 'default' to use ${DEFAULT_CUTS}):"
read -rp "Cut JSON [default]: " _INPUT
if [[ -z "${_INPUT}" || "${_INPUT}" == "default" ]]; then
    CUT_JSON="${DEFAULT_CUTS}"
else
    CUT_JSON="${_INPUT}"
fi

RUN_DIR="${CACHE_BASE}/prad_${RUN_NUMBER}"

# ---------------------------------------------------------------------------
# 2. Check input files
# ---------------------------------------------------------------------------
echo ""
echo "Checking input directory: ${RUN_DIR}"

if [[ ! -d "${RUN_DIR}" ]]; then
    echo "ERROR: directory not found: ${RUN_DIR}"
    exit 1
fi

FILE_COUNT=$(ls "${RUN_DIR}" 2>/dev/null | wc -l)
echo "Found ${FILE_COUNT} file(s) in ${RUN_DIR}:"
ls "${RUN_DIR}"
echo ""

if [[ "${FILE_COUNT}" -eq 0 ]]; then
    echo "No files found in ${RUN_DIR}."
    echo "The data may still be on tape (MSS). Submitting a jcache request to stage them..."
    echo ""
    read -rp "Enter your email address to receive a notification when staging is complete: " USER_EMAIL
    if [[ -z "${USER_EMAIL}" ]]; then
        echo "ERROR: email cannot be empty."
        exit 1
    fi
    MSS_DIR="/mss/clas12/rg-o/data/prad_${RUN_NUMBER}"
    echo ""
    echo "Command: jcache get ${MSS_DIR}/* -e ${USER_EMAIL}"
    jcache get ${MSS_DIR}/* -e "${USER_EMAIL}"
    echo ""
    echo "jcache request submitted."
    echo "You will receive an email at ${USER_EMAIL} when the files are staged to cache."
    echo "Please re-run this script after receiving the notification."
    exit 0
fi

# ---------------------------------------------------------------------------
# 3. Prepare output directory
# ---------------------------------------------------------------------------
OUT_DIR="${OUTPUT_BASE}/prad_${RUN_NUMBER}"
mkdir -p "${OUT_DIR}"
echo "Output directory: ${OUT_DIR}"
echo ""

# ---------------------------------------------------------------------------
# 4. Run replay
# ---------------------------------------------------------------------------
REPLAY_CMD="${PRAD2_BIN}/prad2ana_replay_recon"

if [[ ! -x "${REPLAY_CMD}" ]]; then
    echo "ERROR: executable not found: ${REPLAY_CMD}"
    exit 1
fi

echo "Starting replay..."
echo "Command: ${REPLAY_CMD} ${RUN_DIR} -o ${OUT_DIR} -j ${REPLAY_CORES} -z ${REPLAY_ZERO_SUPPRESS} -f ${REPLAY_MAX_FILES} -m ${REPLAY_MERGE_FILES}"
echo ""

"${REPLAY_CMD}" "${RUN_DIR}" -o "${OUT_DIR}" -j "${REPLAY_CORES}" -z "${REPLAY_ZERO_SUPPRESS}" -f "${REPLAY_MAX_FILES}" -m "${REPLAY_MERGE_FILES}"
REPLAY_EXIT=$?

if [[ "${REPLAY_EXIT}" -ne 0 ]]; then
    echo ""
    echo "ERROR: replay exited with code ${REPLAY_EXIT}."
    exit "${REPLAY_EXIT}"
fi

echo ""
echo "Replay finished successfully."

# ---------------------------------------------------------------------------
# 5. Select reconstructed ROOT files for downstream tools
# ---------------------------------------------------------------------------
mapfile -t RECON_INPUTS < <(find "${OUT_DIR}" -maxdepth 1 -type f -name "prad_${RUN_NUMBER}_recon_*.root" | sort)
if [[ "${#RECON_INPUTS[@]}" -eq 0 ]]; then
    mapfile -t RECON_INPUTS < <(find "${OUT_DIR}" -maxdepth 1 -type f -name "prad_${RUN_NUMBER}.*_recon.root" | sort)
fi
if [[ "${#RECON_INPUTS[@]}" -eq 0 ]]; then
    echo "ERROR: no reconstructed ROOT files found in ${OUT_DIR}."
    exit 1
fi
echo "Downstream input ROOT file(s): ${#RECON_INPUTS[@]}"

# ---------------------------------------------------------------------------
# 6. Run replay filter
# ---------------------------------------------------------------------------
FILTER_CMD="${PRAD2_BIN}/prad2ana_replay_filter"

if [[ ! -x "${FILTER_CMD}" ]]; then
    echo "WARNING: executable not found: ${FILTER_CMD}, skipping replay filter."
else
    if [[ ! -f "${CUT_JSON}" ]]; then
        echo "WARNING: cut JSON not found: ${CUT_JSON}, skipping replay filter."
    else
        REPORT="${OUT_DIR}/prad_${RUN_NUMBER}_filter_report.json"
        SLOW_ROOT="${OUT_DIR}/prad_${RUN_NUMBER}_epics.root"
        echo ""
        echo "Running replay filter..."
        if [[ "${#RECON_INPUTS[@]}" -gt 1 ]]; then
            FILTER_OUT="${OUT_DIR}"
        else
            _BASE="$(basename "${RECON_INPUTS[0]}")"
            if [[ "${_BASE}" =~ ^(.+)_recon_([^/]+)\.root$ ]]; then
                FILTER_OUT="${OUT_DIR}/${BASH_REMATCH[1]}_filter_${BASH_REMATCH[2]}.root"
            elif [[ "${_BASE}" =~ ^(.+)\.evio\.([0-9]+)_recon\.root$ ]]; then
                FILTER_OUT="${OUT_DIR}/${BASH_REMATCH[1]}_filter_${BASH_REMATCH[2]}.root"
            elif [[ "${_BASE}" =~ ^(.+)\.([0-9]+)_recon\.root$ ]]; then
                FILTER_OUT="${OUT_DIR}/${BASH_REMATCH[1]}_filter_${BASH_REMATCH[2]}.root"
            else
                FILTER_OUT="${OUT_DIR}/${_BASE%.root}_filter.root"
            fi
        fi
        echo "Command: ${FILTER_CMD} ${RECON_INPUTS[*]} -o ${FILTER_OUT} -c ${CUT_JSON} -j ${REPORT} -t ${REPLAY_CORES}"
        echo ""

        "${FILTER_CMD}" "${RECON_INPUTS[@]}" -o "${FILTER_OUT}" -c "${CUT_JSON}" -j "${REPORT}" -t "${REPLAY_CORES}"
        FILTER_EXIT=$?

        if [[ "${FILTER_EXIT}" -ne 0 ]]; then
            echo ""
            echo "ERROR: replay filter exited with code ${FILTER_EXIT}."
            exit "${FILTER_EXIT}"
        fi

        echo ""
        echo "Replay filter finished."
        echo "Filtered out  : ${FILTER_OUT}"
        echo "Slow ROOT     : ${SLOW_ROOT}"
        echo "Filter report : ${REPORT}"
    fi
fi

# ---------------------------------------------------------------------------
# 7. Run live charge calculation
# ---------------------------------------------------------------------------
LIVE_CHARGE_CMD="${PRAD2_BIN}/prad2ana_live_charge"

if [[ ! -x "${LIVE_CHARGE_CMD}" ]]; then
    echo "WARNING: executable not found: ${LIVE_CHARGE_CMD}, skipping live charge."
else
    mapfile -t LC_INPUTS < <(find "${OUT_DIR}" -maxdepth 1 -type f -name "prad_${RUN_NUMBER}_filter*.root" | sort)
    if [[ "${#LC_INPUTS[@]}" -eq 0 ]]; then
        echo "ERROR: no filtered ROOT files found in ${OUT_DIR}; live_charge requires prad_${RUN_NUMBER}_filter*.root inputs."
        exit 1
    fi
    LC_JSON="${OUT_DIR}/prad_${RUN_NUMBER}_live_charge.json"
    echo ""
    echo "Running live charge calculation..."
    echo "Command: ${LIVE_CHARGE_CMD} ${LC_INPUTS[*]} -j ${LC_JSON}"
    echo ""

    "${LIVE_CHARGE_CMD}" "${LC_INPUTS[@]}" -j "${LC_JSON}"
    LC_EXIT=$?

    if [[ "${LC_EXIT}" -ne 0 ]]; then
        echo ""
        echo "ERROR: live charge exited with code ${LC_EXIT}."
        exit "${LC_EXIT}"
    fi

    echo ""
    echo "Live charge finished."
    echo "Live charge JSON: ${LC_JSON}"
fi

# ---------------------------------------------------------------------------
# 8. Run quick check
# ---------------------------------------------------------------------------
QUICK_CHECK_CMD="${PRAD2_BIN}/prad2ana_quick_check"

if [[ ! -x "${QUICK_CHECK_CMD}" ]]; then
    echo "WARNING: executable not found: ${QUICK_CHECK_CMD}, skipping quick check."
else
    echo ""
    echo "Running quick check..."
    QC_OUTPUT="${OUT_DIR}/prad_${RUN_NUMBER}_quick_check.root"
    echo "Command: ${QUICK_CHECK_CMD} ${RECON_INPUTS[*]} -o ${QC_OUTPUT} -j ${REPLAY_CORES}"
    echo ""

    "${QUICK_CHECK_CMD}" "${RECON_INPUTS[@]}" -o "${QC_OUTPUT}" -j "${REPLAY_CORES}"
    QC_EXIT=$?

    if [[ "${QC_EXIT}" -ne 0 ]]; then
        echo ""
        echo "ERROR: quick check exited with code ${QC_EXIT}."
        exit "${QC_EXIT}"
    fi

    echo ""
    echo "Quick check finished."
fi
