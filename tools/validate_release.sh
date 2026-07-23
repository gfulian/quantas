#!/usr/bin/env bash

# Validate Quantas from a freshly created development environment.
#
# Run this script from any location. It resolves the repository/source root from
# its own path, recreates a dedicated virtual environment, installs the current
# source tree with all developer dependencies, and stops at the first failure.

set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
VENV_DIR="${QUANTAS_VALIDATION_VENV:-${ROOT}/.venv-validation}"
LOG_DIR="${QUANTAS_VALIDATION_LOG_DIR:-${ROOT}/validation}"
STAGE_TIMEOUT="${QUANTAS_STAGE_TIMEOUT:-1200}"
TIMESTAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_FILE="${LOG_DIR}/quantas-validation-${TIMESTAMP}.log"

mkdir -p "${LOG_DIR}"

on_error() {
    local exit_code=$?
    printf '\nVALIDATION FAILED at line %s\nCommand: %s\nLog: %s\n' \
        "${BASH_LINENO[0]}" "${BASH_COMMAND}" "${LOG_FILE}" >&2
    exit "${exit_code}"
}
trap on_error ERR

exec > >(tee "${LOG_FILE}") 2>&1

cd "${ROOT}"

echo "=== Quantas fresh-source validation ==="
echo "Root: ${ROOT}"
echo "UTC start: $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
echo "Bootstrap Python: $(${PYTHON_BIN} --version 2>&1)"

if [[ ! -f "pyproject.toml" || ! -d "src/quantas" ]]; then
    echo "The resolved directory is not a Quantas source root: ${ROOT}" >&2
    exit 2
fi

rm -rf "${VENV_DIR}"
"${PYTHON_BIN}" -m venv "${VENV_DIR}"
# shellcheck disable=SC1091
source "${VENV_DIR}/bin/activate"

python -m pip install --upgrade pip setuptools wheel
python -m pip install -e ".[dev]"

VERSION="$(python -c 'import quantas; print(quantas.__version__)')"
echo "Quantas version: ${VERSION}"
echo "Validation Python: $(python --version 2>&1)"
python -m pip freeze

rm -rf build dist docs/_build
rm -rf .pytest_cache .mypy_cache .ruff_cache

for directory in src tests tools examples docs/tools; do
    if [[ -d "${directory}" ]]; then
        find "${directory}" -type d -name "__pycache__" -prune -exec rm -rf {} +
        find "${directory}" -type f \( -name "*.pyc" -o -name "*.pyo" \) -delete
    fi
done

python tools/update_examples_manifest.py --check

ruff check src tests tools docs/tools
mypy
python -m compileall -q src tests tools docs/tools

python tools/run_tests.py all --stage-timeout "${STAGE_TIMEOUT}" -- -q
python tools/update_examples_manifest.py --check

python tools/check_architecture.py --root .

python -m sphinx \
    -E -a \
    -b html \
    -W \
    --keep-going \
    docs/source \
    docs/_build/html

python -m build
python -m twine check dist/*
python tools/check_distribution.py dist
python tools/check_repository.py

echo "UTC finish: $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
echo "=== Quantas ${VERSION} validation completed successfully ==="
echo "Log: ${LOG_FILE}"
