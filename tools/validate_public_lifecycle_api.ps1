<#
.SYNOPSIS
Validate the in-development Quantas 2.0.0b9 public lifecycle API on Windows.

.DESCRIPTION
Runs focused public input, execution, persistence, report, plot, export,
registry, and CLI/API-equivalence checks. With -Full, it also runs Ruff, mypy,
the complete staged pytest suite, the architecture audit, warning-as-error
Sphinx documentation, and distribution validation. Execute the script from an
environment installed with .[dev].
#>

[CmdletBinding()]
param(
    [switch]$Full,
    [int]$StageTimeout = 1200
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $PSScriptRoot
Set-Location $Root

$Python = Join-Path $Root ".venv\Scripts\python.exe"
if (-not (Test-Path $Python)) {
    throw @"
The expected virtual-environment interpreter was not found:

  $Python

Create the worktree environment first, for example:

  py -3.10 -m venv .venv
  & ".\.venv\Scripts\python.exe" -m pip install -e ".[dev]"
"@
}

function Invoke-Checked {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Description,
        [Parameter(Mandatory = $true)]
        [scriptblock]$Command
    )
    Write-Host "`n=== $Description ===" -ForegroundColor Cyan
    & $Command
    if ($LASTEXITCODE -ne 0) {
        throw "$Description failed with exit code $LASTEXITCODE"
    }
}

Invoke-Checked "Python environment and required backends" {
    $PythonCheck = @'
import sys
import quantas
import numpy
import spglib
import odrpack

print(sys.version)
print("Python executable:", sys.executable)
print("Quantas:", quantas.__version__)
print("NumPy:", numpy.__version__)
print("spglib:", getattr(spglib, "__version__", "available"))
print("odrpack:", getattr(odrpack, "__version__", "available"))
'@
    $PythonCheck | & $Python -
}

Invoke-Checked "Source compilation" {
    & $Python -m compileall -q src tests tools docs/tools
}

Invoke-Checked "Public lifecycle surface and documentation" {
    & $Python -m pytest -q `
        tests/infrastructure/test_api_surface.py `
        tests/infrastructure/test_public_api.py `
        tests/infrastructure/test_public_lifecycle.py `
        tests/infrastructure/test_api_reference_documentation.py `
        tests/infrastructure/test_public_plotting.py `
        tests/infrastructure/test_public_rendering.py `
        tests/infrastructure/test_plot_inventory.py `
        tests/infrastructure/test_module_contracts.py `
        tests/infrastructure/test_cli_contracts.py `
        tests/infrastructure/test_source_hygiene.py `
        tests/infrastructure/test_version.py
}

Invoke-Checked "Input generation and CLI/API equivalence" {
    & $Python -m pytest -q `
        tests/modules/elasticity/test_cli.py `
        tests/modules/elasticity/test_run_rotations.py `
        tests/modules/ha/test_input_generation.py `
        tests/modules/ha/test_cli.py `
        tests/modules/qha/test_cli.py `
        tests/modules/thermoelasticity/test_input_api.py
}

Invoke-Checked "Public exports and persistence" {
    & $Python -m pytest -q `
        tests/modules/elasticity/test_export.py `
        tests/modules/ha/test_export.py `
        tests/modules/qha/test_export.py `
        tests/modules/qha/test_hdf5_export.py `
        tests/modules/qha/test_table_export.py `
        tests/modules/seismic/test_outputs.py `
        tests/modules/thermoelasticity/test_cli.py `
        tests/modules/thermoelasticity/test_qsa.py `
        tests/modules/eos/test_diagnostics.py `
        tests/modules/eos/test_cli_run.py
}

Invoke-Checked "Frontend-neutral plot contracts" {
    & $Python -m pytest -q `
        tests/modules/elasticity/test_plotting.py `
        tests/modules/seismic/test_outputs.py `
        tests/modules/ha/test_plotting.py `
        tests/modules/qha/test_plotting.py `
        tests/modules/thermoelasticity/test_plotting.py `
        tests/modules/thermoelasticity/test_cli.py `
        tests/modules/eos/test_plot_inventory.py `
        tests/modules/eos/test_plotting.py
}

Invoke-Checked "Architecture audit" {
    & $Python tools/check_architecture.py --root .
}

if ($Full) {
    Invoke-Checked "Ruff" {
        & $Python -m ruff check src tests tools docs/tools
    }
    Invoke-Checked "mypy" {
        & $Python -m mypy
    }
    Invoke-Checked "Complete staged test suite" {
        & $Python tools/run_tests.py all --stage-timeout $StageTimeout -- -q
    }
    Invoke-Checked "Sphinx warning-as-error build" {
        & $Python -m sphinx -E -a -W --keep-going -b html docs/source docs/_build/html
    }
    Invoke-Checked "Distribution build" {
        & $Python -m build
        $Artifacts = Get-ChildItem -Path dist -File
        if (-not $Artifacts) {
            throw "Distribution build produced no files under dist"
        }
        & $Python -m twine check @($Artifacts.FullName)
        & $Python tools/check_distribution.py dist
    }
}

Write-Host "`nAll requested public lifecycle API checks passed." -ForegroundColor Green
