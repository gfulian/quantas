<#
.SYNOPSIS
Validate the in-development Quantas 2.0.0b10 public lifecycle API on Windows.

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
    python -c @"
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
"@
}

Invoke-Checked "Source compilation" {
    python -m compileall -q src tests tools docs/tools
}

Invoke-Checked "Public lifecycle surface and documentation" {
    python -m pytest -q `
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
    python -m pytest -q `
        tests/modules/elasticity/test_cli.py `
        tests/modules/elasticity/test_run_rotations.py `
        tests/modules/ha/test_input_generation.py `
        tests/modules/ha/test_cli.py `
        tests/modules/qha/test_cli.py `
        tests/modules/thermoelasticity/test_input_api.py
}

Invoke-Checked "Public exports and persistence" {
    python -m pytest -q `
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
    python -m pytest -q `
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
    python tools/check_architecture.py --root .
}

if ($Full) {
    Invoke-Checked "Ruff" {
        ruff check src tests tools docs/tools
    }
    Invoke-Checked "mypy" {
        mypy
    }
    Invoke-Checked "Complete staged test suite" {
        python tools/run_tests.py all --stage-timeout $StageTimeout -- -q
    }
    Invoke-Checked "Sphinx warning-as-error build" {
        python -m sphinx -E -a -W --keep-going -b html docs/source docs/_build/html
    }
    Invoke-Checked "Distribution build" {
        python -m build
        $Artifacts = Get-ChildItem -Path dist -File
        if (-not $Artifacts) {
            throw "Distribution build produced no files under dist"
        }
        python -m twine check @($Artifacts.FullName)
        python tools/check_distribution.py dist
    }
}

Write-Host "`nAll requested public lifecycle API checks passed." -ForegroundColor Green
