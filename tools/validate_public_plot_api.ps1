<#
.SYNOPSIS
Validate the in-development Quantas 2.0.0b7 public plotting API on Windows.

.DESCRIPTION
Runs the focused API, registry, CLI, and plotting tests required by the public
plot-contract work. With -Full, it also runs Ruff, mypy, the complete staged
pytest suite, the architecture audit, and the warning-as-error Sphinx build.
The script must be executed from an environment installed with .[dev].
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

Invoke-Checked "Python environment" {
    python -c "import sys, quantas; print(sys.version); print('Quantas', quantas.__version__)"
}

Invoke-Checked "Source compilation" {
    python -m compileall -q src tests tools docs/tools
}

Invoke-Checked "Public plotting API and registry" {
    python -m pytest -q `
        tests/infrastructure/test_api_surface.py `
        tests/infrastructure/test_public_api.py `
        tests/infrastructure/test_public_plotting.py `
        tests/infrastructure/test_public_rendering.py `
        tests/infrastructure/test_plot_inventory.py `
        tests/infrastructure/test_module_contracts.py `
        tests/infrastructure/test_cli_contracts.py `
        tests/infrastructure/test_api_reference_documentation.py `
        tests/infrastructure/test_source_hygiene.py `
        tests/infrastructure/test_version.py
}

Invoke-Checked "Module plotting contracts" {
    python -m pytest -q `
        tests/modules/elasticity/test_plotting.py `
        tests/modules/seismic/test_outputs.py `
        tests/modules/ha/test_plotting.py `
        tests/modules/qha/test_plotting.py `
        tests/modules/thermoelasticity/test_plotting.py `
        tests/modules/thermoelasticity/test_cli.py `
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

Write-Host "`nAll requested public plotting API checks passed." -ForegroundColor Green
