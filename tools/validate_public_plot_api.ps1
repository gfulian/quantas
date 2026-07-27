<#
.SYNOPSIS
Compatibility entry point for the renamed public lifecycle API validation.
#>

[CmdletBinding()]
param(
    [switch]$Full,
    [int]$StageTimeout = 1200
)

& "$PSScriptRoot\validate_public_lifecycle_api.ps1" `
    -Full:$Full `
    -StageTimeout $StageTimeout
exit $LASTEXITCODE
