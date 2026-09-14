# -*- coding: utf-8 -*-

"""Shell-completion support not provided directly by the pinned Click release.

Click 8.1 ships built-in completion backends for Bash, Zsh, and Fish, but not
for PowerShell.  Quantas registers a small PowerShell backend that reuses
Click's own context resolution and every parameter ``shell_complete`` method.
The scientific CLI therefore has one completion source rather than a separate
PowerShell-specific catalogue.
"""

from __future__ import annotations

import json
import os

from click.shell_completion import CompletionItem, ShellComplete, add_completion_class


@add_completion_class
class PowerShellComplete(ShellComplete):
    """Provide native PowerShell completion through ``Register-ArgumentCompleter``.

    Notes
    -----
    The generated script delegates candidate discovery back to Click.  File and
    directory completion markers are handed to PowerShell's own filename
    completer, while ordinary candidates retain the descriptions supplied by
    Quantas parameter types.
    """

    name = "powershell"

    source_template = r'''Register-ArgumentCompleter -Native -CommandName %(prog_name)s -ScriptBlock {
    param($wordToComplete, $commandAst, $cursorPosition)

    $elements = @()
    foreach ($element in $commandAst.CommandElements) {
        if ($element -is [System.Management.Automation.Language.StringConstantExpressionAst]) {
            $elements += $element.Value
        }
        else {
            $elements += $element.Extent.Text
        }
    }

    if ($wordToComplete -ne "" -and $elements.Count -gt 1) {
        $completeArgs = @($elements[0..($elements.Count - 2)])
    }
    else {
        $completeArgs = @($elements)
    }

    try {
        $env:QUANTAS_COMPLETE_ARGS = ConvertTo-Json -Compress -InputObject @($completeArgs)
        $env:QUANTAS_COMPLETE_WORD = $wordToComplete
        $env:%(complete_var)s = "powershell_complete"
        $response = & %(prog_name)s 2>$null
    }
    finally {
        Remove-Item Env:QUANTAS_COMPLETE_ARGS -ErrorAction SilentlyContinue
        Remove-Item Env:QUANTAS_COMPLETE_WORD -ErrorAction SilentlyContinue
        Remove-Item Env:%(complete_var)s -ErrorAction SilentlyContinue
    }

    foreach ($line in $response) {
        if ([string]::IsNullOrWhiteSpace($line)) {
            continue
        }
        $item = $line | ConvertFrom-Json
        if ($item.type -eq "file" -or $item.type -eq "dir") {
            [System.Management.Automation.CompletionCompleters]::CompleteFilename($wordToComplete)
            continue
        }
        $tooltip = if ([string]::IsNullOrWhiteSpace($item.help)) { $item.value } else { $item.help }
        [System.Management.Automation.CompletionResult]::new(
            $item.value,
            $item.value,
            [System.Management.Automation.CompletionResultType]::ParameterValue,
            $tooltip
        )
    }
}
'''

    def get_completion_args(self) -> tuple[list[str], str]:
        """Return complete command arguments and the current PowerShell word.

        Returns
        -------
        tuple of list of str and str
            Arguments after the executable name and the incomplete word.
        """
        raw = os.environ.get("QUANTAS_COMPLETE_ARGS", "[]")
        values = json.loads(raw)
        if not isinstance(values, list):
            values = []
        words = [str(value) for value in values]
        if words and words[0].lower().endswith(("quantas", "quantas.exe")):
            words = words[1:]
        return words, os.environ.get("QUANTAS_COMPLETE_WORD", "")

    def format_completion(self, item: CompletionItem) -> str:
        """Serialize one completion item as one JSON object."""
        return json.dumps(
            {"type": item.type, "value": item.value, "help": item.help},
            ensure_ascii=False,
            separators=(",", ":"),
        )


__all__ = ["PowerShellComplete"]
