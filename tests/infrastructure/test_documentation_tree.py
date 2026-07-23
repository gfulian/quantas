"""Structural tests for the public Sphinx documentation tree."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import posixpath
import re


DOCS_ROOT = Path("docs/source")


def _rst_documents() -> dict[str, Path]:
    """Return source document names mapped to RST paths."""
    return {
        path.relative_to(DOCS_ROOT).with_suffix("").as_posix(): path
        for path in DOCS_ROOT.rglob("*.rst")
    }


def _toctree_targets(path: Path) -> list[str]:
    """Return normalized document targets from all toctrees in one page."""
    targets: list[str] = []
    lines = path.read_text(encoding="utf-8").splitlines()
    index = 0
    while index < len(lines):
        if not lines[index].lstrip().startswith(".. toctree::"):
            index += 1
            continue
        index += 1
        while index < len(lines) and (
            not lines[index].strip() or lines[index].lstrip().startswith(":")
        ):
            index += 1
        while index < len(lines):
            line = lines[index]
            if not line.startswith("   ") or not line.strip():
                break
            target = line.strip()
            if "<" in target and target.endswith(">"):
                target = target.rsplit("<", 1)[1][:-1]
            base = path.relative_to(DOCS_ROOT).parent.as_posix()
            targets.append(posixpath.normpath(posixpath.join(base, target)))
            index += 1
    return targets


def test_root_toctree_exposes_the_approved_manual_sections() -> None:
    """The sidebar presents the approved scientific-manual hierarchy."""
    text = (DOCS_ROOT / "index.rst").read_text(encoding="utf-8")
    captions = re.findall(r"^   :caption: (.+)$", text, flags=re.MULTILINE)
    assert captions == [
        "INTRODUCTION",
        "GETTING STARTED",
        "SCIENTIFIC BACKGROUND",
        "IMPLEMENTATION AND WORKFLOWS",
        "TUTORIALS",
        "INPUT AND OUTPUT FORMATS",
        "COMMAND REFERENCE",
        "API REFERENCE",
        "SCIENTIFIC VALIDATION",
        "DEVELOPMENT GUIDE",
    ]


def test_docs_toctree_targets_are_unique() -> None:
    """Every navigation target exists and belongs to only one toctree."""
    documents = _rst_documents()
    targets = [
        target
        for path in documents.values()
        for target in _toctree_targets(path)
    ]
    assert set(targets) <= set(documents)
    duplicates = {target for target, count in Counter(targets).items() if count > 1}
    assert duplicates == set()


def test_every_rst_page_is_navigated_or_explicitly_orphaned() -> None:
    """No public page silently disappears from the manual navigation."""
    documents = _rst_documents()
    included = {
        target
        for path in documents.values()
        for target in _toctree_targets(path)
    }
    for name, path in documents.items():
        if name == "index" or name in included:
            continue
        assert path.read_text(encoding="utf-8").startswith(":orphan:")


def test_all_doc_roles_resolve_to_existing_pages() -> None:
    """Cross-references remain valid after documentation reorganization."""
    documents = _rst_documents()
    for name, path in documents.items():
        text = path.read_text(encoding="utf-8")
        for raw_target in re.findall(r":doc:`([^`]+)`", text):
            target = raw_target
            if "<" in target and target.endswith(">"):
                target = target.rsplit("<", 1)[1][:-1]
            target = target.strip()
            if target.startswith("/"):
                resolved = target[1:]
            else:
                resolved = posixpath.normpath(
                    posixpath.join(posixpath.dirname(name), target)
                )
            assert resolved in documents, f"{path}: missing :doc: target {raw_target!r}"


def test_ha_qha_tutorial_assets_and_downloads_exist() -> None:
    """Completed HA/QHA tutorials keep all referenced local assets available."""
    expected = [
        DOCS_ROOT / "_downloads/mgo_b3lyp.yaml",
        DOCS_ROOT / "_downloads/tutorials/ha/tutorial_api.py",
        DOCS_ROOT / "_downloads/tutorials/qha/tutorial_api.py",
        DOCS_ROOT / "_static/tutorials/ha/mgo_ha_isochoric_heat_capacity.png",
        DOCS_ROOT / "_static/tutorials/ha/mgo_ha_vibrational_free_energy.png",
        DOCS_ROOT / "_static/tutorials/qha/mgo_qha_equilibrium_volume.png",
        DOCS_ROOT / "_static/tutorials/qha/mgo_qha_heat_capacities.png",
        DOCS_ROOT / "_static/tutorials/qha/mgo_qha_thermal_expansion_map.png",
    ]
    assert all(path.is_file() and path.stat().st_size > 0 for path in expected)


def test_ha_qha_tutorials_are_complete_public_workflows() -> None:
    """HA/QHA tutorials must remain detailed CLI/API examples, not placeholders."""
    for name in ("ha", "qha"):
        text = (DOCS_ROOT / "tutorials" / f"{name}.rst").read_text(encoding="utf-8")
        assert "Work in progress" not in text
        assert "Running the same calculation from Python" in text
        assert f"quantas {name} " in text
        assert "Reproducibility checkpoints" in text


def test_qha_tutorial_contains_method_sensitivity_exercise() -> None:
    """The QHA tutorial compares interpolation and minimization choices."""
    text = (DOCS_ROOT / "tutorials/qha.rst").read_text(encoding="utf-8")
    assert "Exercise: sensitivity to interpolation and minimization" in text
    assert "--scheme freq --minimization poly" in text
    assert "--scheme td --minimization eos --eos BM3" in text
    assert (DOCS_ROOT / "_downloads/tutorials/qha/compare_methods.py").is_file()



def test_elasticity_seismic_tutorial_assets_and_downloads_exist() -> None:
    """Completed Elasticity/SEISMIC tutorials keep local assets available."""
    expected = [
        DOCS_ROOT / "_downloads/tutorials/elasticity/tutorial_api.py",
        DOCS_ROOT / "_downloads/tutorials/seismic/tutorial_api.py",
        DOCS_ROOT / "_static/tutorials/elasticity/calcite_young_modulus_2d.png",
        DOCS_ROOT / "_static/tutorials/elasticity/calcite_young_modulus_3d.png",
        DOCS_ROOT / "_static/tutorials/seismic/hydroxylapatite_summary.png",
        DOCS_ROOT / "_static/tutorials/seismic/hydroxylapatite_vs1_polarization.png",
        DOCS_ROOT / "_static/tutorials/seismic/hydroxylapatite_vp_phase_surface.png",
    ]
    assert all(path.is_file() and path.stat().st_size > 0 for path in expected)


def test_elasticity_seismic_tutorials_are_complete() -> None:
    """Elasticity and SEISMIC tutorials remain detailed CLI/API examples."""
    expectations = {
        "elasticity": (
            "Running Elasticity from the command line",
            "Running the same workflow from Python",
            "quantas elasticity run",
            "Reproducibility checkpoints",
        ),
        "seismic": (
            "Running SEISMIC from the command line",
            "Running the same workflow from Python",
            "quantas seismic run",
            "Reproducibility checkpoints",
        ),
    }
    for name, phrases in expectations.items():
        content = (DOCS_ROOT / "tutorials" / f"{name}.rst").read_text(
            encoding="utf-8"
        )
        assert "Work in progress" not in content
        for phrase in phrases:
            assert phrase in content


def test_elasticity_seismic_tutorial_scripts_use_public_api_only() -> None:
    """Distributed tutorial scripts must not import implementation packages."""
    for path in (
        Path("examples/elasticity/tutorial_api.py"),
        Path("examples/seismic/tutorial_api.py"),
    ):
        content = path.read_text(encoding="utf-8")
        assert "quantas.core" not in content
        assert "quantas.modules" not in content
        assert "quantas.renderers" not in content
        assert "from quantas.api" in content

def test_elasticity_and_seismic_background_are_complete() -> None:
    """Elasticity and SEISMIC theory pages cover the implemented physics."""
    elasticity = (DOCS_ROOT / "theory/elasticity.rst").read_text(encoding="utf-8")
    seismic = (DOCS_ROOT / "theory/seismic.rst").read_text(encoding="utf-8")
    for phrase in (
        "Voigt representation and shear convention",
        "Mechanical stability",
        "Directional elastic properties",
        "Transverse extrema",
        "Tensor rotations and reference frames",
    ):
        assert phrase in elasticity
    for phrase in (
        "Christoffel matrix",
        "Group velocity and ray direction",
        "Degeneracy and near-degeneracy",
        "Acoustic enhancement and phonon focusing",
        "Caustic candidates",
    ):
        assert phrase in seismic


def test_thermoelasticity_tutorial_assets_and_downloads_exist() -> None:
    """The guided thermoelasticity tutorial retains every local artifact."""
    expected = [
        DOCS_ROOT / "_downloads/tutorials/thermoelasticity/tutorial_api.py",
        DOCS_ROOT / "_downloads/tutorials/thermoelasticity/dol_pbe0_qha.yaml",
        DOCS_ROOT / "_downloads/tutorials/thermoelasticity/dol_pbe0_thermoelastic.yaml",
        DOCS_ROOT / "_downloads/tutorials/thermoelasticity/dolomite_continental_profile.yaml",
        DOCS_ROOT / "_downloads/tutorials/thermoelasticity/dolomite_grid.csv",
        DOCS_ROOT / "_downloads/tutorials/thermoelasticity/dolomite_profile.csv",
        DOCS_ROOT / "_downloads/tutorials/thermoelasticity/dolomite_5GPa_800K.dat",
        DOCS_ROOT / "_static/tutorials/thermoelasticity/dolomite_c11_fit.png",
        DOCS_ROOT / "_static/tutorials/thermoelasticity/dolomite_pt_c11_c33.png",
        DOCS_ROOT / "_static/tutorials/thermoelasticity/dolomite_isothermal_adiabatic.png",
        DOCS_ROOT / "_static/tutorials/thermoelasticity/dolomite_profile_relative.png",
        DOCS_ROOT / "_static/tutorials/thermoelasticity/dolomite_pt_domain.png",
    ]
    assert all(path.is_file() and path.stat().st_size > 0 for path in expected)


def test_thermoelasticity_tutorial_is_a_complete_staged_workflow() -> None:
    """The tutorial must guide users through calibration and all analyses."""
    text = (DOCS_ROOT / "tutorials/thermoelasticity.rst").read_text(
        encoding="utf-8"
    )
    assert "Work in progress" not in text
    for phrase in (
        "Stage 1: generating the thermoelastic input",
        "Stage 2: preparing the QHA result",
        "Stage 3: calibrating the cold finite-strain model",
        "Stage 4: evaluating one state",
        "Stage 5: evaluating a pressure-temperature grid",
        "Stage 6: comparing isothermal and adiabatic tensors",
        "Stage 7: evaluating a geological profile",
        "Running the same workflow from Python",
        "Reproducibility checkpoints",
        "quantas thermoelasticity analysis point",
        "quantas thermoelasticity analysis grid",
        "quantas thermoelasticity analysis profile",
    ):
        assert phrase in text


def test_thermoelasticity_tutorial_script_uses_public_api_only() -> None:
    """The executable tutorial must not expose implementation packages."""
    path = Path("examples/thermoelasticity/tutorial_api.py")
    content = path.read_text(encoding="utf-8")
    assert "from quantas.api" in content
    assert "quantas.core" not in content
    assert "quantas.modules" not in content
    assert "quantas.renderers" not in content


def test_thermoelasticity_background_derives_the_quantas_qsa() -> None:
    """The public theory page distinguishes full QHA elasticity from QSA."""
    text = (DOCS_ROOT / "theory/thermoelasticity.rst").read_text(
        encoding="utf-8"
    )
    for phrase in (
        "Thermodynamic foundation",
        "Elastic coefficients under hydrostatic pre-stress",
        "Static and vibrational parts of the free energy",
        "Cold finite-strain elastic coefficients",
        "Explicit quasi-harmonic elastic contribution",
        "The quasi-static approximation",
        "Approximation implemented by Quantas",
        "stixrude_lithgow_bertelloni_2005",
    ):
        if phrase == "stixrude_lithgow_bertelloni_2005":
            assert "Stixrude and C. Lithgow-Bertelloni" in text
        else:
            assert phrase in text
    assert "C^{T,\\mathrm{QSA}}_{IJ}(P,T)" in text
    assert "\\Delta c^{\\mathrm q}_{ijkl}" in text


def test_eos_background_is_scientific_and_follows_angel_order() -> None:
    """EOS theory follows P-V, V-T, P-V-T and excludes solver mechanics."""
    text = (DOCS_ROOT / "theory/eos.rst").read_text(encoding="utf-8")
    headings = (
        "Pressure--volume equations of state",
        "Volume--temperature equations of state",
        "Pressure--volume--temperature equations of state",
    )
    positions = [text.index(heading) for heading in headings]
    assert positions == sorted(positions)
    assert "EosFit7c and a Fortran module" in text
    for phrase in (
        "Murnaghan equation",
        "Tait equation",
        "Birch--Murnaghan equation",
        "Natural-strain or Poirier--Tarantola equation",
        "Vinet equation",
        "Berman equation",
        "Fei equation",
        "Modified Holland--Powell equation",
        "Salje equation",
        "Kroll form of Holland--Powell",
        "Anderson--Gruneisen coupling",
        "Holland--Powell Einstein thermal pressure",
        "Mie--Gruneisen--Debye thermal pressure",
    ):
        assert phrase in text
    for forbidden in (
        "ordinary least squares",
        "weighted least squares",
        "orthogonal distance regression",
        "ParameterConstraint",
        "quantas.core",
        "quantas.modules",
    ):
        assert forbidden not in text


def test_scientific_background_bibliographies_are_canonical() -> None:
    """Theory citations use page-local numbers generated from the registry."""
    import importlib.util

    from quantas.references.registry import get_citation
    from quantas.references.render import render_rst_bibliography

    script = Path("docs/tools/generate_theory_bibliographies.py")
    spec = importlib.util.spec_from_file_location("theory_bibliographies", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    for page, keys in module.THEORY_REFERENCE_KEYS.items():
        theory_path = DOCS_ROOT / "theory" / f"{page}.rst"
        fragment_path = DOCS_ROOT / "_generated" / "references" / f"{page}.inc"
        text = theory_path.read_text(encoding="utf-8")
        inline_keys = tuple(dict.fromkeys(re.findall(r"\[#([A-Za-z0-9_]+)\]_", text)))
        assert inline_keys == keys
        assert fragment_path.read_text(encoding="utf-8") == render_rst_bibliography(keys)
        assert f".. include:: ../_generated/references/{page}.inc" in text
        for key in keys:
            citation = get_citation(key)
            assert f".. [#{key}]" in fragment_path.read_text(encoding="utf-8")
            if citation.doi:
                assert f"https://doi.org/{citation.doi}" in fragment_path.read_text(
                    encoding="utf-8"
                )


def test_theory_pages_do_not_embed_free_form_bibliographies() -> None:
    """Scientific citations remain linked to the canonical registry."""
    for path in (DOCS_ROOT / "theory").glob("*.rst"):
        text = path.read_text(encoding="utf-8")
        assert "doi:" not in text.lower()
        assert "https://doi.org/" not in text


def test_eos_tutorial_assets_and_downloads_exist() -> None:
    """The complete EOS tutorial retains its generated figures and files."""
    expected = [
        DOCS_ROOT / "_static/tutorials/eos/quartz_pv_fit.png",
        DOCS_ROOT / "_static/tutorials/eos/quartz_pv_residuals.png",
        DOCS_ROOT / "_static/tutorials/eos/quartz_pv_normalized_pressure.png",
        DOCS_ROOT / "_static/tutorials/eos/rutile_vt_fit.png",
        DOCS_ROOT / "_static/tutorials/eos/rutile_vt_residuals.png",
        DOCS_ROOT / "_static/tutorials/eos/naf_pvt_coverage.png",
        DOCS_ROOT / "_static/tutorials/eos/naf_pvt_isotherms.png",
        DOCS_ROOT / "_static/tutorials/eos/naf_pvt_isobars.png",
        DOCS_ROOT / "_static/tutorials/eos/naf_pvt_residuals_pressure.png",
        DOCS_ROOT / "_static/tutorials/eos/naf_pvt_residuals_temperature.png",
        DOCS_ROOT / "_downloads/tutorials/eos/quartz_pv_tutorial.spec",
        DOCS_ROOT / "_downloads/tutorials/eos/rutile_vt_tutorial.spec",
        DOCS_ROOT / "_downloads/tutorials/eos/naf_pvt_tutorial.spec",
        DOCS_ROOT / "_downloads/tutorials/eos/run_spec_api.py",
        DOCS_ROOT / "_downloads/tutorials/eos/pv_fit_api.py",
        DOCS_ROOT / "_downloads/tutorials/eos/vt_fit_api.py",
        DOCS_ROOT / "_downloads/tutorials/eos/pvt_fit_api.py",
        DOCS_ROOT / "_downloads/tutorials/eos/quartz_pv_comparison.csv",
        DOCS_ROOT / "_downloads/tutorials/eos/rutile_vt_comparison.csv",
        DOCS_ROOT / "_downloads/tutorials/eos/naf_pvt_coupling_comparison.csv",
    ]
    assert all(path.is_file() and path.stat().st_size > 0 for path in expected)


def test_eos_tutorial_follows_pv_vt_pvt_and_teaches_batch_fitting() -> None:
    """EOS tutorial pages remain detailed scientific CLI/API workflows."""
    landing = (DOCS_ROOT / "tutorials/eos.rst").read_text(encoding="utf-8")
    assert "Work in progress" not in landing
    ordered = ("eos_pv", "eos_vt", "eos_pvt")
    positions = [landing.index(name) for name in ordered]
    assert positions == sorted(positions)

    batch = (DOCS_ROOT / "tutorials/eos_batch.rst").read_text(encoding="utf-8")
    for phrase in (
        "Ordinary least squares",
        "Weighted least squares",
        "Effective variance",
        "Orthogonal distance regression",
        "Dry-run validation",
        "Running the same batch from Python",
        "Exercise: add ODR",
    ):
        assert phrase in batch

    expectations = {
        "eos_pv.rst": (
            "P--V tutorial: quartz compression",
            "Solver comparison",
            "Model-order comparison",
            "Normalized-pressure diagnostic",
            "Exercise 1: is BM4 justified?",
        ),
        "eos_vt.rst": (
            "V--T tutorial: thermal expansion of rutile",
            "Solver sensitivity",
            "Thermal-model sensitivity",
            "Axial thermal expansion",
            "Exercise 1: complete the ODR row",
        ),
        "eos_pvt.rst": (
            "P--V--T tutorial: coupled EOS of NaF",
            "Coupling comparison",
            "Pressure--temperature coverage",
            "Parameter constraints in the MGD job",
            "Exercise 1: release fewer parameters",
        ),
    }
    for filename, phrases in expectations.items():
        text = (DOCS_ROOT / "tutorials" / filename).read_text(encoding="utf-8")
        for phrase in phrases:
            assert phrase in text


def test_eos_tutorial_python_snippets_use_public_rendering_api() -> None:
    """EOS tutorials do not promote concrete renderer backends."""
    for filename in (
        "eos_batch.rst",
        "eos_pv.rst",
        "eos_vt.rst",
        "eos_pvt.rst",
        "eos_plotting.rst",
        "eos_mgd_fit_api.rst",
    ):
        text = (DOCS_ROOT / "tutorials" / filename).read_text(encoding="utf-8")
        assert "quantas.renderers" not in text
        assert "render_plot_collection" not in text
        assert "MatplotlibOptions(" not in text


def test_ha_qha_workflow_pages_are_complete() -> None:
    """HA/QHA implementation pages retain the important numerical decisions."""
    ha = (DOCS_ROOT / "workflows/ha.rst").read_text(encoding="utf-8")
    qha = (DOCS_ROOT / "workflows/qha.rst").read_text(encoding="utf-8")

    for phrase in (
        "Multi-volume HA",
        "q-point weights",
        "Treatment of zero and imaginary frequencies",
        "Performance and memory use",
        "no numerical parameter with a default value of",
    ):
        assert phrase in ha

    for phrase in (
        "Preflight inspection",
        "Frequency and thermodynamic schemes",
        "Polynomial minimization",
        "EOS minimization",
        "Polynomial thermoelastic derivatives",
        "Three routes to volumetric thermal expansion",
        "Failure policies and partial results",
        "Performance and practical acceleration",
        "There is no HA/QHA numerical control with a default value of",
    ):
        assert phrase in qha

    assert "**5 points**" in qha
    assert "**0.05%**" in qha
    assert "mixed derivative — default".lower() in qha.lower()


def test_elasticity_seismic_workflow_pages_are_complete() -> None:
    """Elasticity/SEISMIC implementation pages retain numerical decisions."""
    elasticity = (DOCS_ROOT / "workflows/elasticity.rst").read_text(encoding="utf-8")
    seismic = (DOCS_ROOT / "workflows/seismic.rst").read_text(encoding="utf-8")

    for phrase in (
        "Global directional extrema",
        "Principal-plane 2D fields",
        "Three-dimensional directional fields",
        "Persisted and transient surfaces",
        "Surface batch size",
        "65536",
        "The value ``512`` is **not** an Elasticity",
        "Suggested 3D convergence check",
    ):
        assert phrase in elasticity

    for phrase in (
        "Sampling levels",
        "Christoffel eigenvalue validity",
        "Degeneracy detection",
        "Polarization tracking",
        "Analytical group quantities",
        "Enhancement, pseudoinverse, and caustic candidates",
        "The meaning of ``batch_size=512``",
        "Batch size is not an accuracy or convergence parameter",
        "Angular convergence",
        "16471",
    ):
        assert phrase in seismic


def test_thermoelasticity_workflow_page_is_complete() -> None:
    """Thermoelastic implementation page retains its scientific decisions."""
    text = (DOCS_ROOT / "workflows/thermoelasticity.rst").read_text(
        encoding="utf-8"
    )

    for phrase in (
        "Why CRYSTAL ``PRESSURE`` is required",
        "Reference point and frame normalization",
        "What is and is not fitted from the elastic outputs",
        "Cold reference EOS",
        "Choosing BM2, BM3, or BM4",
        "Second- and third-order finite strain",
        "Scientific support diagnostics",
        "Validation presets",
        "Shared reference-EOS covariance",
        "Two independent extrapolation masks",
        "Isothermal-to-adiabatic conversion",
        "Recommended sensitivity study",
        "Performance and acceleration",
        "No ``512`` convergence parameter",
        "Common interpretation errors",
    ):
        assert phrase in text

    for value in (
        "0.05 GPa",
        "0.5 Å",
        "1.0e-10 GPa",
        "0.005",
        ":math:`10^6`",
        "piecewise-linear interpolation",
        "C^S=C^T",
    ):
        assert value in text


def test_command_reference_pages_are_contextual_and_generated() -> None:
    """CLI pages combine maintained guidance with sphinx-click output."""
    conventions = (DOCS_ROOT / "cli/conventions.rst").read_text(encoding="utf-8")
    for phrase in (
        "Reading a command signature",
        "Option groups",
        "Paths and overwrite policy",
        "Reports, terminal output, and progress",
        "Plot presets",
        "Errors and exit status",
    ):
        assert phrase in conventions

    pages = {
        "ha": "quantas.cli.ha:ha",
        "qha": "quantas.cli.qha:qha",
        "elasticity": "quantas.cli.elasticity:elasticity",
        "seismic": "quantas.cli.seismic:seismic",
        "eos": "quantas.cli.eos:eos",
        "thermoelasticity": "quantas.cli.thermoelastic:thermoelasticity",
    }
    for name, target in pages.items():
        text = (DOCS_ROOT / "cli" / f"{name}.rst").read_text(encoding="utf-8")
        assert "Generated command reference" in text
        assert f".. click:: {target}" in text
        assert ":nested: full" in text
        assert "Recommended sequence" in text
        assert f"../workflows/{name}" in text
        assert f"../tutorials/{name}" in text


def test_command_reference_navigation_includes_shared_conventions() -> None:
    """The shared CLI conventions page is part of the command-reference tree."""
    text = (DOCS_ROOT / "index.rst").read_text(encoding="utf-8")
    command_section = text.split(":caption: COMMAND REFERENCE", maxsplit=1)[1]
    command_section = command_section.split(".. toctree::", maxsplit=1)[0]
    assert "cli/conventions" in command_section


def test_getting_started_pages_are_complete_and_public_api_only() -> None:
    """Getting Started remains a complete, runnable introduction."""
    expectations = {
        "first_calculation.rst": (
            "Run the command-line calculation",
            "Repeat the calculation through Python",
            "Compare the two frontends",
            "83.099307",
        ),
        "python_api.rst": (
            "The single-shot workflow pattern",
            "EOS uses an archive lifecycle",
            "Cross-module transformations",
            "Capability discovery",
        ),
        "command_line.rst": (
            "Discover the command tree",
            "The common command pattern",
            "Command groups are scientifically different",
            "Errors and exit status",
        ),
        "results.rst": (
            "Result envelopes",
            "Native HDF5 results",
            "Tabular exports",
            "What should be kept?",
        ),
    }
    for filename, phrases in expectations.items():
        text = (DOCS_ROOT / "getting_started" / filename).read_text(
            encoding="utf-8"
        )
        assert "Work in progress" not in text
        assert "quantas.models" not in text
        for phrase in phrases:
            assert phrase in text

    script = DOCS_ROOT / "_downloads/getting_started/first_elasticity_api.py"
    script_text = script.read_text(encoding="utf-8")
    assert script.is_file()
    assert "from quantas.api import elasticity, rendering" in script_text
    assert "quantas.core" not in script_text
    assert "quantas.modules" not in script_text
    assert "quantas.renderers" not in script_text


def test_rtd_theme_layout_is_full_width_and_independent() -> None:
    """The desktop manual uses the available column without linked nav scroll."""
    conf = (DOCS_ROOT / "conf.py").read_text(encoding="utf-8")
    css = (DOCS_ROOT / "_static/custom.css").read_text(encoding="utf-8")
    assert '"sticky_navigation": False' in conf
    assert ".wy-nav-content" in css
    assert "max-width: none" in css
    assert ".wy-side-scroll" in css
    assert "overflow-y: auto" in css
    assert "overscroll-behavior: contain" in css


def test_rst_title_underlines_are_not_shorter_than_titles() -> None:
    """Simple reST section adornments do not trigger docutils warnings."""
    adornments = set('=-~^"`:+*#<>_')
    problems: list[str] = []
    for path in DOCS_ROOT.rglob("*.rst"):
        lines = path.read_text(encoding="utf-8").splitlines()
        for index in range(1, len(lines)):
            underline = lines[index]
            title = lines[index - 1]
            stripped = underline.strip()
            if (
                not stripped
                or underline != stripped
                or len(set(stripped)) != 1
                or stripped[0] not in adornments
                or not title.strip()
                or title.startswith("..")
            ):
                continue
            if len(stripped) < len(title):
                problems.append(
                    f"{path}:{index + 1}: {len(stripped)} < {len(title)}: {title}"
                )
    assert problems == []


def test_input_output_formats_are_complete_and_navigated() -> None:
    """The format reference covers every public input and native payload."""
    root = (DOCS_ROOT / "index.rst").read_text(encoding="utf-8")
    pages = {
        "elasticity_input": ("Voigt convention", "hydrostatic pre-stress"),
        "phonon_yaml": ("Q-point weights", "Mode continuity"),
        "eos_input": ("P--V--T tables", "CLI unit overrides"),
        "eos_spec": ("Default precedence", "failure_policy"),
        "thermoelastic_input": ("Top-level structure", "Frame metadata"),
        "earth_profile_spec": ("Composed YAML schema", "Piecewise temperature profiles"),
        "hdf5": ("Common envelope", "Numerical precision"),
        "ha_qha_hdf5": ("HA payload", "QHA payload"),
        "elasticity_seismic_hdf5": ("Elasticity payload", "SEISMIC payload"),
        "eos_hdf5": ("Immutable fit records", "Batch manifest"),
        "thermoelastic_hdf5": ("Component fits", "Depth-profile payloads"),
        "hdf5_inspection": ("Metadata-driven reading with Quantas", "Extracting arrays to NPZ"),
        "tabular_outputs": ("SEISMIC long-form CSV", "Recommended archival set"),
    }
    for name, phrases in pages.items():
        path = DOCS_ROOT / "formats" / f"{name}.rst"
        text = path.read_text(encoding="utf-8")
        assert f"formats/{name}" in root
        assert "Work in progress" not in text
        for phrase in phrases:
            assert phrase in text


def test_hdf5_downloads_are_standalone_scripts() -> None:
    """Distributed HDF5 utilities remain executable and avoid private imports."""
    pairs = (
        (
            Path("examples/io/inspect_hdf5.py"),
            DOCS_ROOT / "_downloads/formats/inspect_hdf5.py",
        ),
        (
            Path("examples/io/extract_hdf5.py"),
            DOCS_ROOT / "_downloads/formats/extract_hdf5.py",
        ),
    )
    for example, download in pairs:
        assert example.is_file() and download.is_file()
        assert example.read_bytes() == download.read_bytes()
        text = example.read_text(encoding="utf-8")
        assert "quantas.core" not in text
        assert "quantas.modules" not in text
        assert "quantas.renderers" not in text
        assert "import h5py" in text


def test_hdf5_guide_preserves_data_context() -> None:
    """Custom-analysis guidance must not reduce native data to bare arrays."""
    text = (DOCS_ROOT / "formats/hdf5_inspection.rst").read_text(
        encoding="utf-8"
    )
    for phrase in (
        "Validity and support masks",
        "source HDF5 filename and SHA-256 checksum",
        "exact HDF5 dataset paths",
        "original units and any conversions",
        "Do not overwrite the original native HDF5 file",
        "EOS archive example",
    ):
        assert phrase in text


def test_format_pages_use_only_public_quantas_namespaces() -> None:
    """Format documentation does not teach imports from implementation layers."""
    forbidden = (
        "quantas.core",
        "quantas.modules",
        "quantas.models",
        "quantas.renderers",
        "quantas.io",
    )
    for path in (DOCS_ROOT / "formats").glob("*.rst"):
        text = path.read_text(encoding="utf-8")
        for name in forbidden:
            assert name not in text, f"{path}: public format page exposes {name}"
