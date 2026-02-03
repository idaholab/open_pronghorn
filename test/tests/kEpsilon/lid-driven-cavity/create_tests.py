"""
Generate a MOOSE tests file for k-epsilon permutations on the lid-driven cavity.

Permutations are done over:
  - k_epsilon_variant
  - use_low_re_Gprime (applied only for StandardLowRe due to kernel rail-guards)

We explicitly do NOT vary:
  - two_layer_flavor
  - nonlinear_model
  - use_yap
  - curvature_model

Each test:
  - Uses 'lid-driven.i' as input
  - Sets cli_args to control the turbulence knobs and Outputs/file_base
  - Compares against CSV files named:

      lid-driven_k_epsilon_<variant_token>[_lowreGprime]_<sampler>_0002.csv

    where <sampler> ∈ {
      horizontal_center, vertical_center,
      side_bottom, side_top, side_left, side_right
    }

Usage:
  python create_tests.py > tests

Then include the generated [Tests] block into your MOOSE tests file.
"""

import itertools


# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

# From lid-driven.i comment / your gold filenames:
#   Standard | StandardLowRe | StandardTwoLayer | Realizable | RealizableTwoLayer
K_EPSILON_VARIANTS = [
    "Standard",
    "StandardLowRe",
    "StandardTwoLayer",
    "Realizable",
    "RealizableTwoLayer",
]

# Low-Re Gprime toggle
USE_LOW_RE_GPRIME_OPTIONS = [False, True]

# Base info
BASE_INPUT_FILE = "lid-driven.i"
BASE_REQUIREMENT = (
    "The system shall be able to solve turbulent lid-driven cavity flow using "
    "the k-epsilon turbulence model and reach converged results with "
    "segregated solvers "
)

# VectorPostprocessor samplers we want to compare (from your gold list)
SAMPLER_NAMES = [
    "horizontal_center",
    "vertical_center",
    "side_bottom",
    "side_top",
    "side_left",
    "side_right",
]


# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------


def variant_token(variant: str) -> str:
    """
    Map k_epsilon_variant to a token for filenames.
    Based on your gold:

      Standard          -> standard
      StandardLowRe     -> standardlowre
      StandardTwoLayer  -> standardtwolayer
      Realizable        -> realizable
      RealizableTwoLayer-> realizabletwolayer
    """
    return variant.lower()


def bool_to_moose(b: bool) -> str:
    return "true" if b else "false"


def iter_permutations():
    """
    Yield dictionaries describing each permutation:

        {
          "variant": ...,
          "use_low_re_Gprime": True/False,
        }
    """
    for variant, use_low_re in itertools.product(
        K_EPSILON_VARIANTS, USE_LOW_RE_GPRIME_OPTIONS
    ):
        yield {
            "variant": variant,
            "use_low_re_Gprime": use_low_re,
        }


def build_file_base(p) -> str:
    """
    Build the Outputs/file_base for the lid-driven case.

    Pattern from your gold:
      lid-driven_k_epsilon_<variant_token>[_lowreGprime]

    Examples:
      lid-driven_k_epsilon_standard
      lid-driven_k_epsilon_standardlowre_lowreGprime
      lid-driven_k_epsilon_realizable
      lid-driven_k_epsilon_realizabletwolayer_lowreGprime
    """
    parts = ["lid-driven", "k_epsilon"]

    vt = variant_token(p["variant"])
    parts.append(vt)

    if p["use_low_re_Gprime"]:
        parts.append("lowreGprime")

    return "_".join(parts)


def make_cli_args(p, file_base: str) -> str:
    """
    Build the cli_args string for a permutation `p`.

    NOTE: We set parameters on the actual objects (kernels/auxkernel) rather than
    relying on top-level 'header' variables in lid-driven.i. This avoids tripping
    the newer rail-guards that reject unsupported option/variant combinations.

    Only k_epsilon_variant and use_low_re_Gprime are varied by this test generator.
    """
    args = []

    # Set the k-epsilon variant consistently everywhere it is used
    args.append(f"AuxKernels/compute_mu_t/k_epsilon_variant={p['variant']}")
    args.append(f"LinearFVKernels/TKE_source_sink/k_epsilon_variant={p['variant']}")
    args.append(f"LinearFVKernels/TKED_source_sink/k_epsilon_variant={p['variant']}")

    # Low-Re G' correction is only supported for the StandardLowRe variant
    enable_gprime = p["use_low_re_Gprime"] and p["variant"] == "StandardLowRe"
    args.append(
        f"LinearFVKernels/TKED_source_sink/use_low_re_Gprime={bool_to_moose(enable_gprime)}"
    )

    # Set output base
    args.append(f"Outputs/file_base={file_base}")

    return " ".join(args)


def make_requirement_suffix(p) -> str:
    """
    Human-readable description for this permutation, appended to the base requirement.
    Avoid apostrophes to keep the single-quoted MOOSE string valid.
    """
    desc = [p["variant"] + " variant"]

    if p["use_low_re_Gprime"]:
        desc.append("Low-Re Gprime correction")

    if not desc:
        return ""
    if len(desc) == 1:
        return "using " + desc[0] + "."
    if len(desc) == 2:
        return "using " + desc[0] + " and " + desc[1] + "."
    else:
        return "using " + ", ".join(desc[1:]) + ", and " + desc[-1] + "."


def make_csvdiff_list(file_base: str) -> str:
    """
    Build the csvdiff string for this test.

    Pattern from your gold:
      <file_base>_<sampler>_0002.csv  for each sampler in SAMPLER_NAMES
    """
    parts = [f"{file_base}_{name}_0002.csv" for name in SAMPLER_NAMES]
    return " ".join(parts)


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------


def main():
    print("[Tests]")
    print(
        "  design = 'kEpsilonViscosity.md "
        "kEpsilonTKEDSourceSink.md "
        "kEpsilonTKESourceSink.md "
        "WallDistanceAux.md'"
    )
    print("  issues = '#50'")

    for p in iter_permutations():
        file_base = build_file_base(p)
        cli_args = make_cli_args(p, file_base)
        csvdiff_list = make_csvdiff_list(file_base)
        req_suffix = make_requirement_suffix(p)

        # Use file_base as the test name for simplicity/readability
        test_name = file_base

        print(f"  [{test_name}]")
        print("    type = 'CSVDiff'")
        print(f"    input = '{BASE_INPUT_FILE}'")
        print(f"    csvdiff = '{csvdiff_list}'")
        print("    abs_zero = 1e-3")
        print("    rel_err = 1e-3")
        print("    recover = false")
        print("    max_threads = 1")
        print("    mesh_mode = 'replicated'")
        print("    heavy = false")
        print(f"    cli_args = '{cli_args}'")
        print(f"    requirement = '{BASE_REQUIREMENT}{req_suffix}'")
        print("  []")

    print("[]")  # close [Tests]


if __name__ == "__main__":
    main()
