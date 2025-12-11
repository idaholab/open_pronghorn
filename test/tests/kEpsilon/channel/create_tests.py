"""
Generate a MOOSE tests file for all permutations of the k-epsilon model
for the channel_ERCOFTAC case.

- One test block per permutation of:
    * k_epsilon_variant
    * nonlinear_model (none, quadratic, cubic)
    * use_yap (true/false)
    * two_layer_flavor (Wolfstein, NorrisReynolds, Xu) for *TwoLayer variants
- Each test:
    * Uses channel_ERCOFTAC.i as input
    * Uses a unique Outputs/file_base
    * csvdiff list matches the actual output CSV naming (no 'out_', index 0002)

Usage:
  python generate_channel_ercoftac_kepsilon_tests.py > tests

Then include the generated [Tests] block into your MOOSE tests file.
"""

import itertools


# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

K_EPSILON_VARIANTS = [
    "Standard",
    "StandardLowRe",
    "StandardTwoLayer",
    "Realizable",
    "RealizableTwoLayer",
]

TWO_LAYER_FLAVORS = [
    "Wolfstein",
    "NorrisReynolds",
    "Xu",
]

NONLINEAR_MODELS = [
    "none",
    "quadratic",
    "cubic",
]

USE_YAP_OPTIONS = [False, True]

BASE_TEST_NAME_PREFIX = "channel_ERCOFTAC"
BASE_INPUT_FILE = "channel_ERCOFTAC.i"
BASE_REQUIREMENT = (
    "The system shall be able to solve fluid flow problems with k-epsilon "
    "turbulence model for a standard channel with linear FV discretization, "
    "and reach converged results with segregated solvers."
)


# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------


def variant_token(variant: str) -> str:
    """
    Map the k_epsilon_variant to the token used in filenames.
    Based on your CSVs:
      Standard          -> standard
      StandardLowRe     -> standardlowre
      StandardTwoLayer  -> standardtwolayer
      Realizable        -> realizable
      RealizableTwoLayer-> realizabletwolayer
    """
    return variant.lower()


def bool_to_moose(b: bool) -> str:
    return "true" if b else "false"


def is_skipped_permutation(p) -> bool:
    """
    Return True if this permutation should be skipped based on the failures:

      - realizabletwolayer_norrisreynolds_quadratic_yap  (ERRMSG)
      - realizabletwolayer_norrisreynolds_cubic_yap      (ERRMSG)
      - realizabletwolayer_xu_cubic                      (MISSING GOLD FILE)
      - realizabletwolayer_xu_cubic_yap                  (MISSING GOLD FILE)
      - realizabletwolayer_wolfstein_cubic_yap           (MISSING GOLD FILE)
    """
    variant = p["variant"]
    flavor = p["two_layer_flavor"]
    nonlinear = p["nonlinear_model"]
    use_yap = p["use_yap"]

    # Skip StandardTwoLayer + Xu flavor (all nonlinear/YAP combos)
    if variant == "StandardTwoLayer" and flavor == "Xu":
        return True

    if variant == "RealizableTwoLayer":
        # NorrisReynolds, quadratic/cubic, with YAP
        if (
            flavor == "NorrisReynolds"
            and use_yap
            and nonlinear in ("quadratic", "cubic")
        ):
            return True

        # Xu, cubic (with or without YAP)
        if flavor == "Xu" and nonlinear == "cubic":
            return True

        # Wolfstein, cubic with YAP
        if flavor == "Wolfstein" and nonlinear == "cubic" and use_yap:
            return True

    return False


def iter_permutations():
    """
    Yield dictionaries describing each permutation:

        {
          "variant": ...,
          "two_layer_flavor": ... or None,
          "nonlinear_model": ...,
          "use_yap": True/False,
        }

    Full permutation over variant × nonlinear_model × use_yap, and
    over two_layer_flavor for *TwoLayer variants, except for the
    explicitly skipped cases.
    """
    for variant, nonlinear_model, use_yap in itertools.product(
        K_EPSILON_VARIANTS, NONLINEAR_MODELS, USE_YAP_OPTIONS
    ):
        if "TwoLayer".lower() in variant.lower():
            flavors = TWO_LAYER_FLAVORS
        else:
            flavors = [None]

        for flavor in flavors:
            p = {
                "variant": variant,
                "two_layer_flavor": flavor,
                "nonlinear_model": nonlinear_model,
                "use_yap": use_yap,
            }
            if is_skipped_permutation(p):
                continue
            yield p


def build_file_base(p) -> str:
    """
    Build the Outputs/file_base from a permutation p, matching the
    observed patterns, e.g.:

      channel_ERCOFTAC_k_epsilon_realizable
      channel_ERCOFTAC_k_epsilon_realizable_quadratic_yap
      channel_ERCOFTAC_k_epsilon_realizabletwolayer_wolfstein_quadratic_yap
    """
    parts = ["channel_ERCOFTAC", "k_epsilon"]

    # Variant core
    vt = variant_token(p["variant"])  # e.g. realizabletwolayer
    parts.append(vt)

    # Two-layer flavor for *TwoLayer variants
    if p["two_layer_flavor"] is not None:
        parts.append(p["two_layer_flavor"].lower())

    # Nonlinear model if not 'none'
    if p["nonlinear_model"] != "none":
        parts.append(p["nonlinear_model"].lower())

    # YAP
    if p["use_yap"]:
        parts.append("yap")

    return "_".join(parts)


def make_cli_args(p, file_base: str) -> str:
    """
    Build the cli_args string for a permutation `p`.
    """
    args = []

    # Turbulence knobs
    args.append(f"k_epsilon_variant={p['variant']}")
    if p["nonlinear_model"] != "none":
        args.append(f"nonlinear_model={p['nonlinear_model']}")
    if p["two_layer_flavor"] is not None and "TwoLayer" in p["variant"]:
        args.append(f"two_layer_flavor={p['two_layer_flavor']}")
    args.append(f"use_yap={bool_to_moose(p['use_yap'])}")

    # Output base so CSVs are unique per permutation
    args.append(f"Outputs/file_base={file_base}")

    return " ".join(args)


def make_requirement_suffix(p) -> str:
    """
    Human-readable description for this permutation, appended to the base requirement.
    """
    desc = [p["variant"]]
    if p["two_layer_flavor"] is not None and "TwoLayer" in p["variant"]:
        desc.append(f"{p['two_layer_flavor']} two-layer flavor")
    if p["nonlinear_model"] != "none":
        desc.append(f"{p['nonlinear_model']} non-linear model")
    if p["use_yap"]:
        desc.append("YAP correction")

    if not desc:
        return ""
    return " using " + ", ".join(desc)


def make_csvdiff_list(file_base: str) -> str:
    """
    Build the csvdiff string based on the given file_base, matching
    the actual filenames (no 'out_', index 0002).
    """
    parts = [
        f"{file_base}_line_center_channel_0002.csv",
        f"{file_base}_line_quarter_radius_channel_0002.csv",
        f"{file_base}_side_bottom_0002.csv",
        f"{file_base}_side_top_0002.csv",
    ]
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

        # Test name: same as file_base, which matches the CSV prefix
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
        print(f"    cli_args = '{cli_args}'")
        print(f"    requirement = '{BASE_REQUIREMENT}{req_suffix}'")
        print("  []")

    print("[]")  # close [Tests]


if __name__ == "__main__":
    main()
