"""Print a reproducible K-006 midpoint-quadrature convergence table."""

import argparse

from kitaev_mf.services.quadrature import midpoint_convergence


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--meshes",
        nargs="+",
        type=int,
        default=[96, 192, 384, 768, 1536, 3072, 6144, 12288],
    )
    parser.add_argument("--row-chunk-size", type=int, default=128)
    args = parser.parse_args()

    estimates = midpoint_convergence(args.meshes, args.row_chunk_size)
    print("N,value,observed_order,estimated_error,richardson_value")
    for item in estimates:
        order = "" if item.observed_order is None else f"{item.observed_order:.8f}"
        error = "" if item.estimated_error is None else f"{item.estimated_error:.6e}"
        richardson = (
            "" if item.richardson_value is None else f"{item.richardson_value:.15f}"
        )
        print(
            f"{item.mesh_size},{item.value:.15f},"
            f"{order},{error},{richardson}"
        )


if __name__ == "__main__":
    main()
