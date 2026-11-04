#!/usr/bin/env python3
"""Print a reproducible K-006 midpoint-quadrature convergence table."""

import click

from kitaev_mf.services.quadrature import midpoint_convergence


@click.command()
@click.option(
    "--meshes",
    default="96,192,384,768,1536,3072,6144,12288",
    help="Comma-separated doubling sequence of mesh sizes (multiples of 3).",
)
@click.option(
    "--row-chunk-size",
    default=128,
    type=int,
    help="Row chunk size to control memory footprint.",
)
def main(meshes: str, row_chunk_size: int) -> None:
    mesh_sizes = tuple(int(m.strip()) for m in meshes.split(",") if m.strip())
    estimates = midpoint_convergence(mesh_sizes, row_chunk_size=row_chunk_size)

    header = (
        f"{'Mesh':>7} | {'U_bar':>18} | {'Order p':>8} | "
        f"{'Est. Error':>12} | {'Richardson U_R':>18}"
    )
    print(header)
    print("-" * len(header))

    for est in estimates:
        order_str = f"{est.observed_order:8.3f}" if est.observed_order is not None else f"{'--':>8}"
        error_str = (
            f"{est.estimated_error:12.3e}" if est.estimated_error is not None else f"{'--':>12}"
        )
        rich_str = (
            f"{est.richardson_value:18.13f}" if est.richardson_value is not None else f"{'--':>18}"
        )
        print(f"{est.mesh_size:7d} | {est.value:18.13f} | {order_str} | {error_str} | {rich_str}")


if __name__ == "__main__":
    main()
