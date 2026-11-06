"""Export clean numeric benchmark tables for Living Lab Book PGFPlots figures."""

from __future__ import annotations

from pathlib import Path

from scripts.m1_6._common import write_dat_benchmark


def generate_plot_tables() -> None:
    # 1. Critical zoom plot table: h, E_fwd, E_bwd, W_fwd, W_bwd, mu_fwd
    zoom_raw = Path("data/processed/benchmarks/m1_6/critical_zoom_111.dat")
    if zoom_raw.is_file():
        fwd_dict = {}
        bwd_dict = {}
        with zoom_raw.open("r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.split()
                h = float(p[0])
                branch = p[1]
                e = float(p[5])
                w = float(p[7])
                mu = float(p[9]) if p[9] != "nan" else 0.0
                if branch == "fwd_qsl":
                    fwd_dict[h] = (e, w, mu)
                elif branch == "bwd_pol":
                    bwd_dict[h] = (e, w, mu)

        common_h = sorted(set(fwd_dict.keys()) & set(bwd_dict.keys()))
        data_zoom = []
        for h in common_h:
            e_f, w_f, mu_f = fwd_dict[h]
            e_b, w_b, _ = bwd_dict[h]
            data_zoom.append(
                f"{h:6.4f}  {e_f:12.6f}  {e_b:12.6f}  {w_f:8.4f}  {w_b:8.4f}  {mu_f:9.4f}"
            )

        headers_zoom = [
            "Description: Clean numeric table for Lab Book Figure 10 (Critical Zoom)",
            "Columns: h  E_fwd  E_bwd  W_fwd  W_bwd  mu_fwd",
        ]
        write_dat_benchmark("critical_zoom_plot.dat", headers_zoom, data_zoom)

    # 2. Majorana gap plot table: h, gap_K, gap_global, gap_cubic_ref
    gap_raw = Path("data/processed/benchmarks/m1_6/majorana_gap_scaling.dat")
    if gap_raw.is_file():
        data_gap = []
        with gap_raw.open("r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.split()
                h = float(p[0])
                gk = float(p[1])
                gg = float(p[2])
                ref = 2.5921 * (h**3)
                data_gap.append(f"{h:6.3f}  {gk:12.6e}  {gg:12.6e}  {ref:12.6e}")

        headers_gap = [
            "Description: Clean numeric table for Lab Book Figure 11 (Majorana Gap Scaling)",
            "Columns: h  gap_K  gap_global  ref_cubic",
        ]
        write_dat_benchmark("majorana_gap_plot.dat", headers_gap, data_gap)

    # 3. Supercell comparison plot table: h, E_2site_qsl, E_4site_qsl, E_2site_pol, E_4site_pol
    super_raw = Path("data/processed/benchmarks/m1_6/supercell_comparison.dat")
    if super_raw.is_file():
        records: dict[float, dict[tuple[str, str], tuple[float, float]]] = {}
        with super_raw.open("r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.split()
                h = float(p[0])
                cell = p[1]
                seed = p[2]
                e = float(p[7])
                nem = float(p[9])
                records.setdefault(h, {})[(cell, seed)] = (e, nem)

        data_super = []
        for h in sorted(records.keys()):
            e_2q, n_2q = records[h].get(("2site", "qsl"), (float("nan"), 0.0))
            e_4q, n_4q = records[h].get(("4site", "qsl"), (float("nan"), 0.0))
            e_2p, _ = records[h].get(("2site", "pol"), (float("nan"), 0.0))
            e_4p, _ = records[h].get(("4site", "pol"), (float("nan"), 0.0))
            data_super.append(
                f"{h:6.3f}  {e_2q:12.6f}  {e_4q:12.6f}  {e_2p:12.6f}  {e_4p:12.6f}  {n_4q:12.4e}"
            )

        headers_super = [
            "Description: Clean numeric table for Lab Book Figure 12 (Supercell Parity)",
            "Columns: h  E_2q  E_4q  E_2p  E_4p  nem_4q",
        ]
        write_dat_benchmark("supercell_plot.dat", headers_super, data_super)


if __name__ == "__main__":
    generate_plot_tables()
