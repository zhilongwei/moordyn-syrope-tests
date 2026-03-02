import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from typing import Any, Dict
from math import isclose

from syropepy import Syrope

# rc('text', usetex=True) # Requires LaTeX installation, so not enabled by default
rc('font', family='serif', size=11)

kG = 9.80665
kMBL = 885.0e3 * kG  # Minimul Breaking Load in N

Tdur = 3 * 60 * 60   # Duration of one sea state in seconds (3 hours)

Tmax0  = 0.16 * kMBL  # Initial Tmax guess (16% of MBL)
Tmean0 = 0.05 * kMBL  # Initial Tmean guess (6% of MBL)

phase_list = ['I', 'II', 'III', 'IV', 'V', 'VI']

# Parse Syrope inputs
def parse_syrope_inputs(file_path):
    """Parse the SYROPE-related bits from a MoorDyn input file.

    Extracts:
    - SYROPE `owc.dat` reference
    - Dynamic stiffness parameters: a, b
    - Damping parameters: C2, C1
    - Working curve: WCType, k1, k2
    - Initial positions of points
    - SYROPE IC (Tmax, Tmean)

    Returns a dict with keys:
    - line_types, working_curves, points, lines, syrope_ic (raw section data)
    - line_summary: compact per-line extracted values
    """

    def _as_float(value: str) -> float:
        return float(value.replace('D', 'e').replace('d', 'e'))

    def _is_data_row(line: str) -> bool:
        s = line.strip()
        if not s:
            return False
        if s.startswith('-'):
            return False
        if set(s) <= {'-'}:
            return False
        if s.startswith('(-)'):
            return False
        return True

    current_section = None
    results: Dict[str, Any] = {
        'line_types': {},
        'working_curves': {},
        'points': {},
        'lines': {},
        'syrope_ic': {},
        'line_summary': {},
    }

    with open(file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    for raw in lines:
        line = raw.strip('\n')
        upper = line.upper()

        if 'LINE TYPES' in upper and 'SYROPE WORKING CURVES' not in upper:
            current_section = 'line_types'
            continue
        if 'SYROPE WORKING CURVES' in upper:
            current_section = 'working_curves'
            continue
        if upper.strip().startswith('----------------------- POINTS') or upper.strip().startswith('-----------------------POINTS') or ' POINTS ' in upper:
            if 'POINTS' in upper and 'LINES' not in upper:
                current_section = 'points'
                continue
        if upper.strip().startswith('-------------------------- LINES') or upper.strip().startswith('--------------------------LINES'):
            current_section = 'lines'
            continue
        if 'SYROPE IC' in upper:
            current_section = 'syrope_ic'
            continue
        if 'SOLVER OPTIONS' in upper:
            current_section = None
            continue

        if not _is_data_row(line):
            continue

        tokens = line.split()
        if current_section == 'line_types':
            # Expected columns:
            # LineType Diam MassDenInAir EA BA/-zeta EI Can Cat Cdn Cdt
            if not tokens or tokens[0].lower() in {'linetype', 'line', 'node'}:
                continue
            if len(tokens) < 5:
                continue

            line_type = tokens[0]
            ea_token = tokens[3]
            ba_token = tokens[4]

            if ea_token.startswith('SYROPE:'):
                spec = ea_token[len('SYROPE:'):]
                parts = spec.split('|')
                if len(parts) >= 3:
                    owc_path = parts[0]
                    a = _as_float(parts[1])
                    b = _as_float(parts[2])
                else:
                    owc_path = parts[0] if parts else ''
                    a = float('nan')
                    b = float('nan')

                damp_parts = ba_token.split('|')
                if len(damp_parts) >= 2:
                    C2 = _as_float(damp_parts[0])
                    C1 = _as_float(damp_parts[1])
                else:
                    C2 = float('nan')
                    C1 = float('nan')

                results['line_types'][line_type] = {
                    'owc_path': owc_path,
                    'a': a,
                    'b': b,
                    'C2': C2,
                    'C1': C1,
                }

        elif current_section == 'working_curves':
            # Expected columns: LineType WCType k1 k2
            if not tokens or tokens[0].lower() in {'linetype', 'wctype'}:
                continue
            if len(tokens) < 4:
                continue
            line_type = tokens[0]
            wc_type = tokens[1]
            k1 = _as_float(tokens[2])
            k2 = _as_float(tokens[3])
            results['working_curves'][line_type] = {
                'WCType': wc_type,
                'k1': k1,
                'k2': k2,
            }

        elif current_section == 'points':
            # Expected columns:
            # Node Type X Y Z M V CdA CA
            if not tokens or tokens[0].lower() in {'node', 'line'}:
                continue
            if len(tokens) < 5:
                continue
            node = int(tokens[0])
            ptype = tokens[1]
            x = _as_float(tokens[2])
            y = _as_float(tokens[3])
            z = _as_float(tokens[4])
            results['points'][node] = {
                'type': ptype,
                'xyz': (x, y, z),
            }

        elif current_section == 'lines':
            # Expected columns:
            # Line LineType NodeA NodeB UnstrLen NumSegs Flags/Outputs
            if not tokens or tokens[0].lower() in {'line', 'nodea'}:
                continue
            if len(tokens) < 6:
                continue
            line_id = int(tokens[0])
            line_type = tokens[1]
            node_a = int(tokens[2])
            node_b = int(tokens[3])
            unstr_len = _as_float(tokens[4])
            num_segs = int(tokens[5])
            flags = tokens[6] if len(tokens) >= 7 else ''
            results['lines'][line_id] = {
                'line_type': line_type,
                'nodeA': node_a,
                'nodeB': node_b,
                'unstrLen': unstr_len,
                'numSegs': num_segs,
                'flags': flags,
            }

        elif current_section == 'syrope_ic':
            # Expected columns:
            # Line Tmax Tmean
            if not tokens or tokens[0].lower() in {'line', 'tmax'}:
                continue
            if len(tokens) < 3:
                continue
            try:
                line_id = int(tokens[0])
                tmax = _as_float(tokens[1])
                tmean = _as_float(tokens[2])
            except ValueError:
                continue
            results['syrope_ic'][line_id] = {
                'Tmax': tmax,
                'Tmean': tmean,
            }

    for line_id, line_info in results['lines'].items():
        line_type = line_info['line_type']
        node_a = line_info['nodeA']
        node_b = line_info['nodeB']

        lt_info = results['line_types'].get(line_type, {})
        wc_info = results['working_curves'].get(line_type, {})
        ic_info = results['syrope_ic'].get(line_id, {})
        owc_path = lt_info.get('owc_path')
        owc_path_resolved = None
        if owc_path:
            owc_path_resolved = os.path.normpath(os.path.join(os.path.dirname(file_path), owc_path))

        results['line_summary'][line_id] = {
            'owc_path': owc_path,
            'owc_path_resolved': owc_path_resolved,
            'a': lt_info.get('a'),
            'b': lt_info.get('b'),
            'C2': lt_info.get('C2'),
            'C1': lt_info.get('C1'),
            'WCType': wc_info.get('WCType'),
            'k1': wc_info.get('k1'),
            'k2': wc_info.get('k2'),
            'eps0': ((results['points'].get(node_a, {}).get('xyz', (0, 0, 0))[0] - results['points'].get(node_b, {}).get('xyz', (0, 0, 0))[0]) / line_info.get('unstrLen', 1)) - 1.0,
            'Tmax': ic_info.get('Tmax'),
            'Tmean': ic_info.get('Tmean'),
        }

    return results

def read_owc_dat(owc_file_path):
    """Read OWC strain-tension data from owc.dat.

    Expected format:
    - line 1: column names
    - line 2: units
    - remaining lines: <strain> <tension>
    """
    if not owc_file_path:
        raise ValueError("OWC file path is empty")
    if not os.path.exists(owc_file_path):
        raise FileNotFoundError(f"OWC file not found: {owc_file_path}")

    data = np.loadtxt(owc_file_path, skiprows=2)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError(f"Invalid OWC data format in: {owc_file_path}")

    return {
        'path': owc_file_path,
        'strain': data[:, 0],
        'tension': data[:, 1],
    }

def mean_strain_and_mean_tension_history(rope, eps0, Tmax0, dt, Tdur, n_phase=6):
    Tmax1 = Tmax0
    Tmax2 = np.interp(3.0 * eps0, rope.owc_strain, rope.owc_tension)

    rope.set_Tmax(Tmax1)
    wc_strain_1 = rope.wc_strain.copy()
    wc_tension_1 = rope.wc_tension.copy()

    rope.set_Tmax(Tmax2)
    wc_strain_2 = rope.wc_strain.copy()
    wc_tension_2 = rope.wc_tension.copy()

    t = np.linspace(0, Tdur * n_phase, int(Tdur * n_phase / dt) + 1)
    strain = np.zeros_like(t)
    tension = np.zeros_like(t)
    # First phase
    for i, tt in enumerate(t):
        if tt <= 1.0 * Tdur:
            strain[i] = eps0
            tension[i] = Tmean0
        elif tt <= 2.0 * Tdur:
            strain[i] = eps0 + (3.0 * eps0 - eps0) * (tt - 1.0 * Tdur) / Tdur
            tension[i] = np.interp(strain[i], wc_strain_1, wc_tension_1)
            if isclose(tension[i], Tmax1):
                tension[i] = np.interp(strain[i], rope.owc_strain, rope.owc_tension)
        elif tt <= 3.0 * Tdur:
            strain[i] = 3.0 * eps0
            tension[i] = Tmax2
        elif tt <= 4.0 * Tdur:
            strain[i] = 3.0 * eps0 - eps0 * (tt - 3.0 * Tdur) / Tdur
            tension[i] = np.interp(strain[i], wc_strain_2, wc_tension_2)
        elif tt <= 5.0 * Tdur:
            strain[i] = 2.0 * eps0
            tension[i] = np.interp(strain[i], wc_strain_2, wc_tension_2)
        else:
            strain[i] = 2.0 * eps0 + 1.5 * eps0 * (tt - 5.0 * Tdur) / Tdur
            tension[i] = np.interp(strain[i], wc_strain_2, wc_tension_2)
            if isclose(tension[i], Tmax2):
                tension[i] = np.interp(strain[i], rope.owc_strain, rope.owc_tension)     
    return t, strain, tension, wc_strain_1, wc_tension_1, wc_strain_2, wc_tension_2, Tmax1, Tmax2

def read_moordyn_output(file_path):
    data = np.loadtxt(file_path, skiprows=2)

    time    = data[:, 0]
    tension = data[:, 1]
    damping = data[:, 2]
    strain  = data[:, 8]

    return time, tension, damping, strain

def get_first_line_summary(parsed):
    if not parsed['line_summary']:
        raise ValueError("No parsed line summary found")
    first_line_id = sorted(parsed['line_summary'].keys())[0]
    return first_line_id, parsed['line_summary'][first_line_id]

def build_syrope_from_summary(summary):
    a = summary.get('a')
    b = summary.get('b')
    c1 = summary.get('C1')
    c2 = summary.get('C2')
    alpha = summary.get('k1')
    A = summary.get('k2')
    wc_type = summary.get('WCType')
    eps0 = summary.get('eps0')
    tmax = summary.get('Tmax')

    owc_data = read_owc_dat(summary.get('owc_path_resolved'))
    rope = Syrope(
        owc_data['strain'],
        owc_data['tension'],
        wc_type.lower(),
        alpha,
        A,
        a,
        b,
        c1,
        c2,
    )
    return rope, eps0, tmax

def compute_relative_l2_error(reference, numerical):
    return np.sqrt(np.sum((numerical - reference) ** 2) / np.sum(reference ** 2))

def describe_working_curve(summary):
    wc_type = summary.get('WCType')
    k1 = summary.get('k1')
    k2 = summary.get('k2')
    return f"{wc_type} (k1={k1}, k2={k2})"

def create_comparison_plot(
    time_num,
    strain_num,
    tension_num,
    damping_num,
    time_ref,
    strain_ref,
    mean_tension_ref,
    rope,
    wc_strain_1,
    wc_tension_1,
    wc_strain_2,
    wc_tension_2,
    eps0,
    tmax2,
    fast_loading=False,
    output_path='tension_strain.png',
):
    fig = plt.figure(figsize=(9, 5))
    gs = fig.add_gridspec(
        nrows=2,
        ncols=2,
        width_ratios=[0.75, 1],
        wspace=0.25,
        hspace=0.10,
    )

    ax_topleft = fig.add_subplot(gs[0, 0])
    ax_bottomleft = fig.add_subplot(gs[1, 0], sharex=ax_topleft)
    ax_right = fig.add_subplot(gs[:, 1])

    if fast_loading:
        ax_topleft.plot(time_num / 3600, strain_num / eps0, color='black', alpha=0.60, linewidth=0.80)
    else:
        ax_topleft.plot(time_num / 3600, strain_num / eps0, color='black')
    ax_topleft.plot(time_ref / 3600, strain_ref / eps0, '--', color='red', linewidth=1.0)
    ax_topleft.tick_params(labelbottom=False)
    ax_topleft.set_ylabel(r'$\epsilon/\epsilon_0$', fontsize=10)
    ax_topleft.set_xlim(0, time_ref[-1] / 3600)
    ax_topleft.set_ylim(0, 4)
    ax_topleft.set_xticks([0, 3, 6, 9, 12, 15, 18])
    for i in range(1, 6):
        ax_topleft.axvline(
            x=i * Tdur / 3600, color='k', linestyle='--', linewidth=0.5, alpha=0.5
        )
    for i in range(6):
        ax_topleft.text(
            1.5 + 3 * i,
            3.5,
            f"{phase_list[i]}",
            ha='center',
            va='center',
            fontsize=11,
            fontweight='bold',
        )

    if fast_loading:
        ax_bottomleft.plot(time_num / 3600, (tension_num - damping_num) / kMBL, color='black', alpha=0.60, linewidth=0.80)
        ax_bottomleft.plot(time_num / 3600, tension_num / kMBL, color='black')
    else:
        ax_bottomleft.plot(time_num / 3600, (tension_num - damping_num) / kMBL, color='black')
    ax_bottomleft.plot(time_ref / 3600, mean_tension_ref / kMBL, '--', color='red', linewidth=1.0)
    ax_bottomleft.set_xlabel(r'$t$ [h]', fontsize=10)
    ax_bottomleft.set_ylabel(r'$T/\mathrm{MBL}$', fontsize=10)
    ax_bottomleft.set_xlim(0, time_ref[-1] / 3600)
    ax_bottomleft.set_ylim(0, 0.40)
    ax_bottomleft.set_xticks([0, 3, 6, 9, 12, 15, 18])
    for i in range(1, 6):
        ax_bottomleft.axvline(
            x=i * Tdur / 3600, color='k', linestyle='--', linewidth=0.5, alpha=0.5
        )
    for i in range(6):
        ax_bottomleft.text(
            1.5 + 3 * i,
            0.35,
            f"{phase_list[i]}",
            ha='center',
            va='center',
            fontsize=11,
            fontweight='bold',
        )

    ax_right.plot(
        rope.owc_strain,
        rope.owc_tension / kMBL,
        color='g',
        label='Original working curve',
        alpha=1.00,
    )

    if fast_loading:
        ax_right.plot(
            strain_num,
            (tension_num - damping_num) / kMBL,
            '-',
            color='black',
            linewidth=0.80,
            alpha=0.60,
        )
    ax_right.plot(wc_strain_1, wc_tension_1 / kMBL, color='b', label='Working curves')
    ax_right.plot(wc_strain_2, wc_tension_2 / kMBL, color='b')

    if fast_loading:
        ax_right.plot(
            strain_ref,
            tension_num / kMBL,
            '-',
            color='black',
            label='Mean tension and strain',
        )
    else:
        ax_right.plot(
            strain_num,
            (tension_num - damping_num) / kMBL,
            '--',
            color='black',
            alpha=1.0,
            label='MoorDyn tension and strain',
        )

    ax_right.plot([0, 0.06], [0.16, 0.16], '--', linewidth=0.80, color='black', alpha=0.5)
    ax_right.plot(
        [0, 0.06],
        [tmax2 / kMBL, tmax2 / kMBL],
        '--',
        linewidth=0.80,
        color='black',
        alpha=0.5,
    )
    ax_right.text(0.002, 0.17, r'$0.16$', fontsize=10, color='black', alpha=0.8)
    ax_right.text(
        0.002,
        tmax2 / kMBL + 0.01,
        f"{tmax2 / kMBL:.2f}",
        fontsize=10,
        color='black',
        alpha=0.8,
    )

    ax_right.set_xlabel(r'$\epsilon$', fontsize=10)
    ax_right.set_ylabel(r'$T/\mathrm{MBL}$', fontsize=10)
    ax_right.set_xlim(0, 0.06)
    ax_right.set_ylim(0, 0.50)

    leg = ax_right.legend(
        fontsize=10,
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        edgecolor='black',
        facecolor='white',
    )
    leg.get_frame().set_linewidth(0.60)
    leg.get_frame().set_linestyle('-')

    plt.savefig(output_path, dpi=300, bbox_inches='tight')

def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python mean_tension_verification.py <all|linear|exp|quadratic> <True|False>")

    wc_mode = sys.argv[1].strip().lower()
    has_fast_loading_arg = sys.argv[2].strip()
    if has_fast_loading_arg not in {'True', 'False'}:
        raise SystemExit("Invalid has_fast_loading value. Use True or False")
    has_fast_loading = has_fast_loading_arg == 'True'

    wc_dir_map = {
        'linear': 'linear_wc',
        'exp': 'exp_wc',
        'quadratic': 'quadratic_wc',
    }
    if wc_mode == 'all':
        modes = ['linear', 'exp', 'quadratic']
    elif wc_mode in wc_dir_map:
        modes = [wc_mode]
    else:
        raise SystemExit(f"Invalid mode '{wc_mode}'. Use one of: all, linear, exp, quadratic")

    summary_rows = []

    for mode in modes:
        file_path = os.path.join(
            os.path.dirname(__file__),
            '..',
            'input',
            wc_dir_map[mode],
            'line.txt',
        )

        moordyn_output = os.path.join(os.path.dirname(file_path), 'line_Line1.out')
        time_num, tension_num, damping_num, strain_num = read_moordyn_output(moordyn_output)

        parsed = parse_syrope_inputs(file_path)

        print(f"Input ({mode}): {file_path}")
        try:
            first_line_id, summary = get_first_line_summary(parsed)
        except ValueError:
            print("No lines parsed.")
            continue
        print(f"Line {first_line_id} summary: {summary}")

        rope, eps0, tmax = build_syrope_from_summary(summary)

        # Check if Tmax from input file is consistent with target Tmax
        assert isclose(abs(tmax - Tmax0) / Tmax0, 0.0, abs_tol=1e-3), f"Input Tmax {tmax} is not close to target Tmax {Tmax0}"

        time, strain_ref, mean_tension_ref, wc_strain_1, wc_tension_1, wc_strain_2, wc_tension_2, Tmax1, Tmax2 = mean_strain_and_mean_tension_history(rope, eps0, tmax, 0.1, Tdur)

        # L2 error between numerical and reference mean tension
        error = compute_relative_l2_error(mean_tension_ref, tension_num)
        print(f"L2-error between tension and tension_ref: {error:.4e}")

        summary_rows.append({
            'mode': mode,
            'working_curve': describe_working_curve(summary),
            'fast_loading': has_fast_loading,
            'fast_loading_where': 'Right panel strain-tension traces in create_comparison_plot()',
            'l2_error': error,
        })

        create_comparison_plot(
            time_num=time_num,
            strain_num=strain_num,
            tension_num=tension_num,
            damping_num=damping_num,
            time_ref=time,
            strain_ref=strain_ref,
            mean_tension_ref=mean_tension_ref,
            rope=rope,
            wc_strain_1=wc_strain_1,
            wc_tension_1=wc_tension_1,
            wc_strain_2=wc_strain_2,
            wc_tension_2=wc_tension_2,
            eps0=eps0,
            tmax2=Tmax2,
            fast_loading=has_fast_loading,
            output_path=os.path.join(os.path.dirname(file_path), 'tension_strain.png'),
        )

    print("\n=== Summary ===")
    if not summary_rows:
        print("No cases were processed.")
        return

    for row in summary_rows:
        print(f"Case: {row['mode']}")
        print(f"  Working curve formula: {row['working_curve']}")
        print(f"  Fast loading applied: {row['fast_loading']}")
        print(f"  Where applied: {row['fast_loading_where']}")
        print(f"  L2 error: {row['l2_error']:.4e}")

if __name__ == "__main__":
    main()