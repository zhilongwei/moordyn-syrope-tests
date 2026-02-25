import numpy as np
from pathlib import Path

from syropepy import Syrope

kG = 9.80665
kMBL = 885.0e3 * kG

Tmax0 = 0.16 * kMBL
Tmean0 = 0.05 * kMBL

print(f"Using MBL = {kMBL:.4e} N, Tmax0 = {Tmax0:.4e} N, Tmean0 = {Tmean0:.4e} N")

alpha = 17.667 * kMBL
beta = 0.2313 * 100

c1 = 1.0e5

CASES = [
	("Linear", "linear", 1.25e8, 0.0, 5.0e10),
	("Quadratic", "quadratic", 0.25, 1.0, 2.5e10),
	("Exponential", "exp", 0.2, 1.1, 2.5e10),
]


def evaluate_case(owc, case_name, wc_type, p1, p2, c2):
	syrope = Syrope(owc[:, 0], owc[:, 1], wc_type, p1, p2, alpha, beta, c1, c2)
	syrope.set_Tmax(Tmax0)

	eps, eps_1, eps_2 = syrope.find_strains(Tmean0)
	k1, k2, omega_critical = syrope.critical_frequency(Tmean0)
	period = (2.0 * np.pi / omega_critical) if omega_critical > 0.0 else np.inf

	return {
		"case": case_name,
		"eps": eps,
		"eps_1": eps_1,
		"eps_2": eps_2,
		"k1_mbl": k1 / kMBL,
		"k2_mbl": k2 / kMBL,
		"omega_critical": omega_critical,
		"period": period,
	}


def print_results(results):
	print("\nSYROPE INITIAL CONDITIONS SUMMARY")
	print(f"Tmean = {Tmean0 / kMBL:.3f} * MBL, Tmax = {Tmax0 / kMBL:.3f} * MBL")
	print()
	print(
		f"{'WC':<12}"
		f"{'eps':>12}"
		f"{'eps_1':>12}"
		f"{'eps_2':>12}"
		f"{'k1/MBL':>14}"
		f"{'k2/MBL':>14}"
		f"{'omega_c [rad/s]':>18}"
		f"{'period [s]':>12}"
	)
	print("-" * 106)

	for r in results:
		print(
			f"{r['case']:<12}"
			f"{r['eps']:>12.5e}"
			f"{r['eps_1']:>12.5e}"
			f"{r['eps_2']:>12.5e}"
			f"{r['k1_mbl']:>14.5e}"
			f"{r['k2_mbl']:>14.5e}"
			f"{r['omega_critical']:>18.5e}"
			f"{r['period']:>12.3f}"
		)


def main():
	script_dir = Path(__file__).parent
	owc = np.loadtxt((script_dir / "../input/owc.dat").resolve(), skiprows=2)

	results = [evaluate_case(owc, *case) for case in CASES]
	print_results(results)


if __name__ == "__main__":
	main()