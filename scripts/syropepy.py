import numpy as np
from cmath import isclose


class Syrope:
    def __init__(
        self,
        owc_strain,
        owc_tension,
        wc_mod,
        p1,
        p2,
        alpha,
        beta,
        c1,
        c2,
        N=30,
    ):
        assert wc_mod in ["linear", "quadratic", "exp"], (
            "wc_mod must be 'linear', 'quadratic', or 'exp'"
        )

        self.owc_strain = owc_strain
        self.owc_tension = owc_tension
        self.wc_mod = wc_mod
        self.p1 = p1
        self.p2 = p2
        self.alpha = alpha
        self.beta = beta
        self.c1 = c1
        self.c2 = c2
        self.N = N

        self.owc_fast_spring_strain_static = (
            1.0 / self.beta * np.log(1.0 + self.beta / self.alpha * self.owc_tension)
        )
        self.owc_slow_spring_strain_static = (
            self.owc_strain - self.owc_fast_spring_strain_static
        )

        # Check if owc_slow_spring_strain_static is strictly increasing
        if not np.all(np.diff(self.owc_slow_spring_strain_static) > 0):
            raise ValueError(
                "owc_slow_spring_strain_static must be strictly increasing"
            )

    def set_Tmax(self, Tmax):
        self.Tmax = Tmax

        eps_max = np.interp(self.Tmax, self.owc_tension, self.owc_strain)
        eps_min = self.owc_strain[0] + self.p1 * (eps_max - self.owc_strain[0])

        if self.wc_mod == "linear" and self.p1 >= 1.0:
            eps_min = eps_max - self.Tmax / self.p1

        self.wc_strain = np.linspace(eps_min, eps_max, self.N)
        xi = (self.wc_strain - eps_min) / (eps_max - eps_min)  # Normalized strain

        if self.wc_mod == "linear":
            self.wc_tension = self.Tmax * xi
        elif self.wc_mod == "quadratic":
            self.wc_tension = self.Tmax * xi * (self.p2 * xi + (1.0 - self.p2))
        elif self.wc_mod == "exp":
            self.wc_tension = (
                self.Tmax * (1.0 - np.exp(self.p2 * xi)) / (1.0 - np.exp(self.p2))
            )

        self.wc_fast_spring_strain_static = (
            1.0 / self.beta * np.log(1.0 + self.beta / self.alpha * self.wc_tension)
        )
        self.wc_slow_spring_strain_static = (
            self.wc_strain - self.wc_fast_spring_strain_static
        )

        # Check if wc_slow_spring_strain_static is strictly increasing
        if not np.all(np.diff(self.wc_slow_spring_strain_static) > 0):
            raise ValueError(
                "owc_slow_spring_strain_static must be strictly increasing"
            )

    def find_strains(self, Tmean):
        if Tmean < self.Tmax:  # On the working curve
            strain_static = np.interp(Tmean, self.wc_tension, self.wc_strain)
            fast_spring_strain_static = (
                1.0 / self.beta * np.log(1.0 + self.beta / self.alpha * Tmean)
            )
            slow_spring_strain_static = strain_static - fast_spring_strain_static
        else:  # On the overworking curve
            strain_static = np.interp(Tmean, self.owc_tension, self.owc_strain)
            fast_spring_strain_static = (
                1.0 / self.beta * np.log(1.0 + self.beta / self.alpha * Tmean)
            )
            slow_spring_strain_static = strain_static - fast_spring_strain_static

        return strain_static, fast_spring_strain_static, slow_spring_strain_static

    def find_tension(self, strain, slow_spring_strain_static):
        # Check if the mean tension is below Tmax
        Tmean = np.interp(
            slow_spring_strain_static,
            self.wc_slow_spring_strain_static,
            self.wc_tension,
        )

        if Tmean < self.Tmax:
            strain_static = np.interp(Tmean, self.wc_tension, self.wc_strain)
            fast_spring_strain_static = (
                1.0 / self.beta * np.log(1.0 + self.beta / self.alpha * Tmean)
            )
            slow_spring_strain_static = strain_static - fast_spring_strain_static
        else:  # On the original working curve
            Tmean = np.interp(
                slow_spring_strain_static,
                self.owc_slow_spring_strain_static,
                self.owc_tension,
            )
            strain_static = np.interp(Tmean, self.owc_tension, self.owc_strain)
            fast_spring_strain_static = (
                1.0 / self.beta * np.log(1.0 + self.beta / self.alpha * Tmean)
            )
            slow_spring_strain_static = strain_static - fast_spring_strain_static

        if isclose(strain, strain_static):
            print("The provided strain corresponds to the calculated tension.")
        else:
            print("Warning: Tension may not correspond to the given strain.")

        return strain_static, fast_spring_strain_static, slow_spring_strain_static

    def critical_frequency(self, Tmean, delta=0.01):
        EA_d = self.alpha + self.beta * Tmean
        EA_1 = EA_d

        dTmean = delta * Tmean
        if Tmean < self.Tmax:
            eps_slow_spring_1 = np.interp(
                Tmean - dTmean, self.wc_tension, self.wc_slow_spring_strain_static
            )
            eps_slow_spring_2 = np.interp(
                Tmean + dTmean, self.wc_tension, self.wc_slow_spring_strain_static
            )
        else:
            eps_slow_spring_1 = np.interp(
                Tmean - dTmean, self.owc_tension, self.owc_slow_spring_strain_static
            )
            eps_slow_spring_2 = np.interp(
                Tmean + dTmean, self.owc_tension, self.owc_slow_spring_strain_static
            )

        EA_2 = (2 * dTmean) / (
            eps_slow_spring_2 - eps_slow_spring_1
        )  # Local slow spring stiffness

        return  EA_1, EA_2, (
            (EA_1 + EA_2)
            / self.c2
            * np.sqrt((EA_1 + 4.0 * EA_2) / (3.0 * EA_1 + 4.0 * EA_2))
        )
