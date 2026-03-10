import numpy as np

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
                "wc_slow_spring_strain_static must be strictly increasing"
            )

    def find_strains(self, Tmean):
        if Tmean < self.Tmax:  # On the working curve
            strain_static = np.interp(Tmean, self.wc_tension, self.wc_strain)
            fast_spring_strain_static = (
                1.0 / self.beta * np.log(1.0 + self.beta / self.alpha * Tmean)
            )
            slow_spring_strain_static = strain_static - fast_spring_strain_static
        else:  # On the original working curve
            strain_static = np.interp(Tmean, self.owc_tension, self.owc_strain)
            fast_spring_strain_static = (
                1.0 / self.beta * np.log(1.0 + self.beta / self.alpha * Tmean)
            )
            slow_spring_strain_static = strain_static - fast_spring_strain_static

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
    
    def set_Tmean(self, Tmean):
        self.Tmean = Tmean
        self.strain_static, self.fast_spring_strain_static, self.slow_spring_strain_static = self.find_strains(Tmean)

    def slow_spring_strain_rate_from_instantaneous_strain(self, total_eps, dtotal_eps, slow_eps):
        if self.Tmean < self.Tmax:
            # On the working curve
            Tmean = np.interp(slow_eps, self.wc_slow_spring_strain_static, self.wc_tension)
            K1 = self.alpha + self.beta * Tmean
            dslow_eps = (K1 * (total_eps - np.interp(Tmean, self.wc_tension, self.wc_strain)) + self.c1 * dtotal_eps) / (self.c1 + self.c2)
            tension = Tmean + K1 * (total_eps - np.interp(Tmean, self.wc_tension, self.wc_strain))
        else:
            # On the original working curve
            Tmean = np.interp(slow_eps, self.owc_slow_spring_strain_static, self.owc_tension)
            K1 = self.alpha + self.beta * Tmean
            dslow_eps = (K1 * (total_eps - np.interp(Tmean, self.owc_tension, self.owc_strain)) + self.c1 * dtotal_eps) / (self.c1 + self.c2)
            tension = Tmean + K1 * (total_eps - np.interp(Tmean, self.owc_tension, self.owc_strain))
        
        if (Tmean > self.Tmax):
            self.Tmax = Tmean
            self.set_Tmax(self.Tmax)
        
        return dslow_eps, tension

    def slow_spring_strain_rate_from_instantaneous_tension(self, tension, slow_eps):
        if self.Tmean < self.Tmax:
            # On the working curve
            self.Tmean = np.interp(slow_eps, self.wc_slow_spring_strain_static, self.wc_tension)
            K1 = self.alpha + self.beta * self.Tmean
            dslow_eps = (tension - self.Tmean) / self.c2
            total_eps = (tension - self.Tmean) / K1 + np.interp(self.Tmean, self.wc_tension, self.wc_strain)
        else:
            # On the original working curve
            self.Tmean = np.interp(slow_eps, self.owc_slow_spring_strain_static, self.owc_tension)
            K1 = self.alpha + self.beta * self.Tmean
            dslow_eps = (tension - self.Tmean) / self.c2
            total_eps = (tension - self.Tmean) / K1 + np.interp(self.Tmean, self.owc_tension, self.owc_strain)
        
        if (self.Tmean > self.Tmax):
            self.Tmax = self.Tmean
            self.set_Tmax(self.Tmax)
        
        return dslow_eps, total_eps