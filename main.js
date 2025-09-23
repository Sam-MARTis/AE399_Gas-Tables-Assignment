"use strict";
class Isentropic {
    static Tt_by_T(M, gamma) {
        return 1 + ((gamma - 1) * 0.5) * M * M;
    }
    static Pt_by_P(M, gamma) {
        const ttr = Isentropic.Tt_by_T(M, gamma);
        return Math.pow(ttr, gamma / (gamma - 1));
    }
    static rhot_by_rho(M, gamma) {
        const ttr = Isentropic.Tt_by_T(M, gamma);
        return Math.pow(ttr, 1 / (gamma - 1));
    }
    static at_by_a(M, gamma) {
        const ttr = Isentropic.Tt_by_T(M, gamma);
        return Math.sqrt(ttr);
    }
    tstar_by_t(M, gamma) {
        return 2 / (gamma + 1);
    }
    pstar_by_p(M, gamma) {
        const tstar = this.tstar_by_t(M, gamma);
        return Math.pow(tstar, gamma / (gamma - 1));
    }
    rhostar_by_rho(M, gamma) {
        const tstar = this.tstar_by_t(M, gamma);
        return Math.pow(tstar, 1 / (gamma - 1));
    }
    a_by_astar_pressure(M, gamma, p, p0) {
        const gp = gamma + 1;
        const gn = gamma - 1;
        const pr = p / p0;
        const num = Math.pow(2 / gp, 0.5 * gp / gn) * Math.pow(0.5 * gn, 0.5);
        const den = Math.pow(Math.pow(pr, 2 / gamma) - Math.pow(pr, (gp) / gamma), 0.5);
        return num / den;
    }
    a_by_astar_Mach(M, gamma) {
        const pow = (gamma + 1) / (2 * (gamma - 1));
        const t1 = (2 / (gamma + 1));
        const t2 = (1 + 0.5 * (gamma - 1) * M * M);
        const i1 = Math.pow(t1 * t2, pow);
        return i1 / M;
    }
}
class NormalShock {
    static DownStreamMachNumber(M1, gamma) {
        const t1 = (gamma - 1);
        const m1sq = M1 * M1;
        const m2sq = (t1 * m1sq + 2) / (2 * gamma * m1sq - t1);
        return Math.sqrt(m2sq);
    }
    static P2_by_P1(M1, gamma) {
        const m1sq = M1 * M1;
        return (2 * gamma * m1sq - (gamma - 1)) / (gamma + 1);
    }
    static RHO2_by_RHO1(M1, gamma) {
        const m1sq = M1 * M1;
        return ((gamma + 1) * m1sq) / ((gamma - 1) * m1sq + 2);
    }
    static T2_by_T1(M1, gamma) {
        const p2byb1 = NormalShock.P2_by_P1(M1, gamma);
        const rho2byrho1 = NormalShock.RHO2_by_RHO1(M1, gamma);
        return p2byb1 / rho2byrho1;
    }
    static A2_by_A1(M1, gamma) {
        return Math.sqrt(NormalShock.T2_by_T1(M1, gamma));
    }
}
