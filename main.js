"use strict";
const PRECISION = 5;
const TEST_ISENTROPIC = false;
const TEST_NORMAL_SHOCK = true;
const TEST_OBLIQUE_SHOCK = false;
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
    static tstar_by_t(M, gamma) {
        return 2 / (gamma + 1);
    }
    static pstar_by_p(M, gamma) {
        const tstar = this.tstar_by_t(M, gamma);
        return Math.pow(tstar, gamma / (gamma - 1));
    }
    static rhostar_by_rho(M, gamma) {
        const tstar = this.tstar_by_t(M, gamma);
        return Math.pow(tstar, 1 / (gamma - 1));
    }
    static A_by_Astar_pressure(M, gamma, p, p0) {
        const gp = gamma + 1;
        const gn = gamma - 1;
        const pr = p / p0;
        const num = Math.pow(2 / gp, 0.5 * gp / gn) * Math.pow(0.5 * gn, 0.5);
        const den = Math.pow(Math.pow(pr, 2 / gamma) - Math.pow(pr, (gp) / gamma), 0.5);
        return num / den;
    }
    static A_by_Astar_mach(M, gamma) {
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
    static a2_by_a1(M1, gamma) {
        return Math.sqrt(NormalShock.T2_by_T1(M1, gamma));
    }
    static Pt2_by_Pt1(M1, gamma) {
        const t1 = (gamma + 1) * M1 * M1 * 0.5 / (1 + ((gamma - 1) * 0.5) * M1 * M1);
        const pow1 = gamma / (gamma - 1);
        const t2_1 = 2 * gamma * M1 * M1 / (gamma + 1);
        const t2_2 = -(gamma - 1) / (gamma + 1);
        const t2 = t2_1 + t2_2;
        const pow2 = -1 / (gamma - 1);
        const t1f = Math.pow(t1, pow1);
        const t2f = Math.pow(t2, pow2);
        return t1f * t2f;
    }
    static P1_by_Pt2(M1, gamma) {
        return (1 / (NormalShock.P2_by_P1(M1, gamma)) / Isentropic.Pt_by_P(NormalShock.DownStreamMachNumber(M1, gamma), gamma));
    }
}
if (TEST_ISENTROPIC) {
    console.log("\nMach2 Normal Shock test");
    console.log("M2 = ", NormalShock.DownStreamMachNumber(2, 1.4).toFixed(PRECISION));
    console.log("P2/P1 = ", NormalShock.P2_by_P1(2, 1.4).toFixed(PRECISION));
    console.log("rho2/rho1 = ", NormalShock.RHO2_by_RHO1(2, 1.4).toFixed(PRECISION));
    console.log("T2/T1 = ", NormalShock.T2_by_T1(2, 1.4).toFixed(PRECISION));
    console.log("a2/a1 = ", NormalShock.a2_by_a1(2, 1.4).toFixed(PRECISION));
    console.log("Pt2/Pt1 = ", NormalShock.Pt2_by_Pt1(2, 1.4).toFixed(PRECISION));
    console.log("P1/Pt2 = ", NormalShock.P1_by_Pt2(2, 1.4).toFixed(PRECISION));
    console.log("Mach2 Normal Shock test complete\n");
}
