"use strict";
const PRECISION = 5;
const TEST_ISENTROPIC = false;
const TEST_NORMAL_SHOCK = true;
const TEST_OBLIQUE_SHOCK = false;
const deflection_tolerance = 0.001;
const epsilon = 0.00001;
const epsilon_derivative = 0.001;
const step = 0.1;
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
    // static A_by_Astar_pressure(M:number, gamma:number, p:number, p0:number):number{
    //   const gp = gamma+1;
    //   const gn = gamma-1;
    //   const pr = p/p0;
    //   const num = Math.pow(2/gp, 0.5 * gp/gn) * Math.pow(0.5 * gn, 0.5);
    //   const den = Math.pow(Math.pow(pr, 2/gamma) - Math.pow(pr, (gp)/gamma), 0.5);
    //   return num/den;
    // }
    static A_by_Astar_mach(M, gamma) {
        const pow = (gamma + 1) / (2 * (gamma - 1));
        const t1 = (2 / (gamma + 1));
        const t2 = (1 + 0.5 * (gamma - 1) * M * M);
        const i1 = Math.pow(t1 * t2, pow);
        return i1 / M;
    }
    // Now for the inverse functions
    static findMFrom_Tt_by_T(ratio, gamma) {
        const m1sq = (2 / (gamma - 1)) * (ratio - 1);
        return Math.sqrt(m1sq);
    }
    static findMFrom_Pt_by_P(ratio, gamma) {
        const pow = (gamma - 1) / gamma;
        const t_ratio = Math.pow(ratio, pow);
        return Isentropic.findMFrom_Tt_by_T(t_ratio, gamma);
    }
    static findMFrom_rhot_by_rho(ratio, gamma) {
        const pow = gamma - 1;
        const t_ratio = Math.pow(ratio, pow);
        return Isentropic.findMFrom_Tt_by_T(t_ratio, gamma);
    }
    static findMFrom_A_by_Astar_mach(ratio, gamma) {
        // Numerical solution
        let M = 0.9999; // Initial guess
        for (let i = 0; i < 10000; i++) {
            const f1 = Isentropic.A_by_Astar_mach(M, gamma);
            const f2 = Isentropic.A_by_Astar_mach(M - epsilon_derivative, gamma);
            const deriv_A_by_Astar = (f2 - f1) / epsilon_derivative;
            const M_new = M - step * (ratio - f1) / deriv_A_by_Astar;
            const new_ratio = Isentropic.A_by_Astar_mach(M_new, gamma);
            if (Math.abs(new_ratio - ratio) < 0.00001) {
                return M_new;
            }
            M = M_new;
        }
        throw new Error("Could not find Mach number from A/A*");
    }
}
class NormalShock {
    static downStreamMachNumber(M1, gamma) {
        const t1 = (gamma - 1);
        const m1sq = M1 * M1;
        const m2sq = (t1 * m1sq + 2) / (2 * gamma * m1sq - t1);
        return Math.sqrt(m2sq);
    }
    static upStreamMachNumber(M2, gamma) {
        const num = (gamma - 1) * M2 * M2 + 2;
        const den = 2 * gamma * M2 * M2 - (gamma - 1);
        return Math.sqrt(num / den);
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
        return (1 / (NormalShock.P2_by_P1(M1, gamma)) / Isentropic.Pt_by_P(NormalShock.downStreamMachNumber(M1, gamma), gamma));
    }
    // to allow for user flexibility, we find M1 by various other things. Inverse functions, basically
    // Need to manually calculate these
    static findM1From_P2_by_P1(ratio, gamma) {
        const num = (gamma - 1) + ratio * (gamma + 1);
        const den = 2 * gamma;
        return Math.sqrt(num / den);
    }
    static findM1From_RHO2_by_RHO1(ratio, gamma) {
        const num = 2 * ratio;
        const den = (gamma + 1) - ratio * (gamma - 1);
        return Math.sqrt(num / den);
    }
    static findM1From_T2_by_T1(ratio, gamma) {
        const a = 2 * gamma * (gamma - 1);
        const b = 4 * gamma - (gamma - 1) * (gamma - 1) - (gamma + 1) * (gamma + 1) * ratio;
        const c = -2 * (gamma - 1);
        const msq = (-b + Math.sqrt(b * b - 4 * a * c)) / (2 * a);
        return Math.sqrt(msq);
    }
}
if (TEST_NORMAL_SHOCK) {
    console.log("\nMach2 Normal Shock test");
    // console.log("M2 = ", NormalShock.downStreamMachNumber(2, 1.4).toFixed(PRECISION));
    // console.log("P2/P1 = ", NormalShock.P2_by_P1(2, 1.4).toFixed(PRECISION));
    // console.log("rho2/rho1 = ", NormalShock.RHO2_by_RHO1(2, 1.4).toFixed(PRECISION));
    // console.log("T2/T1 = ", NormalShock.T2_by_T1(2, 1.4).toFixed(PRECISION));
    // console.log("a2/a1 = ", NormalShock.a2_by_a1(2, 1.4).toFixed(PRECISION));
    // console.log("Pt2/Pt1 = ", NormalShock.Pt2_by_Pt1(2, 1.4).toFixed(PRECISION));
    // console.log("P1/Pt2 = ", NormalShock.P1_by_Pt2(2, 1.4).toFixed(PRECISION));
    // console.log("Mach2 Normal Shock test complete\n")
    const m1 = 2.0;
    const gamma = 1.4;
    const p2by1 = NormalShock.P2_by_P1(m1, gamma);
    const rho2by1 = NormalShock.RHO2_by_RHO1(m1, gamma);
    const t2by1 = NormalShock.T2_by_T1(m1, gamma);
    const m2 = NormalShock.downStreamMachNumber(m1, gamma);
    console.log("Given M1 = ", m1);
    console.log("P2/P1 = ", p2by1.toFixed(PRECISION));
    console.log("rho2/rho1 = ", rho2by1.toFixed(PRECISION));
    console.log("T2/T1 = ", t2by1.toFixed(PRECISION));
    console.log("M2 = ", m2.toFixed(PRECISION));
    console.log("Testing inverse functions...");
    console.log("M1 from P2/P1 = ", NormalShock.findM1From_P2_by_P1(p2by1, gamma).toFixed(PRECISION));
    console.log("M1 from rho2/rho1 = ", NormalShock.findM1From_RHO2_by_RHO1(rho2by1, gamma).toFixed(PRECISION));
    console.log("M1 from T2/T1 = ", NormalShock.findM1From_T2_by_T1(t2by1, gamma).toFixed(PRECISION));
    console.log("Mach2 Normal Shock test complete\n");
}
class ObliqueShock {
    static deflectionAngle(M1, beta, gamma) {
        const num = 2 * (Math.cos(beta) / Math.sin(beta)) * (Math.pow(M1 * Math.sin(beta), 2) - 1);
        const den = 2 + M1 * M1 * (gamma + Math.cos(2 * beta));
        return Math.atan(num / den);
    }
    static downStreamMachNumber(M1, beta, gamma) {
        const deflection = ObliqueShock.deflectionAngle(M1, beta, gamma);
        const t0 = Math.pow(Math.sin(beta - deflection), 2);
        const num = Math.pow(M1 * Math.sin(beta), 2) + (2 / (gamma - 1));
        const den = (2 * gamma / (gamma - 1)) * Math.pow(M1 * Math.sin(beta), 2) - 1;
        return Math.sqrt(num / (den * t0));
    }
    static P2_by_P1(M1, beta, gamma) {
        const num = 2 * gamma * Math.pow(M1 * Math.sin(beta), 2) - (gamma - 1);
        const den = gamma + 1;
        return num / den;
    }
    static RHO2_by_RHO1(M1, beta, gamma) {
        const num = (gamma + 1) * Math.pow(M1 * Math.sin(beta), 2);
        const den = (gamma - 1) * Math.pow(M1 * Math.sin(beta), 2) + 2;
        return num / den;
    }
    static T2_by_T1(M1, beta, gamma) {
        const p2byb1 = ObliqueShock.P2_by_P1(M1, beta, gamma);
        const rho2byrho1 = ObliqueShock.RHO2_by_RHO1(M1, beta, gamma);
        return p2byb1 / rho2byrho1;
    }
    static find_max_deflection_angle(M1, gamma) {
        const t1 = (gamma + 1) / (4 * gamma);
        const t2_0 = -1 / (gamma * M1 * M1);
        const t2_11 = gamma + 1;
        const t2_12 = 1 + 0.5 * (gamma - 1) * M1 * M1 + (gamma + 1) * Math.pow(M1, 4) / 16;
        const t2_1 = Math.sqrt(t2_11 * t2_12);
        const t = t1 + t2_0 * (1 - t2_1);
        return Math.asin(Math.sqrt(t));
    }
    static find_shock_angle_solutions(M1, gamma, deflection) {
        // const lhs = Math.tan(deflection);
        // const max_deflection = ObliqueShock.find_max_deflection_angle(M1, gamma);
        // if( )
        let P1 = Math.PI / 2 - epsilon;
        let P2 = Math.asin(1 / M1) + epsilon;
        let P1_final = P1;
        let P2_final = P2;
        for (let i = 0; i < 10000; i++) {
            const f1 = ObliqueShock.deflectionAngle(M1, P1, gamma);
            const f2 = ObliqueShock.deflectionAngle(M1, P1 - epsilon, gamma);
            const deriv_deflection = (f2 - f1) / epsilon;
            const P1_new = P1 - step * (deflection - f1) / deriv_deflection;
            const new_deflection = ObliqueShock.deflectionAngle(M1, P1_new, gamma);
            if (Math.abs(new_deflection - deflection) < deflection_tolerance) {
                P1_final = P1;
                break;
            }
            P1 = P1_new;
        }
        for (let i = 0; i < 10000; i++) {
            const f1 = ObliqueShock.deflectionAngle(M1, P2, gamma);
            const f2 = ObliqueShock.deflectionAngle(M1, P2 + epsilon, gamma);
            const deriv_deflection = (f2 - f1) / epsilon;
            const P2_new = P2 + step * (deflection - f1) / deriv_deflection;
            const new_deflection = ObliqueShock.deflectionAngle(M1, P2_new, gamma);
            if (Math.abs(new_deflection - deflection) < deflection_tolerance) {
                P2_final = P2;
                break;
            }
            P2 = P2_new;
        }
        if (Math.abs(P1_final - P2_final) < deflection_tolerance) {
            const average = (P1_final + P2_final) / 2;
            P1_final = average;
            P2_final = average;
        }
        return [P1_final, P2_final];
    }
    static findStrongWeakSolutions(M1, gamma, deflection) {
        const max_deflection = ObliqueShock.find_max_deflection_angle(M1, gamma);
        if (deflection > max_deflection) {
            throw new Error("Deflection angle exceeds maximum deflection angle for given Mach number");
        }
        const betas = ObliqueShock.find_shock_angle_solutions(M1, gamma, deflection);
        const M2s = [ObliqueShock.downStreamMachNumber(M1, betas[1], gamma), ObliqueShock.downStreamMachNumber(M1, betas[0], gamma)];
        const P2_by_P1s = [ObliqueShock.P2_by_P1(M1, betas[1], gamma), ObliqueShock.P2_by_P1(M1, betas[0], gamma)];
        const rho2_by_rho1s = [ObliqueShock.RHO2_by_RHO1(M1, betas[1], gamma), ObliqueShock.RHO2_by_RHO1(M1, betas[0], gamma)];
        const T2_by_T1s = [ObliqueShock.T2_by_T1(M1, betas[1], gamma), ObliqueShock.T2_by_T1(M1, betas[0], gamma)];
        return [betas, M2s, P2_by_P1s, rho2_by_rho1s, T2_by_T1s];
    }
    static normalMachs(M1, gamma, deflection) {
        const betas = ObliqueShock.find_shock_angle_solutions(M1, gamma, deflection);
        const M1n = [M1 * Math.sin(betas[0]), M1 * Math.sin(betas[1])];
        const M2n = [ObliqueShock.downStreamMachNumber(M1, betas[0], gamma) * Math.sin(betas[0] - deflection), ObliqueShock.downStreamMachNumber(M1, betas[1], gamma) * Math.sin(betas[1] - deflection)];
        return [M1n, M2n];
    }
}
if (TEST_OBLIQUE_SHOCK) {
    console.log("\nOblique Shock test");
    const M1 = 2.0;
    const gamma = 1.4;
    const deflection = 20 * Math.PI / 180;
    const [beta_weak, beta_strong] = ObliqueShock.find_shock_angle_solutions(M1, gamma, deflection);
    console.log("M1 = ", M1);
    console.log("Deflection (deg) = ", (deflection * 180 / Math.PI).toFixed(PRECISION));
    console.log("Weak Shock Angle (deg) = ", (beta_weak * 180 / Math.PI).toFixed(PRECISION));
    console.log("Strong Shock Angle (deg) = ", (beta_strong * 180 / Math.PI).toFixed(PRECISION));
}
