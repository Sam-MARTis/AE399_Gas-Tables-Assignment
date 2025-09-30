let PRECISION = 10;
let γ = 1.4;
let tolerance = 0.0001;
const TEST_ISENTROPIC = true;
const TEST_NORMAL_SHOCK = false;
const TEST_OBLIQUE_SHOCK = false;
const epsilon = 0.00001;
const epsilon_derivative = 0.001;
const step = 0.1;
const max_M_step = 0.04;
const MAX_ITER = 100000;

class Isentropic {
  static getMachAngle(M: number): number {
    if (M >= 1) {
      return Math.asin(1 / M);
    }
    return NaN;
  }

  static Tt_by_T(M: number, gamma: number): number {
    return 1 + (gamma - 1) * 0.5 * M * M;
  }
  static Pt_by_P(M: number, gamma: number): number {
    const ttr = Isentropic.Tt_by_T(M, gamma);
    return Math.pow(ttr, gamma / (gamma - 1));
  }
  static rhot_by_rho(M: number, gamma: number): number {
    const ttr = Isentropic.Tt_by_T(M, gamma);
    return Math.pow(ttr, 1 / (gamma - 1));
  }
  static at_by_a(M: number, gamma: number): number {
    const ttr = Isentropic.Tt_by_T(M, gamma);
    return Math.sqrt(ttr);
  }
  static tstar_by_t(M: number, gamma: number): number {
    return 2 / (gamma + 1);
  }
  static pstar_by_p(M: number, gamma: number): number {
    const tstar = this.tstar_by_t(M, gamma);
    return Math.pow(tstar, gamma / (gamma - 1));
  }
  static rhostar_by_rho(M: number, gamma: number): number {
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
  static A_by_Astar_mach(M: number, gamma: number): number {
    const pow = (gamma + 1) / (2 * (gamma - 1));
    const t1 = 2 / (gamma + 1);
    const t2 = 1 + 0.5 * (gamma - 1) * M * M;
    const i1 = Math.pow(t1 * t2, pow);
    return i1 / M;
  }
  static prandtlMeyerAngle(M: number, gamma: number): number {
    const t1 = Math.sqrt((gamma + 1) / (gamma - 1));
    const t2 = Math.atan(Math.sqrt(M * M - 1) / t1);
    const t3 = Math.atan(Math.sqrt(M * M - 1));
    return t1 * t2 - t3;
  }

  // Now for the inverse functions
  static findMFromMachAngle(angle: number, gamma: number): number {
    return 1 / Math.sin(angle);
  }
  static findMFrom_Tt_by_T(ratio: number, gamma: number): number {
    const m1sq = (2 / (gamma - 1)) * (ratio - 1);
    return Math.sqrt(m1sq);
  }
  static findMFrom_Pt_by_P(ratio: number, gamma: number): number {
    const pow = (gamma - 1) / gamma;
    const t_ratio = Math.pow(ratio, pow);
    return Isentropic.findMFrom_Tt_by_T(t_ratio, gamma);
  }
  static findMFrom_rhot_by_rho(ratio: number, gamma: number): number {
    const pow = gamma - 1;
    const t_ratio = Math.pow(ratio, pow);
    return Isentropic.findMFrom_Tt_by_T(t_ratio, gamma);
  }

  static findMFrom_A_by_Astar_mach_subsonic(
    ratio: number,
    gamma: number
  ): number {
    // Numerical solution
    let M = 0.9999; // Initial guess
    for (let i = 0; i < 10000; i++) {
      const f1 = Isentropic.A_by_Astar_mach(M, gamma);
      const f2 = Isentropic.A_by_Astar_mach(M - epsilon_derivative, gamma);
      const deriv_A_by_Astar = (f2 - f1) / epsilon_derivative;
      let dM = (step * (ratio - f1)) / deriv_A_by_Astar;
      dM = Math.min(dM, max_M_step);
      const M_new = M - dM;
      const new_ratio = Isentropic.A_by_Astar_mach(M_new, gamma);
      if (Math.abs(new_ratio - ratio) < tolerance) {
        return M_new;
      }
      M = M_new;
    }

    throw new Error("Could not find Mach number from A/A*");
  }
  static findMFrom_A_by_Astar_mach_supersonic(
    ratio: number,
    gamma: number
  ): number {
    // Numerical solution
    let M = 1.0001; // Initial guess
    for (let i = 0; i < MAX_ITER; i++) {
      const f1 = Isentropic.A_by_Astar_mach(M, gamma);
      const f2 = Isentropic.A_by_Astar_mach(M + epsilon_derivative, gamma);
      const deriv_A_by_Astar = (f2 - f1) / epsilon_derivative;
      let dM = (step * (ratio - f1)) / deriv_A_by_Astar;
      dM = Math.min(dM, max_M_step);
      const M_new = M + dM;
      const new_ratio = Isentropic.A_by_Astar_mach(M_new, gamma);
      if (Math.abs(new_ratio - ratio) < tolerance) {
        return M_new;
      }
      M = M_new;
    }
    throw new Error("Could not find Mach number from A/A*");
  }

  static get_ouputs(
    M: number,
    gamma: number
  ): {
    mach_angle: number;
    Tt_by_T: number;
    Pt_by_P: number;
    rhot_by_rho: number;
    at_by_a: number;
    tstar_by_t: number;
    pstar_by_p: number;
    rhostar_by_rho: number;
    A_by_Astar_mach: number;
    prandtl_meyer_angle: number;
  } {
    // let mach_angle = NaN
    // if(M >= 1){
    //   mach_angle = Math.asin(1/M);
    // }

    const mach_angle = Isentropic.getMachAngle(M);
    const Tt_by_t = Isentropic.Tt_by_T(M, gamma);
    const Pt_by_p = Isentropic.Pt_by_P(M, gamma);
    const rhot_by_rho = Isentropic.rhot_by_rho(M, gamma);
    const at_by_a = Isentropic.at_by_a(M, gamma);
    const tstar_by_t = Isentropic.tstar_by_t(M, gamma);
    const pstar_by_p = Isentropic.pstar_by_p(M, gamma);
    const rhostar_by_rho = Isentropic.rhostar_by_rho(M, gamma);
    const A_by_Astar_mach = Isentropic.A_by_Astar_mach(M, gamma);
    const prandtl_meyer_angle = Isentropic.prandtlMeyerAngle(M, gamma);
    return {
      mach_angle: mach_angle,
      Tt_by_T: Tt_by_t,
      Pt_by_P: Pt_by_p,
      rhot_by_rho: rhot_by_rho,
      at_by_a: at_by_a,
      tstar_by_t: tstar_by_t,
      pstar_by_p: pstar_by_p,
      rhostar_by_rho: rhostar_by_rho,
      A_by_Astar_mach: A_by_Astar_mach,
      prandtl_meyer_angle: prandtl_meyer_angle,
    };
  }
}

if (TEST_ISENTROPIC) {
  console.log("\nIsentropic test");
  const m1 = 0.4;
  console.log("Tt/T = ", Isentropic.Tt_by_T(m1, 1.4).toFixed(PRECISION));
  console.log("Pt/P = ", Isentropic.Pt_by_P(m1, 1.4).toFixed(PRECISION));
  console.log(
    "rhot/rho = ",
    Isentropic.rhot_by_rho(m1, 1.4).toFixed(PRECISION)
  );
  console.log("at/a = ", Isentropic.at_by_a(m1, 1.4).toFixed(PRECISION));
  // console.log("A/A* = ", Isentropic.A_by_Astar_mach(m1, 1.4).toFixed(PRECISION));
  const a_by_a_star = Isentropic.A_by_Astar_mach(m1, 1.4);
  console.log("A/A* = ", a_by_a_star.toFixed(PRECISION));
  console.log(
    "M from Tt/T = ",
    Isentropic.findMFrom_Tt_by_T(Isentropic.Tt_by_T(m1, 1.4), 1.4).toFixed(
      PRECISION
    )
  );
  console.log(
    "M from Pt/P = ",
    Isentropic.findMFrom_Pt_by_P(Isentropic.Pt_by_P(m1, 1.4), 1.4).toFixed(
      PRECISION
    )
  );
  console.log(
    "M from rhot/rho = ",
    Isentropic.findMFrom_rhot_by_rho(
      Isentropic.rhot_by_rho(m1, 1.4),
      1.4
    ).toFixed(PRECISION)
  );
  const supersonic_mach_solution =
    Isentropic.findMFrom_A_by_Astar_mach_supersonic(a_by_a_star, 1.4);
  console.log(
    "M from A/A* subsonic = ",
    Isentropic.findMFrom_A_by_Astar_mach_subsonic(a_by_a_star, 1.4).toFixed(
      PRECISION
    )
  );
  console.log(
    "M from A/A* supersonic = ",
    supersonic_mach_solution.toFixed(PRECISION)
  );
  console.log(
    "A/A* for supersonic solution of A/A* = ",
    Isentropic.A_by_Astar_mach(supersonic_mach_solution, 1.4).toFixed(PRECISION)
  );
  console.log("Isentropic test complete\n");
}

class NormalShock {
  static downStreamMachNumber(M1: number, gamma: number): number {
    const t1 = gamma - 1;
    const m1sq = M1 * M1;
    const m2sq = (t1 * m1sq + 2) / (2 * gamma * m1sq - t1);
    return Math.sqrt(m2sq);
  }
  static upStreamMachNumber(M2: number, gamma: number): number {
    const num = (gamma - 1) * M2 * M2 + 2;
    const den = 2 * gamma * M2 * M2 - (gamma - 1);
    return Math.sqrt(num / den);
  }
  static P2_by_P1(M1: number, gamma: number): number {
    const m1sq = M1 * M1;
    return (2 * gamma * m1sq - (gamma - 1)) / (gamma + 1);
  }
  static RHO2_by_RHO1(M1: number, gamma: number): number {
    const m1sq = M1 * M1;
    return ((gamma + 1) * m1sq) / ((gamma - 1) * m1sq + 2);
  }
  static T2_by_T1(M1: number, gamma: number): number {
    const p2byb1 = NormalShock.P2_by_P1(M1, gamma);
    const rho2byrho1 = NormalShock.RHO2_by_RHO1(M1, gamma);
    return p2byb1 / rho2byrho1;
  }
  static a2_by_a1(M1: number, gamma: number): number {
    return Math.sqrt(NormalShock.T2_by_T1(M1, gamma));
  }
  static Pt2_by_Pt1(M1: number, gamma: number): number {
    const t1 =
      ((gamma + 1) * M1 * M1 * 0.5) / (1 + (gamma - 1) * 0.5 * M1 * M1);
    const pow1 = gamma / (gamma - 1);
    const t2_1 = (2 * gamma * M1 * M1) / (gamma + 1);
    const t2_2 = -(gamma - 1) / (gamma + 1);
    const t2 = t2_1 + t2_2;
    const pow2 = -1 / (gamma - 1);
    const t1f = Math.pow(t1, pow1);
    const t2f = Math.pow(t2, pow2);
    return t1f * t2f;
  }

  static P1_by_Pt2(M1: number, gamma: number): number {
    return (
      1 /
      NormalShock.P2_by_P1(M1, gamma) /
      Isentropic.Pt_by_P(NormalShock.downStreamMachNumber(M1, gamma), gamma)
    );
  }

  // to allow for user flexibility, we find M1 by various other things. Inverse functions, basically
  // Need to manually calculate these
  static findM1From_P2_by_P1(ratio: number, gamma: number): number {
    const num = gamma - 1 + ratio * (gamma + 1);
    const den = 2 * gamma;
    return Math.sqrt(num / den);
  }
  static findM1From_RHO2_by_RHO1(ratio: number, gamma: number): number {
    const num = 2 * ratio;
    const den = gamma + 1 - ratio * (gamma - 1);
    return Math.sqrt(num / den);
  }
  static findM1From_T2_by_T1(ratio: number, gamma: number): number {
    const a = 2 * gamma * (gamma - 1);
    const b =
      4 * gamma - (gamma - 1) * (gamma - 1) - (gamma + 1) * (gamma + 1) * ratio;
    const c = -2 * (gamma - 1);
    const msq = (-b + Math.sqrt(b * b - 4 * a * c)) / (2 * a);
    return Math.sqrt(msq);
  }
  static findM1FromPt2_by_Pt1(ratio: number, gamma: number): number {
    let M = 1.0;
    let dM = max_M_step;
    let last_exceeded = false;
    for (let i = 0; i < MAX_ITER; i++) {
      const f1 = NormalShock.Pt2_by_Pt1(M, gamma);
      if (Math.abs(f1 - ratio) < tolerance) {
        return M;
      }
      if (f1 > ratio) {
        if (last_exceeded == false) {
          dM *= -1 / 2;
        }
        last_exceeded = true;
      } else {
        if (last_exceeded == true) {
          dM *= -1 / 2;
        }
        last_exceeded = false;
      }
      M -= dM;
    }
    return M;
  }

  static findM1FromP1_by_Pt2(ratio: number, gamma: number): number {
    let M = 1.0;
    let dM = max_M_step;
    let last_exceeded = false;
    for (let i = 0; i < MAX_ITER; i++) {
      const f1 = NormalShock.P1_by_Pt2(M, gamma);
      if (Math.abs(f1 - ratio) < tolerance) {
        return M;
      }
      if (f1 > ratio) {
        if (last_exceeded == false) {
          dM *= -1 / 2;
        }
        last_exceeded = true;
      } else {
        if (last_exceeded == true) {
          dM *= -1 / 2;
        }
        last_exceeded = false;
      }
      M -= dM;
    }
    return M;
  }

  /*
  I'm not doing p02_by_p01. That's too much work - requires numerical solution.
  */

  static get_ouputs(
    M1: number,
    gamma: number
  ): {
    M2: number;
    P2_by_P1: number;
    RHO2_by_RHO1: number;
    T2_by_T1: number;
    a2_by_a1: number;
    Pt2_by_Pt1: number;
    P1_by_Pt2: number;
  } {
    const M2 = NormalShock.downStreamMachNumber(M1, gamma);
    const P2_by_P1 = NormalShock.P2_by_P1(M1, gamma);
    const RHO2_by_RHO1 = NormalShock.RHO2_by_RHO1(M1, gamma);
    const T2_by_T1 = NormalShock.T2_by_T1(M1, gamma);
    const a2_by_a1 = NormalShock.a2_by_a1(M1, gamma);
    const Pt2_by_Pt1 = NormalShock.Pt2_by_Pt1(M1, gamma);
    const P1_by_Pt2 = NormalShock.P1_by_Pt2(M1, gamma);
    return {
      M2: M2,
      P2_by_P1: P2_by_P1,
      RHO2_by_RHO1: RHO2_by_RHO1,
      T2_by_T1: T2_by_T1,
      a2_by_a1: a2_by_a1,
      Pt2_by_Pt1: Pt2_by_Pt1,
      P1_by_Pt2: P1_by_Pt2,
    };
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
  console.log(
    "M1 from P2/P1 = ",
    NormalShock.findM1From_P2_by_P1(p2by1, gamma).toFixed(PRECISION)
  );
  console.log(
    "M1 from rho2/rho1 = ",
    NormalShock.findM1From_RHO2_by_RHO1(rho2by1, gamma).toFixed(PRECISION)
  );
  console.log(
    "M1 from T2/T1 = ",
    NormalShock.findM1From_T2_by_T1(t2by1, gamma).toFixed(PRECISION)
  );
  console.log("Mach2 Normal Shock test complete\n");
}

class ObliqueShock {
  static deflectionAngle(M1: number, beta: number, gamma: number): number {
    const num =
      2 *
      (Math.cos(beta) / Math.sin(beta)) *
      (Math.pow(M1 * Math.sin(beta), 2) - 1);
    const den = 2 + M1 * M1 * (gamma + Math.cos(2 * beta));
    return Math.atan(num / den);
  }
  static downStreamMachNumber(M1: number, beta: number, gamma: number): number {
    const deflection = ObliqueShock.deflectionAngle(M1, beta, gamma);
    const t0 = Math.pow(Math.sin(beta - deflection), 2);
    const num = Math.pow(M1 * Math.sin(beta), 2) + 2 / (gamma - 1);
    const den =
      ((2 * gamma) / (gamma - 1)) * Math.pow(M1 * Math.sin(beta), 2) - 1;
    return Math.sqrt(num / (den * t0));
  }
  static P2_by_P1(M1: number, beta: number, gamma: number): number {
    const num = 2 * gamma * Math.pow(M1 * Math.sin(beta), 2) - (gamma - 1);
    const den = gamma + 1;
    return num / den;
  }
  static RHO2_by_RHO1(M1: number, beta: number, gamma: number): number {
    const num = (gamma + 1) * Math.pow(M1 * Math.sin(beta), 2);
    const den = (gamma - 1) * Math.pow(M1 * Math.sin(beta), 2) + 2;
    return num / den;
  }
  static T2_by_T1(M1: number, beta: number, gamma: number): number {
    const p2byb1 = ObliqueShock.P2_by_P1(M1, beta, gamma);
    const rho2byrho1 = ObliqueShock.RHO2_by_RHO1(M1, beta, gamma);
    return p2byb1 / rho2byrho1;
  }
  static find_beta_from_M1_and_M1n1(M1: number, Mn1: number): number {
    return Math.asin(Mn1 / M1);
  }
  static find_max_deflection_angle(M1: number, gamma: number): number {
    const t1 = (gamma + 1) / (4 * gamma);
    const t2_0 = -1 / (gamma * M1 * M1);
    const t2_11 = gamma + 1;
    const t2_12 =
      1 + 0.5 * (gamma - 1) * M1 * M1 + ((gamma + 1) * Math.pow(M1, 4)) / 16;
    const t2_1 = Math.sqrt(t2_11 * t2_12);
    const t = t1 + t2_0 * (1 - t2_1);
    return Math.asin(Math.sqrt(t));
  }

  static find_shock_angle_solutions(
    M1: number,
    gamma: number,
    deflection: number
  ): [number, number] {
    // const lhs = Math.tan(deflection);
    // const max_deflection = ObliqueShock.find_max_deflection_angle(M1, gamma);
    // if( )

    let P1 = Math.PI / 2 - epsilon;
    let P2 = Math.asin(1 / M1) + epsilon;
    let P1_final = P1;
    let P2_final = P2;
    for (let i = 0; i < MAX_ITER; i++) {
      const f1 = ObliqueShock.deflectionAngle(M1, P1, gamma);
      const f2 = ObliqueShock.deflectionAngle(M1, P1 - epsilon, gamma);
      const deriv_deflection = (f2 - f1) / epsilon;
      const P1_new = P1 - (step * (deflection - f1)) / deriv_deflection;
      const new_deflection = ObliqueShock.deflectionAngle(M1, P1_new, gamma);
      if (Math.abs(new_deflection - deflection) < tolerance) {
        P1_final = P1;
        break;
      }
      P1 = P1_new;
    }
    for (let i = 0; i < MAX_ITER; i++) {
      const f1 = ObliqueShock.deflectionAngle(M1, P2, gamma);
      const f2 = ObliqueShock.deflectionAngle(M1, P2 + epsilon, gamma);
      const deriv_deflection = (f2 - f1) / epsilon;
      const P2_new = P2 + (step * (deflection - f1)) / deriv_deflection;
      const new_deflection = ObliqueShock.deflectionAngle(M1, P2_new, gamma);
      if (Math.abs(new_deflection - deflection) < tolerance) {
        P2_final = P2;
        break;
      }
      P2 = P2_new;
    }
    if (Math.abs(P1_final - P2_final) < tolerance) {
      const average = (P1_final + P2_final) / 2;
      P1_final = average;
      P2_final = average;
    }
    return [P1_final, P2_final];
  }
  static findStrongWeakSolutions(
    M1: number,
    gamma: number,
    deflection: number
  ): [number[], number[], number[], number[], number[]] {
    const max_deflection = ObliqueShock.find_max_deflection_angle(M1, gamma);
    if (deflection > max_deflection) {
      throw new Error(
        "Deflection angle exceeds maximum deflection angle for given Mach number"
      );
    }
    const betas = ObliqueShock.find_shock_angle_solutions(
      M1,
      gamma,
      deflection
    );
    const M2s = [
      ObliqueShock.downStreamMachNumber(M1, betas[1], gamma),
      ObliqueShock.downStreamMachNumber(M1, betas[0], gamma),
    ];
    const P2_by_P1s = [
      ObliqueShock.P2_by_P1(M1, betas[1], gamma),
      ObliqueShock.P2_by_P1(M1, betas[0], gamma),
    ];
    const rho2_by_rho1s = [
      ObliqueShock.RHO2_by_RHO1(M1, betas[1], gamma),
      ObliqueShock.RHO2_by_RHO1(M1, betas[0], gamma),
    ];
    const T2_by_T1s = [
      ObliqueShock.T2_by_T1(M1, betas[1], gamma),
      ObliqueShock.T2_by_T1(M1, betas[0], gamma),
    ];

    return [betas, M2s, P2_by_P1s, rho2_by_rho1s, T2_by_T1s];
  }
  static normalMachs(
    M1: number,
    gamma: number,
    deflection: number
  ): [[number, number], [number, number]] {
    const betas = ObliqueShock.find_shock_angle_solutions(
      M1,
      gamma,
      deflection
    );
    const M1n: [number, number] = [
      M1 * Math.sin(betas[0]),
      M1 * Math.sin(betas[1]),
    ];
    const M2n: [number, number] = [
      ObliqueShock.downStreamMachNumber(M1, betas[0], gamma) *
        Math.sin(betas[0] - deflection),
      ObliqueShock.downStreamMachNumber(M1, betas[1], gamma) *
        Math.sin(betas[1] - deflection),
    ];
    return [M1n, M2n];
  }

  static get_outputs(
    M1: number,
    beta: number,
    gamma: number
  ): {
    M2: number;
    deflection: number;
    P2_by_P1: number;
    RHO2_by_RHO1: number;
    T2_by_T1: number;
    max_deflection: number;
  } {
    const M2 = ObliqueShock.downStreamMachNumber(M1, beta, gamma);

    const deflection = ObliqueShock.deflectionAngle(M1, beta, gamma);
    const P2_by_P1 = ObliqueShock.P2_by_P1(M1, beta, gamma);
    const RHO2_by_RHO1 = ObliqueShock.RHO2_by_RHO1(M1, beta, gamma);
    const T2_by_T1 = ObliqueShock.T2_by_T1(M1, beta, gamma);
    const max_deflection = ObliqueShock.find_max_deflection_angle(M1, gamma);

    return { M2, deflection, P2_by_P1, RHO2_by_RHO1, T2_by_T1, max_deflection };
  }
}

if (TEST_OBLIQUE_SHOCK) {
  console.log("\nOblique Shock test");
  const M1 = 2.0;
  const gamma = 1.4;
  const deflection = (20 * Math.PI) / 180;
  const [beta_weak, beta_strong] = ObliqueShock.find_shock_angle_solutions(
    M1,
    gamma,
    deflection
  );
  console.log("M1 = ", M1);
  console.log(
    "Deflection (deg) = ",
    ((deflection * 180) / Math.PI).toFixed(PRECISION)
  );
  console.log(
    "Weak Shock Angle (deg) = ",
    ((beta_weak * 180) / Math.PI).toFixed(PRECISION)
  );
  console.log(
    "Strong Shock Angle (deg) = ",
    ((beta_strong * 180) / Math.PI).toFixed(PRECISION)
  );
}

const updatePrecision = (): void => {
  PRECISION =
    parseInt(
      (document.getElementById("precision_input") as HTMLInputElement).value
    ) || 4;
};
const updateGamma = (): void => {
  γ =
    parseFloat(
      (document.getElementById("gamma_input") as HTMLInputElement).value
    ) || 1.4;
};
const updateTolerance = (): void => {
  tolerance =
    parseFloat(
      (document.getElementById("tolerance_input") as HTMLInputElement).value
    ) || 0.0001;
};

function calculateIsentropic() {
  const input = (
    document.getElementById("input_isentropic_number") as HTMLInputElement
  ).value;
  const type = (
    document.getElementById("isentropic_input_type") as HTMLSelectElement
  ).value;
  const output = document.getElementById(
    "isentropic_output_area"
  ) as HTMLDivElement;

  if (!input) {
    output.textContent = "Please enter a value";
    return;
  }

  const value_input = parseFloat(input);
  updateGamma();
  updatePrecision();
  updateTolerance();

  try {
    let mach: number;
    switch (type) {
      case "Mach":
        mach = value_input;
        break;
      case "Mach_Angle":
        mach = Isentropic.findMFromMachAngle((value_input * Math.PI) / 180, γ);
        break;
      case "Tt_by_T":
        mach = Isentropic.findMFrom_Tt_by_T(value_input, γ);
        break;
      case "Pt_by_P":
        mach = Isentropic.findMFrom_Pt_by_P(value_input, γ);
        break;
      case "rhot_by_rho":
        mach = Isentropic.findMFrom_rhot_by_rho(value_input, γ);
        break;
      case "A_by_Astar_sub":
        mach = Isentropic.findMFrom_A_by_Astar_mach_subsonic(value_input, γ);
        break;
      case "A_by_Astar_sup":
        mach = Isentropic.findMFrom_A_by_Astar_mach_supersonic(value_input, γ);
        break;
      default:
        throw new Error("Invalid input type");
    }

    const results = Isentropic.get_ouputs(mach, γ);
    output.textContent = `Mach Number: ${mach.toFixed(PRECISION)}
    
All Isentropic Relations:
Mach Angle = ${
      isNaN(results.mach_angle)
        ? "N/A (Subsonic)"
        : ((results.mach_angle * 180) / Math.PI).toFixed(PRECISION) + "°"
    }
Tt/T = ${results.Tt_by_T.toFixed(PRECISION)}
Pt/P = ${results.Pt_by_P.toFixed(PRECISION)}
ρt/ρ = ${results.rhot_by_rho.toFixed(PRECISION)}
at/a = ${results.at_by_a.toFixed(PRECISION)}
T*/T = ${results.tstar_by_t.toFixed(PRECISION)}
P*/P = ${results.pstar_by_p.toFixed(PRECISION)}
ρ*/ρ = ${results.rhostar_by_rho.toFixed(PRECISION)}
A/A* = ${results.A_by_Astar_mach.toFixed(PRECISION)}
Prandtl-Meyer Angle = ${((results.prandtl_meyer_angle * 180) / Math.PI).toFixed(
      PRECISION
    )}°`;
  } catch (error) {
    output.textContent = `Error: ${
      error instanceof Error ? error.message : String(error)
    }`;
  }
}

function calculateNormalShock() {
  const input = (
    document.getElementById("input_normalshock_number") as HTMLInputElement
  ).value;
  const type = (
    document.getElementById("normalshock_input_type") as HTMLSelectElement
  ).value;
  const output = document.getElementById(
    "normalshock_output_area"
  ) as HTMLDivElement;

  if (!input) {
    output.textContent = "Please enter a value";
    return;
  }

  const user_input = parseFloat(input);
  updateGamma();
  updatePrecision();
  updateTolerance();

  try {
    let m1: number;
    switch (type) {
      case "Mach_Up":
        m1 = user_input;
        break;
      case "Mach_Down":
        m1 = NormalShock.upStreamMachNumber(user_input, γ);
        break;
      case "P2_by_P1":
        m1 = NormalShock.findM1From_P2_by_P1(user_input, γ);
        break;
      case "Pt2_by_Pt1":
        m1 = NormalShock.findM1FromPt2_by_Pt1(user_input, γ);
        break;
      case "P1_by_Pt2":
        m1 = NormalShock.findM1FromP1_by_Pt2(user_input, γ);
        break;
      case "RHO2_by_RHO1":
        m1 = NormalShock.findM1From_RHO2_by_RHO1(user_input, γ);
        break;
      case "T2_by_T1":
        m1 = NormalShock.findM1From_T2_by_T1(user_input, γ);
        break;
      default:
        throw new Error("Invalid input type");
    }

    const results = NormalShock.get_ouputs(m1, γ);
    output.textContent = `Upstream Mach (M1): ${m1.toFixed(PRECISION)}

All Normal Shock Relations:
M2 = ${results.M2.toFixed(PRECISION)}
P2/P1 = ${results.P2_by_P1.toFixed(PRECISION)}
ρ2/ρ1 = ${results.RHO2_by_RHO1.toFixed(PRECISION)}
T2/T1 = ${results.T2_by_T1.toFixed(PRECISION)}
a2/a1 = ${results.a2_by_a1.toFixed(PRECISION)}
Pt2/Pt1 = ${results.Pt2_by_Pt1.toFixed(PRECISION)}
P1/Pt2 = ${results.P1_by_Pt2.toFixed(PRECISION)}`;
  } catch (error) {
    output.textContent = `Error: ${
      error instanceof Error ? error.message : String(error)
    }`;
  }
}

function calculateObliqueShock() {
  const m1Input = (
    document.getElementById("input_obliqueshock_M") as HTMLInputElement
  ).value;
  const secondInput = (
    document.getElementById("input_obliqueshock_2") as HTMLInputElement
  ).value;
  const type = (
    document.getElementById("obliqueshock_input_2_type") as HTMLSelectElement
  ).value;
  const output = document.getElementById(
    "obliqueshock_output_area"
  ) as HTMLDivElement;

  if (!m1Input || !secondInput) {
    output.textContent = "Please enter both values";
    return;
  }

  const m1 = parseFloat(m1Input);
  const secondValue = parseFloat(secondInput);
  updateGamma();
  updatePrecision();
  updateTolerance();
  try {
    if (type === "beta") {
      const beta = (secondValue * Math.PI) / 180;
      const results = ObliqueShock.get_outputs(m1, beta, γ);
      output.textContent = `Oblique Shock Results:
M1 = ${m1.toFixed(PRECISION)}
Beta (shock angle) = ${secondValue.toFixed(2)}°
Deflection angle = ${((results.deflection * 180) / Math.PI).toFixed(PRECISION)}°

M2 = ${results.M2.toFixed(PRECISION)}
P2/P1 = ${results.P2_by_P1.toFixed(PRECISION)}
ρ2/ρ1 = ${results.RHO2_by_RHO1.toFixed(PRECISION)}
T2/T1 = ${results.T2_by_T1.toFixed(PRECISION)}
Max deflection = ${((results.max_deflection * 180) / Math.PI).toFixed(
        PRECISION
      )}°`;
    } else if (type === "deflection") {
      const deflection = (secondValue * Math.PI) / 180;
      const [betas, M2s, P2_by_P1s, rho2_by_rho1s, T2_by_T1s] =
        ObliqueShock.findStrongWeakSolutions(m1, γ, deflection);
      output.textContent = `Oblique Shock Solutions:
M1 = ${m1.toFixed(PRECISION)}
Deflection angle = ${secondValue.toFixed(2)}°

Weak Solution:
Beta = ${((betas[0] * 180) / Math.PI).toFixed(PRECISION)}°
M2 = ${M2s[0].toFixed(PRECISION)}
P2/P1 = ${P2_by_P1s[0].toFixed(PRECISION)}
ρ2/ρ1 = ${rho2_by_rho1s[0].toFixed(PRECISION)}
T2/T1 = ${T2_by_T1s[0].toFixed(PRECISION)}

Strong Solution:
Beta = ${((betas[1] * 180) / Math.PI).toFixed(PRECISION)}°
M2 = ${M2s[1].toFixed(PRECISION)}
P2/P1 = ${P2_by_P1s[1].toFixed(PRECISION)}
ρ2/ρ1 = ${rho2_by_rho1s[1].toFixed(PRECISION)}
T2/T1 = ${T2_by_T1s[1].toFixed(PRECISION)}`;
    }
  } catch (error) {
    output.textContent = `Error: ${
      error instanceof Error ? error.message : String(error)
    }`;
  }
}
