const PRECISION = 5;
const TEST_ISENTROPIC = false;
const TEST_NORMAL_SHOCK = true;
const TEST_OBLIQUE_SHOCK = false;
class Isentropic{
  static Tt_by_T(M:number, gamma:number):number {
    return 1 + ((gamma - 1) * 0.5)*M*M;
  }
  static Pt_by_P(M:number, gamma:number):number {
    const ttr = Isentropic.Tt_by_T(M, gamma);
    return Math.pow(ttr, gamma/(gamma - 1));
  }
  static rhot_by_rho(M:number, gamma:number):number {
    const ttr = Isentropic.Tt_by_T(M, gamma);
    return Math.pow(ttr, 1/(gamma - 1));
  }
  static at_by_a(M:number, gamma:number):number {
    const ttr = Isentropic.Tt_by_T(M, gamma);
    return Math.sqrt(ttr);
  }
  static tstar_by_t(M:number, gamma:number):number{
    return 2/(gamma + 1);
  }
  static pstar_by_p(M:number, gamma:number):number{
    const tstar = this.tstar_by_t(M, gamma);
    return Math.pow(tstar, gamma/(gamma - 1));
  }
  static rhostar_by_rho(M:number, gamma:number):number{
    const tstar = this.tstar_by_t(M, gamma);
    return Math.pow(tstar, 1/(gamma - 1));
  }
  static A_by_Astar_pressure(M:number, gamma:number, p:number, p0:number):number{
    const gp = gamma+1;
    const gn = gamma-1;
    const pr = p/p0;
    const num = Math.pow(2/gp, 0.5 * gp/gn) * Math.pow(0.5 * gn, 0.5);
    const den = Math.pow(Math.pow(pr, 2/gamma) - Math.pow(pr, (gp)/gamma), 0.5);
    return num/den;
  }
  static A_by_Astar_mach(M:number, gamma: number):number{
    const pow = (gamma + 1)/(2*(gamma - 1));
    const t1 = (2/(gamma + 1));
    const t2 = (1 + 0.5*(gamma - 1)*M*M);
    const i1 = Math.pow(t1 * t2, pow);
    return i1/M;
  }
}


class NormalShock{
  static downStreamMachNumber(M1:number, gamma:number):number {

    const t1 = (gamma - 1);
    const m1sq = M1*M1;
    const m2sq = (t1*m1sq + 2)/(2*gamma*m1sq - t1);
    return Math.sqrt(m2sq);
  }
  static P2_by_P1(M1:number, gamma:number):number {
    const m1sq = M1*M1;
    return (2*gamma*m1sq - (gamma - 1))/(gamma + 1);
  }
  static RHO2_by_RHO1(M1:number, gamma:number):number {
    const m1sq = M1*M1;
    return ((gamma + 1)*m1sq)/((gamma - 1)*m1sq + 2);
  }
  static T2_by_T1(M1:number, gamma:number):number {
    const p2byb1 = NormalShock.P2_by_P1(M1, gamma);
    const rho2byrho1 = NormalShock.RHO2_by_RHO1(M1, gamma);
    return p2byb1/rho2byrho1;
  }
  static a2_by_a1(M1:number, gamma: number):number {
    return Math.sqrt(NormalShock.T2_by_T1(M1, gamma));
  }
  static Pt2_by_Pt1(M1:number, gamma:number):number {
    const t1 = (gamma +1)*M1*M1*0.5/ (1 + ((gamma -1)*0.5)*M1*M1);
    const pow1 = gamma/(gamma -1);
    const t2_1 = 2*gamma*M1*M1/(gamma+1);
    const t2_2 = -(gamma-1)/(gamma+1);
    const t2 = t2_1 + t2_2;
    const pow2 = -1/(gamma -1);
    const t1f = Math.pow(t1, pow1)
    const t2f = Math.pow(t2, pow2);
    return t1f * t2f;
  }

  static P1_by_Pt2(M1:number, gamma:number):number {
    return (1/(NormalShock.P2_by_P1(M1, gamma)) / Isentropic.Pt_by_P(NormalShock.downStreamMachNumber(M1, gamma), gamma));
  }
  
}

if(TEST_ISENTROPIC){
  console.log("\nMach2 Normal Shock test");
  console.log("M2 = ", NormalShock.downStreamMachNumber(2, 1.4).toFixed(PRECISION));
  console.log("P2/P1 = ", NormalShock.P2_by_P1(2, 1.4).toFixed(PRECISION));
  console.log("rho2/rho1 = ", NormalShock.RHO2_by_RHO1(2, 1.4).toFixed(PRECISION));
  console.log("T2/T1 = ", NormalShock.T2_by_T1(2, 1.4).toFixed(PRECISION));
  console.log("a2/a1 = ", NormalShock.a2_by_a1(2, 1.4).toFixed(PRECISION));
  console.log("Pt2/Pt1 = ", NormalShock.Pt2_by_Pt1(2, 1.4).toFixed(PRECISION));
  console.log("P1/Pt2 = ", NormalShock.P1_by_Pt2(2, 1.4).toFixed(PRECISION));
  console.log("Mach2 Normal Shock test complete\n")
}


class ObliqueShock{
  static deflectionAngle(M1:number, beta:number, gamma: number):number{
    const num = 2* (Math.cos(beta)/Math.sin(beta)) * (Math.pow(M1 * Math.sin(beta), 2) - 1);
    const den = 2 + M1*M1 * (gamma + Math.cos(2*beta));
    return Math.atan(num/den);
  }
  static downStreamMachNumber(M1: number, beta:number, gamma: number):number{
    const deflection = ObliqueShock.deflectionAngle(M1, beta, gamma);
    const t0 = Math.pow(Math.sin(beta-deflection), 2);
    const num = Math.pow(M1 * Math.sin(beta), 2) + (2/(gamma - 1));
    const den = (2*gamma/(gamma - 1)) * Math.pow(M1 * Math.sin(beta), 2) - 1;
    return Math.sqrt(num/(den*t0));
  }
  static P2_by_P1(M1:number, beta:number, gamma:number):number {
    const num = 2 * gamma * Math.pow(M1 * Math.sin(beta), 2) - (gamma - 1);
    const den = gamma + 1;
    return num/den;
  }
  static RHO2_by_RHO1(M1:number, beta:number, gamma:number):number {
    const num = (gamma + 1) * Math.pow(M1 * Math.sin(beta), 2);
    const den = (gamma - 1) * Math.pow(M1 * Math.sin(beta), 2) + 2;
    return num/den;
  }
  static T2_by_T1(M1:number, beta:number, gamma:number):number {
    const p2byb1 = ObliqueShock.P2_by_P1(M1, beta, gamma);
    const rho2byrho1 = ObliqueShock.RHO2_by_RHO1(M1, beta, gamma);
    return p2byb1/rho2byrho1;
  }
}



