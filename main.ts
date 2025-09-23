
class NormalShock{
  static DownStreamMachNumber(M1:number, gamma:number):number {
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
  static A2_by_A1(M1:number, gamma: number):number {
    return Math.sqrt(NormalShock.T2_by_T1(M1, gamma));
  }


}