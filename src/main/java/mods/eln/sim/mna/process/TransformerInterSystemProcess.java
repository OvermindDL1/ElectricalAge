package mods.eln.sim.mna.process;

import mods.eln.sim.mna.SubSystem.Th;
import mods.eln.sim.mna.component.CurrentSource;
import mods.eln.sim.mna.component.VoltageSource;
import mods.eln.sim.mna.misc.IRootSystemPreStepProcess;
import mods.eln.sim.mna.state.State;

public class TransformerInterSystemProcess implements IRootSystemPreStepProcess {
    State primaryState, secondaryState;
    VoltageSource primarySource, secondarySource;

    double ratio = 1;

    public TransformerInterSystemProcess(State primaryState, State secondaryState, VoltageSource primarySource, VoltageSource secondarySource) {
        this.primaryState = primaryState;
        this.secondaryState = secondaryState;
        this.primarySource = primarySource;
        this.secondarySource = secondarySource;
    }

    @Override
    public void rootSystemPreStepProcess() {
        /*  Math equations for the following code:

        Primary Voltage: V1
        Primary Resistance: R1
        Secondary Voltage: V2
        Secondary Resistance: R2
        Voltage midpoint at transformer voltage sources: Vm
        Conversion ratio: r

        Circuit:
            (V1)---[ R1 ]---(Vm) | (r Vm)---[ R2 ]---(V2)

        Circuit, scaling the primary side so that the voltages on each side match:
            (rV1)---[ r^2 R1 ]---(r Vm) | (r Vm)---[ R2 ]---(V2)

        Total Voltage Drop (Vd) = r V1 - V2
        Total Resistance (Rt) = r^2 R1 + R2
        Current (i) = Vd / Rt

        Voltage Drop across primary resistor (r Vr1) = i r^2 R1
            Vri = i r R1
        Voltage Drop across secondary resistor (Vr2) = i R2

        Vm on primary side = V1 - Vr1 = V1 - i r R1
        Vm on secondary side = Vr2 + V2 = V2 + i R2
        */

        Th primaryEquivalent = primaryState.getSubSystem().getTh(primaryState, primarySource);
        Th secondaryEquivalent = secondaryState.getSubSystem().getTh(secondaryState, secondarySource);

        if (primaryEquivalent.isHighImpedance() || secondaryEquivalent.isHighImpedance()) {
            primarySource.setU(0);
            secondarySource.setU(0);
            return;
        }

        double current = (primaryEquivalent.U * ratio - secondaryEquivalent.U) /
                (primaryEquivalent.R * ratio * ratio + secondaryEquivalent.R);
        primarySource.setU(primaryEquivalent.U - current * ratio * primaryEquivalent.R);
        secondarySource.setU(secondaryEquivalent.U + current * secondaryEquivalent.R);
    }

    public void setRatio(double ratio) {
        this.ratio = ratio;
    }

    public double getRatio() {
        return this.ratio;
    }
}
