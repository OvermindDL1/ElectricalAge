package mods.eln.sim.mna.process;

import mods.eln.sim.mna.SubSystem.Th;
import mods.eln.sim.mna.component.CurrentSource;
import mods.eln.sim.mna.component.VoltageSource;
import mods.eln.sim.mna.misc.IRootSystemPreStepProcess;
import mods.eln.sim.mna.state.State;

public class TransformerInterSystemProcess implements IRootSystemPreStepProcess {
    State primaryState, secondaryState;
    CurrentSource primarySource;
    VoltageSource secondarySource;

    double ratio = 1;

    public TransformerInterSystemProcess(State primaryState, State secondaryState, CurrentSource primarySource, VoltageSource secondarySource) {
        this.primaryState = primaryState;
        this.secondaryState = secondaryState;
        this.primarySource = primarySource;
        this.secondarySource = secondarySource;
    }

    @Override
    public void rootSystemPreStepProcess() {
        Th primaryEquivalent = primaryState.getSubSystem().getTh(primaryState, primarySource);
        Th secondaryEquivalent = secondaryState.getSubSystem().getTh(secondaryState, secondarySource);

        if (primaryEquivalent.isHighImpedance() && secondaryEquivalent.isHighImpedance()) {
            primarySource.setCurrent(0);
            secondarySource.setU(0);
        } else if (primaryEquivalent.isHighImpedance()) {
            primarySource.setCurrent(0);
            secondarySource.setU(secondaryEquivalent.U);
        } else if (secondaryEquivalent.isHighImpedance()) {
            secondarySource.setU(primaryEquivalent.U * ratio);
            primarySource.setCurrent(secondarySource.getSubSystem().solve(secondarySource.getCurrentState()) * ratio);
        } else {
            primaryEquivalent.U *= ratio;
            primaryEquivalent.R *= ratio * ratio;

            double Vd = secondaryEquivalent.U - primaryEquivalent.U;
            double Rt = primaryEquivalent.R + secondaryEquivalent.R;
            primarySource.setCurrent(Vd / Rt * ratio);
            secondarySource.setU(-secondaryEquivalent.R * Vd / Rt + secondaryEquivalent.U);
        }
    }

    public void setRatio(double ratio) {
        this.ratio = ratio;
    }

    public double getRatio() {
        return this.ratio;
    }
}
