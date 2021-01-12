package mods.eln.sim.mna;

import mods.eln.misc.Profiler;
import mods.eln.misc.Utils;
import mods.eln.sim.mna.component.*;
import mods.eln.sim.mna.misc.*;
import mods.eln.sim.mna.process.TransformerInterSystemProcess;
import mods.eln.sim.mna.state.State;
import mods.eln.sim.mna.state.VoltageState;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class SubSystem {
    public ArrayList<Component> component = new ArrayList<Component>();
    public List<State> states = new ArrayList<State>();
    public LinkedList<IDestructor> breakDestructor = new LinkedList<IDestructor>();
    public ArrayList<SubSystem> interSystemConnectivity = new ArrayList<SubSystem>();
    ArrayList<ISubSystemProcessI> processI = new ArrayList<ISubSystemProcessI>();
    State[] statesTab;

    RootSystem root;

    double dt;
    boolean matrixValid = false;

    int stateCount;
    RealMatrix A;
    //RealMatrix I;
    boolean singularMatrix;

    double[][] AInvdata;
    double[] Idata;

    double[] XtempData;

    boolean breaked = false;

    ArrayList<ISubSystemProcessFlush> processF = new ArrayList<ISubSystemProcessFlush>();

    public RootSystem getRoot() {
        return root;
    }

    public SubSystem(RootSystem root, double dt) {
        this.dt = dt;
        this.root = root;
    }

    public void invalidate() {
        matrixValid = false;
    }

    public void addComponent(Component c) {
        component.add(c);
        c.addedTo(this);
        invalidate();
    }

    public void addState(State s) {
        states.add(s);
        s.addedTo(this);
        invalidate();
    }

    public void removeComponent(Component c) {
        component.remove(c);
        c.quitSubSystem();
        invalidate();
    }

    public void removeState(State s) {
        states.remove(s);
        s.quitSubSystem();
        invalidate();
    }

	/*public void removeAll() {
		for (Component c : component) {
			c.disconnectFromSubSystem();
		}
		for (State s : states) {
			s.disconnectFromSubSystem();
		}	
		invalidate();
	}*/

    public void removeProcess(ISubSystemProcessI p) {
        processI.remove(p);
        invalidate();
    }

    public void addComponent(Iterable<Component> i) {
        for (Component c : i) {
            addComponent(c);
        }
    }

    public void addState(Iterable<State> i) {
        for (State s : i) {
            addState(s);
        }
    }

    public void addProcess(ISubSystemProcessI p) {
        processI.add(p);
    }

    //double[][] getDataRef()

    public void generateMatrix() {
        stateCount = states.size();

        Profiler p = new Profiler();
        p.add("Inversse with " + stateCount + " state : ");

        A = MatrixUtils.createRealMatrix(stateCount, stateCount);
        //Adata = ((Array2DRowRealMatrix) A).getDataRef();
        // X = MatrixUtils.createRealMatrix(stateCount, 1); Xdata =
        // ((Array2DRowRealMatrix)X).getDataRef();
        //I = MatrixUtils.createRealMatrix(stateCount, 1);
        //Idata = ((Array2DRowRealMatrix) I).getDataRef();
        Idata = new double[stateCount];
        XtempData = new double[stateCount];
        {
            int idx = 0;
            for (State s : states) {
                s.setId(idx++);
            }
        }

        for (Component c : component) {
            c.applyTo(this);
        }

        //	org.apache.commons.math3.linear.

        try {
            //FieldLUDecomposition QRDecomposition  LUDecomposition RRQRDecomposition
            RealMatrix Ainv = new QRDecomposition(A).getSolver().getInverse();
            AInvdata = Ainv.getData();
            singularMatrix = false;
        } catch (Exception e) {
            singularMatrix = true;
            if (stateCount > 1) {
                int idx = 0;
                idx++;
                Utils.println("//////////SingularMatrix////////////");
            }
        }

        statesTab = new State[stateCount];
        statesTab = states.toArray(statesTab);

        matrixValid = true;

        p.stop();
        Utils.println(p);
    }

    public void addToA(State a, State b, double v) {
        if (a == null || b == null)
            return;
        A.addToEntry(a.getId(), b.getId(), v);
        //Adata[a.getId()][b.getId()] += v;
    }

    public void addToI(State s, double v) {
        if (s == null) return;
        Idata[s.getId()] = v;
        //Idata[s.getId()][0] += v;
    }

	/*
	 * public void pushX(){
	 * 
	 * }
	 */
	/*
	 * public void popX(){
	 * 
	 * }
	 */

    public void step() {
        stepCalc();
        stepFlush();
    }

    public void stepCalc() {
        //Profiler profiler = new Profiler();
        //profiler.add("generateMatrix");
        if (!matrixValid) {
            generateMatrix();
        }

        if (!singularMatrix) {
            //profiler.add("generateMatrix");
            for (int y = 0; y < stateCount; y++) {
                Idata[y] = 0;
            }
            //profiler.add("generateMatrix");
            for (ISubSystemProcessI p : processI) {
                p.simProcessI(this);
            }
            //profiler.add("generateMatrix");

            for (int idx2 = 0; idx2 < stateCount; idx2++) {
                double stack = 0;
                for (int idx = 0; idx < stateCount; idx++) {
                    stack += AInvdata[idx2][idx] * Idata[idx];
                }
                XtempData[idx2] = stack;
            }
            //Xtemp = Ainv.multiply(I);
        }
        //profiler.stop();
        //Utils.println(profiler);
    }

    public double solve(State pin) {
        if (!matrixValid) {
            generateMatrix();
        }

        if (!singularMatrix) {
            for (int y = 0; y < stateCount; y++)
                Idata[y] = 0;

            for (ISubSystemProcessI p : processI)
                p.simProcessI(this);

            int idx2 = pin.getId();
            double stack = 0;
            for (int idx = 0; idx < stateCount; idx++)
                stack += AInvdata[idx2][idx] * Idata[idx];

            return stack;
        }
        return 0;
    }

    //RealMatrix Xtemp;
    public void stepFlush() {
        if (!singularMatrix) {
            for (int idx = 0; idx < stateCount; idx++) {
                //statesTab[idx].state = Xtemp.getEntry(idx, 0);
                statesTab[idx].state = XtempData[idx];

            }
        } else {
            for (int idx = 0; idx < stateCount; idx++) {
                statesTab[idx].state = 0;
            }
        }

        for (ISubSystemProcessFlush p : processF) {
            p.simProcessFlush();
        }
    }

    public static void main(String[] args) {
        RootSystem root = new RootSystem(0.05, 1000);
        SubSystem ss1 = new SubSystem(root, 0.05);
        SubSystem ss2 = new SubSystem(root, 0.05);

        root.systems.add(ss1);
        root.systems.add(ss2);

        State s1 = new State(),
            s2 = new State(),
            s3 = new State(),
            s4 = new State();
        ss1.addState(s1);
        ss1.addState(s2);
        ss2.addState(s3);
        ss2.addState(s4);

        VoltageSource e1 = new VoltageSource("e1", s1, null).setU(5);
        VoltageSource e2 = new VoltageSource("e2", null, s4).setU(0);
        Resistor r1 = new Resistor().setR(100);
        Resistor r2 = new Resistor().setR(300);

        ss1.addComponent(e1);
        ss1.addComponent(r1.connectTo(s1,s2));
        ss2.addComponent(e2);
        ss2.addComponent(r2.connectTo(s3,s4));

        CurrentSource magicIn = new CurrentSource("magicIn", s2, null).setCurrent(0);
        VoltageSource magicOut = new VoltageSource("magicOut", s3, null).setU(0);

        ss1.addComponent(magicIn);
        ss2.addComponent(magicOut);

        double ratio = 1;
        root.addProcess((IRootSystemPreStepProcess) () -> {
            Th thIn = ss1.getTh(s2, magicIn);
            Th thOut = ss2.getTh(s3, magicOut);

            if (thIn.isHighImpedance() || thOut.isHighImpedance()) {
                magicIn.setCurrent(0);
                magicOut.setU(0);
            } else {
                thIn.U *= ratio;
                thIn.R *= ratio * ratio;

                double Vd = thOut.U - thIn.U;
                double Rt = thIn.R + thOut.R;
                magicIn.setCurrent(Vd / Rt * ratio);
                magicOut.setU(-thOut.R * Vd / Rt + thOut.U);
            }
        });

        root.step();

        System.out.println("e1: V = " + e1.getU() + ", I = " + e1.getI());
        System.out.println("s1: V = " + s1.state);
        System.out.println("r1: Vd = " + r1.getU() + ", I = " + r1.getI());
        System.out.println("s2: V = " + s2.state);
        System.out.println("magicIn: V = " + magicIn.getU() + ", I = " + magicIn.getCurrent());
        System.out.println("magicOut: V = " + magicOut.getU() + ", I = " + magicOut.getCurrent());
        System.out.println("s3: V = " + s3.state);
        System.out.println("r2: Vd = " + r2.getU() + ", I = " + r2.getI());
        System.out.println("s4: V = " + s4.state);
        System.out.println("e2: V = " + e2.getU() + ", I = " + e2.getI());
    }

    public boolean containe(State state) {
        return states.contains(state);
    }

    public void setX(State s, double value) {
        s.state = value;
    }

    public double getX(State s) {
        return s.state;
    }

    public double getXSafe(State bPin) {
        return bPin == null ? 0 : getX(bPin);
    }

    public boolean breakSystem() {
        if (breaked) return false;
        while (!breakDestructor.isEmpty()) {
            breakDestructor.pop().destruct();
        }

        for (Component c : component) {
            c.quitSubSystem();
        }
        for (State s : states) {
            s.quitSubSystem();
        }

        if (root != null) {
            for (Component c : component) {
                c.returnToRootSystem(root);
            }
            for (State s : states) {
                s.returnToRootSystem(root);
            }
        }
        root.systems.remove(this);

        invalidate();

        breaked = true;
        return true;
    }

    public void addProcess(ISubSystemProcessFlush p) {
        processF.add(p);
    }

    public void removeProcess(ISubSystemProcessFlush p) {
        processF.remove(p);
    }

    public double getDt() {
        return dt;
    }

    static public class Th {
        public double R, U;

        public boolean isHighImpedance() {
            return R > 1e8;
        }

        @Override
        public String toString() {
            return "Th{" +
                "R=" + R +
                ", U=" + U +
                '}';
        }
    }

    public Th getTh(State d, VoltageSource voltageSource) {
        Th th = new Th();
        double originalU = d.state;

        double otherU = originalU + 5;
        voltageSource.setU(otherU);
        double otherI = solve(voltageSource.getCurrentState());

        voltageSource.setU(originalU);
        double originalI = solve(voltageSource.getCurrentState());

        double Rth = (otherU - originalU) / (originalI - otherI);
        double Uth;
        //if(Double.isInfinite(d.Rth)) d.Rth = Double.MAX_VALUE;
        if (Rth > 10000000000000000000.0 || Rth < 0) {
            Uth = 0;
            Rth = 10000000000000000000.0;
        } else {
            Uth = otherU + Rth * otherI;
        }
        voltageSource.setU(originalU);

        th.R = Rth;
        th.U = Uth;

        if(Double.isNaN(th.U)) {
            th.U = originalU;
            th.R = MnaConst.highImpedance;
        }
        if (Double.isNaN(th.R)) {
            th.U = originalU;
            th.R = MnaConst.highImpedance;
        }

        return th;
    }

    public Th getTh(State d, CurrentSource currentSource) {
        Th th = new Th();
        double originalI = currentSource.getCurrent();

        currentSource.setCurrent(1);
        double n2p = d.getSubSystem().solve(d);

        currentSource.setCurrent(-1);
        double n2m = d.getSubSystem().solve(d);

        th.U = (n2p + n2m) / 2;
        th.R = n2p - th.U;

        if(Double.isNaN(th.U)) {
            th.U = 0;
            th.R = MnaConst.highImpedance;
        }
        if (Double.isNaN(th.R)) {
            th.U = 0;
            th.R = MnaConst.highImpedance;
        }

        currentSource.setCurrent(originalI);
        return th;
    }

    public String toString() {
        String str = "";
        for (Component c: component) {
            if (c != null)
                str += c.toString();
        }
        //str = component.size() + "components";
        return str;
    }

    public int componentSize() {
        return component.size();
    }
}
